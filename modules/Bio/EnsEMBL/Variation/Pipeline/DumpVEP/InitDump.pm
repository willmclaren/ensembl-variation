=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::InitDump;

use strict;
use warnings;

use File::Path qw(make_path);

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

use DBI;

my $DEBUG = 0;

sub param_defaults {
  return {
    'refseq' => 0,
    'merged' => 0,
    'convert' => 0,
    'include_pattern' => '',
    'exclude_pattern' => '',
  };
}

sub fetch_input {
  my $self = shift;

  my $servers = $self->required_param('dump_servers');
  
  my @jobs;
  
  foreach my $server(@$servers) {
    push @jobs, @{$self->get_all_jobs_by_server($server)};
  }

  # this will contain the join jobs
  my $pre_joins = {};

  # make some lists
  foreach my $type(qw(core otherfeatures variation regulation)) {
    my @type_jobs = grep {$_->{type} eq $type} @jobs;
    $self->param($type, \@type_jobs);

    foreach my $job(@type_jobs) {
      my $pre_join = $pre_joins->{$job->{species}}->{$job->{assembly}} ||= {};
      $pre_join->{$job->{type}} = 1;
      $pre_join->{dir_suffix} = $job->{dir_suffix} || '';
      $pre_join->{$_} = $job->{$_} for qw(host user pass port dbname is_multispecies species_id);
    }
  }

  # now create the actual join jobs
  my @join_jobs;

  foreach my $species(keys %$pre_joins) {
    foreach my $assembly(keys %{$pre_joins->{$species}}) {
      my $pre_join = $pre_joins->{$species}->{$assembly};

      my %base_job = (
        species    => $species,
        assembly   => $assembly,
      );
      $base_job{$_} = $pre_join->{$_} for qw(dir_suffix host user pass port dbname is_multispecies species_id);
      
      map {$base_job{$_} = 1} grep {$pre_joins->{$species}->{$assembly}->{$_}} qw(variation regulation);

      # core job
      if($pre_joins->{$species}->{$assembly}->{core}) {
        my %job = %base_job;
        $job{type} = 'core';
        push @join_jobs, \%job;
      }

      # refseq job
      if($pre_joins->{$species}->{$assembly}->{otherfeatures}) {
        my %job = %base_job;
        $job{type} = 'refseq';
        push @join_jobs, \%job;

        # merge job
        if($self->param('merged')) {
          my %job = %base_job;
          $job{type} = 'merged';
          push @join_jobs, \%job;
        }
      }
    }
  }

  $self->param('merges', [grep {$_->{type} eq 'merged'} @join_jobs]);

  $self->param('joins', \@join_jobs);

  return;
}

sub write_output {
  my $self = shift;

  # distribute
  $self->dataflow_output_id({}, 1);

  $self->dataflow_output_id($self->param('core'), 2);
  $self->dataflow_output_id($self->param('otherfeatures'), 3);
  $self->dataflow_output_id($self->param('variation'), 4);
  $self->dataflow_output_id($self->param('regulation'), 5);

  # merge ens+refseq
  $self->dataflow_output_id($self->param('merges'), 6);

  # join
  $self->dataflow_output_id($self->param('joins'), 7);

  # qc
  $self->dataflow_output_id($self->param('joins'), 8);

  # finish (rm dirs)
  $self->dataflow_output_id($self->param('joins'), 9);
  
  return;
}

sub get_all_jobs_by_server {
  my $self = shift;
  my $server = shift;

  my $connection_string = sprintf(
    "DBI:mysql(RaiseError=>1):host=%s;port=%s",
    $server->{host},
    $server->{port}
  );
  
  # connect to DB
  my $dbc = DBI->connect(
    $connection_string, $server->{user}, $server->{pass}
  );
  
  my $version = $self->param('eg_version') || $self->required_param('ensembl_release');

  my $sth = $dbc->prepare(qq{
    SHOW DATABASES LIKE '%\_core\_$version%'
  });
  $sth->execute();
  
  my $db;
  $sth->bind_columns(\$db);
  
  my @dbs;
  push @dbs, $db while $sth->fetch;
  $sth->finish;
  
  # refseq?
  if($self->param('refseq')) {
    $sth = $dbc->prepare(qq{
      SHOW DATABASES LIKE '%\_otherfeatures\_$version%'
    });
    $sth->execute();
    $sth->bind_columns(\$db);
    
    push @dbs, $db while $sth->fetch;
    $sth->finish;
  }

  # remove master and coreexpression
  @dbs = grep {$_ !~ /master|express/} @dbs;

  # filter on pattern if given
  my $pattern = exists($server->{include_pattern}) ? $server->{include_pattern} : $self->param('include_pattern');
  my $exclude = exists($server->{exclude_pattern}) ? $server->{exclude_pattern} : $self->param('exclude_pattern');
  @dbs = grep {$_ =~ /$pattern/i} @dbs if $pattern;
  @dbs = grep {$_ !~ /$exclude/i} @dbs if $exclude;

  my @return;

  foreach my $current_db_name (@dbs) {

    next if $self->is_strain($dbc, $current_db_name);
    
    # special case otherfeatures
    if($current_db_name =~ /otherfeatures/) {

      # check it has refseq transcripts
      $sth = $dbc->prepare(qq{
        SELECT COUNT(*)
        FROM $current_db_name\.transcript
        WHERE stable_id LIKE 'NM%'
        OR source = 'refseq'
      });
      $sth->execute;
      
      my $count;
      $sth->bind_columns(\$count);
      $sth->fetch;
      $sth->finish();
      next unless $count;

      my $species_ids = $self->get_species_id_hash($dbc, $current_db_name);
      
      foreach my $species_id(keys %$species_ids) {
        my $assembly = $self->get_assembly($dbc, $current_db_name, $species_id);
        next unless $assembly;

        # copy server details
        my %species_hash = %$server;
      
        $species_hash{species} = $species_ids->{$species_id};
        $species_hash{species_id} = $species_id;
        $species_hash{assembly} = $assembly;
        $species_hash{dbname} = $current_db_name;
        $species_hash{is_multispecies} = scalar keys %$species_ids > 1 ? 1 : 0;

        # do we have SIFT or PolyPhen?
        if(my $var_db_name = $self->has_var_db($dbc, $current_db_name)) {
          my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name, $species_id);
          $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
        }

        push @return, @{$self->get_all_jobs_by_species_hash(\%species_hash, undef, undef, 'otherfeatures')};
      }
    }
    
    else {
      my $species_ids = $self->get_species_id_hash($dbc, $current_db_name);
      
      # do we have a variation DB?
      my $var_db_name = $self->has_var_db($dbc, $current_db_name);
      
      # do we have a regulation DB?
      my $reg_db_name = $self->has_reg_build($dbc, $current_db_name);

      my $species_count = 0;
      
      foreach my $species_id(keys %$species_ids) {
        my $assembly = $self->get_assembly($dbc, $current_db_name, $species_id);
        next unless $assembly;
        
        # copy server details
        my %species_hash = %$server;
        
        $species_hash{species} = $species_ids->{$species_id};
        $species_hash{species_id} = $species_id;
        $species_hash{assembly} = $assembly;
        $species_hash{dbname} = $current_db_name;
        $species_hash{is_multispecies} = scalar keys %$species_ids > 1 ? 1 : 0;
        
        # do we have SIFT or PolyPhen?
        if($var_db_name) {
          my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name, $species_id);
          $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
        }

        push @return, @{$self->get_all_jobs_by_species_hash(\%species_hash, $var_db_name, $reg_db_name)};

        $species_count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $species_count;
    }
  }
  
  return \@return;
}

sub is_strain {
  my ($self, $dbc, $current_db_name) = @_;

  my $sth = $dbc->prepare("select meta_value from ".$current_db_name.".meta where meta_key = 'species.strain';");
  $sth->execute();
  my $strain_value;
  $sth->bind_columns(\$strain_value);
  $sth->execute();
  $sth->fetch();
  $sth->finish();

  return 1 if $strain_value && $strain_value !~ /^reference/;

  return $current_db_name =~ /mus_musculus_.+?_(core|otherfeatures)/;
}

sub get_species_id_hash {
  my ($self, $dbc, $current_db_name) = @_;

  # get species names by id
  my $sth = $dbc->prepare("select species_id, meta_value from ".$current_db_name.".meta where meta_key = 'species.production_name';");
  $sth->execute();
  
  my ($species_id, $value, $species_ids);
  $sth->bind_columns(\$species_id, \$value);
  
  $species_ids->{$species_id} = $value while $sth->fetch();
  $sth->finish();

  return $species_ids;
}

sub get_assembly {
  my ($self, $dbc, $current_db_name, $species_id) = @_;

  my $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
  $sth->execute();
  my $assembly;
  $sth->bind_columns(\$assembly);
  $sth->execute();
  $sth->fetch();
  $sth->finish();

  return $assembly;
}

sub has_var_db {
  my $self = shift;
  my $dbc = shift;
  my $current_db_name = shift;

  my $var_db_name = $current_db_name;
  $var_db_name =~ s/core|otherfeatures/variation/;
  my $has_var_db;
  my $sth = $dbc->prepare("SHOW DATABASES LIKE '$var_db_name';");
  $sth->execute();
  $sth->bind_columns(\$has_var_db);
  $sth->fetch;
  $sth->finish;

  return $has_var_db ? $var_db_name : undef;
}

sub has_sift_poly {
  my $self = shift;
  my $dbc = shift;
  my $var_db_name = shift;
  my $species_id = shift;
  $species_id ||= 0;

  my $sth = $dbc->prepare(qq{
    SELECT meta_key, meta_value
    FROM $var_db_name\.meta
    WHERE meta_key in ('sift_version','polyphen_version')
    AND (species_id = $species_id OR species_id IS NULL)
  });
  $sth->execute();
  
  my ($key, $val);
  $sth->bind_columns(\$key, \$val);
  
  my %data;

  while($sth->fetch) {
    $key =~ s/\_version//;
    $data{$key} = 'b';
  }
  $sth->finish();

  return \%data;
}

sub has_reg_build {
  my $self = shift;
  my $dbc = shift;
  my $current_db_name = shift;

  my $reg_db_name = $current_db_name;
  $reg_db_name =~ s/core|otherfeatures/funcgen/;
  my $has_reg_db;
  my $sth = $dbc->prepare("SHOW DATABASES LIKE '$reg_db_name';");
  $sth->execute();
  $sth->bind_columns(\$has_reg_db);
  $sth->fetch;
  $sth->finish;
  my $has_reg_build;

  if($has_reg_db) {
    $sth = $dbc->prepare("SELECT version FROM $reg_db_name.regulatory_build");
    $sth->execute();
    $sth->bind_columns(\$has_reg_build);
    $sth->fetch;
    $sth->finish;
  }

  return $has_reg_build ? $reg_db_name : undef;
}

sub get_all_jobs_by_species_hash {
  my $self = shift;
  my $species_hash = shift;
  my $has_var_db = shift;
  my $has_reg_db = shift;
  my $group = shift || 'core';

  # clearing the registry prevents a warning when we connect to
  # mutiple core DBs of the same species (e.g. core, otherfeatures)
  Bio::EnsEMBL::Registry->clear();

  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -group   => 'core',
    -species => $species_hash->{species},
    -port    => $species_hash->{port},
    -host    => $species_hash->{host},
    -user    => $species_hash->{user},
    -pass    => $species_hash->{pass},
    -dbname  => $species_hash->{dbname},
    -MULTISPECIES_DB => $species_hash->{is_multispecies},
    -species_id => $species_hash->{species_id},
  );

  if($self->param('eg')) {
    my $mca = $dba->get_MetaContainerAdaptor;

    if($mca->is_multispecies == 1) {
      my $collection_db = $1 if($mca->dbc->dbname()=~/(.+)\_core/);
      $species_hash->{dir_suffix} = "/".$collection_db;
    }

    $species_hash->{assembly}   = $mca->single_value_by_key('assembly.default');
    $species_hash->{db_version} = $mca->schema_version();
  }
  else {
    $species_hash->{db_version} = $self->param('ensembl_release');
  }

  # get slices
  my $sa = $dba->get_SliceAdaptor;
  my @slices = @{$sa->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$sa->fetch_all('lrg', undef, 1, undef, 1)} if $self->param('lrg');

  # remove/sort out duplicates, in human you get 3 Y slices
  my %by_name;
  $by_name{$_->seq_region_name}++ for @slices;
  if(my @to_fix = grep {$by_name{$_} > 1} keys %by_name) {

    foreach my $name(@to_fix) {

      # remove those with duplicate name
      @slices = grep {$_->seq_region_name ne $name} @slices;

      # add a standard-fetched slice
      push @slices, $sa->fetch_by_region(undef, $name);
    }
  }

  # dumps synonyms
  make_path($self->param('pipeline_dir').'/synonyms');

  open SYN, sprintf(
    ">%s/synonyms/%s_%s_chr_synonyms.txt",
    $self->param('pipeline_dir'),
    $species_hash->{species},
    $species_hash->{assembly}
  ) or die "ERROR: Could not write to synonyms file\n";

  foreach my $slice(@slices) {
    print SYN $slice->seq_region_name."\t".$_->name."\n" for @{$slice->get_all_synonyms};
    delete($slice->{synonym});
  }
  close SYN;

  # remove slices with no transcripts or variants on them
  # otherwise the dumper will create a load of "empty" cache files
  # in species with lots of unplaced scaffolds this means we create
  # masses of pointless directories which take ages to process
  my $ta = $dba->get_TranscriptAdaptor();
  @slices = grep {$ta->count_all_by_Slice($_) || $self->count_vars($has_var_db, $species_hash, $_)} @slices;
  delete($self->{_var_dba}) if $self->{_var_dba};

  # now distribute the slices into jobs
  # jobs can contain multiple slices
  my @jobs;
  my $min_length = 10e6;
  my $added_length = 0;
  my %hash;

  foreach my $slice(sort {$b->end - $b->start <=> $a->end - $a->start} @slices) {
    unless(%hash) {
      %hash = %$species_hash;
      $hash{type} = $group;
      $hash{added_length} = 0;
    }

    push @{$hash{regions}}, {
      chr => $slice->seq_region_name,
      seq_region_id => $slice->get_seq_region_id,
      start => $slice->start,
      end => $slice->end,
    };

    $hash{added_length} += ($slice->end - $slice->start);

    if($hash{added_length} > $min_length) {
      $self->add_to_jobs(\@jobs, \%hash, $has_var_db, $has_reg_db);
      $added_length = 0;
      %hash = ();
    }
  }

  $self->add_to_jobs(\@jobs, \%hash, $has_var_db, $has_reg_db);

  # sort by type then length
  @jobs = sort {$a->{type} cmp $b->{type} || $b->{added_length} <=> $b->{added_length}} @jobs;

  return \@jobs;
}

sub count_vars {
  my $self = shift;
  my $var_db_name = shift;
  my $species_hash = shift;
  my $slice = shift;

  return 0 unless $var_db_name;

  if(!exists($self->{_var_dba})) {
    $self->{_var_dba} = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
      -group   => 'variation',
      -species => $species_hash->{species},
      -port    => $species_hash->{port},
      -host    => $species_hash->{host},
      -user    => $species_hash->{user},
      -pass    => $species_hash->{pass},
      -dbname  => $var_db_name,
      -MULTISPECIES_DB => $species_hash->{is_multispecies},
      -species_id => $species_hash->{species_id},
    );
  }

  my $sth = $self->{_var_dba}->dbc->prepare("SELECT COUNT(*) FROM variation_feature WHERE seq_region_id = ?");
  $sth->execute($slice->get_seq_region_id);

  my $v_count;
  $sth->bind_columns(\$v_count);
  $sth->fetch;
  $sth->finish;

  return $v_count;
}

sub add_to_jobs {
  my $self = shift;
  my $jobs = shift;
  my $hash = shift;
  my $has_var_db = shift;
  my $has_reg_db = shift;

  my %copy = %$hash;
  push @$jobs, \%copy;

  if($has_var_db) {
    my %var = %{$hash};
    $var{type} = 'variation';
    push @$jobs, \%var;
  }

  if($has_reg_db) {
    my %reg = %{$hash};
    $reg{type} = 'regulation';
    push @$jobs, \%reg;
  }
}

1;

