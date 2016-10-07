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

  # make some lists
  foreach my $type(qw(core otherfeatures variation regulation)) {
    $self->param($type, [grep {$_->{type} eq $type} @jobs]);
  }

  return;
}

sub write_output {
  my $self = shift;

  # 1 = distribute dumps (not set here)
  # 2 = normal
  $self->dataflow_output_id($self->param('core'), 2);
  $self->dataflow_output_id($self->param('otherfeatures'), 3);
  $self->dataflow_output_id($self->param('variation'), 4);
  $self->dataflow_output_id($self->param('regulation'), 5);
  
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
  
  my $version = $self->required_param('ensembl_release');

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
    
    # special case otherfeatures
    # check it has refseq transcripts
    if($current_db_name =~ /otherfeatures/) {
    
      # get assembly name
      $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system ORDER BY rank LIMIT 1;");
      $sth->execute();
      my $assembly = $sth->fetchall_arrayref()->[0]->[0];
      die("ERROR: Could not get assembly name from meta table for $current_db_name\n") unless $assembly;
      
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
      
      if($count) {
        
        my $species = $current_db_name;
        $species =~ s/^([a-z]+\_[a-z,1-9]+)(\_[a-z]+)?(.+)/$1$2/;
        $species =~ s/\_otherfeatures$//;
        
        # copy server details
        my %species_hash = %$server;
      
        $species_hash{species} = $species;
        $species_hash{assembly} = $assembly;
        $species_hash{dbname} = $current_db_name;

        # do we have SIFT or PolyPhen?
        if(my $var_db_name = $self->has_var_db($dbc, $current_db_name)) {
          my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name);
          $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
        }

        push @return, @{$self->get_all_jobs_by_species_hash(\%species_hash, undef, undef, 'otherfeatures')};
      }
    }
    
    else {
      # get assembly and species names
      $sth = $dbc->prepare("select species_id, meta_value from ".$current_db_name.".meta where meta_key = 'species.production_name';");
      $sth->execute();
      
      my ($species_id, $value, $species_ids);
      $sth->bind_columns(\$species_id, \$value);
      
      $species_ids->{$species_id} = $value while $sth->fetch();
      $sth->finish();
      
      my $count = 0;
      
      # do we have a variation DB?
      my $var_db_name = $self->has_var_db($dbc, $current_db_name);
      
      # do we have a regulation DB?
      my $reg_db_name = $self->has_reg_build($dbc, $current_db_name);
      
      foreach $species_id(keys %$species_ids) {
        $sth = $dbc->prepare("SELECT version FROM ".$current_db_name.".coord_system WHERE species_id = ".$species_id." ORDER BY rank LIMIT 1;");
        $sth->execute();
        my $assembly;
        $sth->bind_columns(\$assembly);
        $sth->execute();
        $sth->fetch();
        $sth->finish();
        
        next unless $assembly;
        
        # copy server details
        my %species_hash = %$server;
        
        $species_hash{species} = $species_ids->{$species_id};
        $species_hash{assembly} = $assembly;
        $species_hash{dbname} = $current_db_name;
        
        # do we have SIFT or PolyPhen?
        if($var_db_name) {
          my $has_sift_poly = $self->has_sift_poly($dbc, $var_db_name, $species_id);
          $species_hash{$_} = $has_sift_poly->{$_} for keys %$has_sift_poly;
        }

        push @return, @{$self->get_all_jobs_by_species_hash(\%species_hash, $var_db_name, $reg_db_name)};
        $count++;
      }
      
      die("ERROR: Problem getting species and assembly names from $current_db_name; check coord_system table\n") unless $count;
    }
  }
  
  return \@return;
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

  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -group   => 'core',
    -species => $species_hash->{species},
    -port    => $species_hash->{port},
    -host    => $species_hash->{host},
    -user    => $species_hash->{user},
    -pass    => $species_hash->{pass},
    -dbname  => $species_hash->{dbname},
  );

  my $sa = $dba->get_SliceAdaptor;

  my @slices = @{$sa->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$sa->fetch_all('lrg', undef, 1, undef, 1)} if $self->param('lrg');

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

