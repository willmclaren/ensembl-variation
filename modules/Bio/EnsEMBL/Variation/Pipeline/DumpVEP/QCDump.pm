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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::QCDump;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::BaseVEP);

use File::Path qw(make_path rmtree);
use Digest::MD5 qw(md5_hex);
use JSON;
use FileHandle;
use Data::Dumper;
$Data::Dumper::Indent = 1;

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::CacheDir;
use Bio::EnsEMBL::VEP::Runner;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;

sub param_defaults {
  return {
    'variation'  => 0,
    'regulation' => 0,
    'dir_suffix' => '',
    'convert'    => 0,
  };
}

sub run {
  my $self = shift;

  $self->qc();
  $self->qc('_tabixconverted') if $self->param('convert');

  return 1;
}

sub qc {
  my ($self, $mod) = @_;

  my $type      = $self->param('type');
  my $has_var   = $self->param('variation');
  my $has_reg   = $self->param('regulation');
  my $converted = $mod && $mod =~ /tabix/;
  my $species   = $self->required_param('species');
  my $assembly  = $self->required_param('assembly');

  my $dump_dir = $self->dump_dir();
  my $tar_file = $self->get_tar_file_name(
    $dump_dir,
    $species.($type eq 'core' ? '' : "_$type"),
    $assembly,
    $mod
  );

  die("ERROR: Tar file $tar_file not found\n") unless -e $tar_file;

  my $qc_dir = $dump_dir.'/qc/'.md5_hex($self->input_id);
  rmtree($qc_dir) if -d $qc_dir;
  make_path($qc_dir);

  # untar
  $self->run_cmd(sprintf('tar -C %s -xzf %s', $qc_dir, $tar_file));

  # check dir exists
  my $method_name = ($type eq 'core' ? '' : $type.'_').'species_suffix';
  my $extracted_dir = $qc_dir.'/'.$self->$method_name;
  die("ERROR: Expected to find $extracted_dir\n") unless -d $extracted_dir;

  # check info.txt exists
  die("ERROR: Expected to find $extracted_dir/info.txt") unless -e "$extracted_dir/info.txt";

  # create objects
  my $config_obj = Bio::EnsEMBL::VEP::Config->new({
    dir => $extracted_dir,
    offline => 1,
    species => $species,
    assembly => $assembly,
    check_existing => 1,
    regulatory => 1,
  });
  my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({dir => $extracted_dir, config => $config_obj});

  # these subs check different aspects of the cache
  $self->check_info($cache_dir_obj, $converted);

  $self->check_annotation_sources($cache_dir_obj, $converted);

  $self->check_dirs($cache_dir_obj, $converted);

  $self->run_test_set($qc_dir) if $has_var && $type ne 'otherfeatures';

  # clean up
  rmtree($qc_dir);
}

sub check_info {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $species  = $self->required_param('species');
  my $assembly = $self->required_param('assembly');

  my $info = $cache_dir_obj->info;
  die("ERROR: No keys in info\n") unless keys %$info;

  die("ERROR: species info key does not exist\n") unless $info->{species};
  die("ERROR: species info key does not match job species name\n")
    unless $info->{species} eq $species;

  die("ERROR: assembly info key does not exist\n") unless $info->{assembly};
  die("ERROR: assembly info key does not match job assembly name\n")
    unless $info->{assembly} eq $assembly;

  # var_type should be tabix if converted
  if($converted) {
    die("ERROR: var_type info key is not set to tabix\n")
      unless $info->{var_type} && $info->{var_type} eq 'tabix';
  }

  # check version_data
  die("ERROR: version_data key not found\n") unless $info->{version_data} && ref($info->{version_data}) eq 'HASH';

  # check variation_cols defined
  if($self->param('variation')) {
    die("ERROR: variation_cols info key not found\n") unless $info->{variation_cols};

    die("ERROR: variation_cols info value looks wrong: ".$info->{variation_cols}."\n")
      unless ref($info->{variation_cols}) eq 'ARRAY';

    if($converted) {
      die("ERROR: variation_cols value does not include chr\n")
        unless $info->{variation_cols}->[0] eq 'chr';
    }
  }
}

sub check_annotation_sources {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $annotation_sources = $cache_dir_obj->get_all_AnnotationSources;

  die("ERROR: No annotation sources retrieved\n") unless ref($annotation_sources) eq 'ARRAY' && scalar @$annotation_sources;

  # expect core always
  die("ERROR: No Transcript annotation source found\n")
    unless grep {ref($_) eq 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Transcript'} @$annotation_sources;

  if($self->param('variation')) {
    my $t = 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation'.($converted ? 'Tabix' : '');
    die("ERROR: No Variation annotation source found\n")
      unless grep {ref($_) eq $t} @$annotation_sources;
  }
  if($self->param('regulation')) {
    die("ERROR: No Regulation annotation source found\n")
      unless grep {ref($_) eq 'Bio::EnsEMBL::VEP::AnnotationSource::Cache::RegFeat'} @$annotation_sources;
  }
}

sub check_dirs {
  my ($self, $cache_dir_obj, $converted) = @_;

  my $dir = $cache_dir_obj->dir;
  opendir DIR, $dir or die $!;
  my %dirs_in_cache;
  $dirs_in_cache{$_} = 1 for grep {!/^\./ && -d "$dir/$_" } readdir(DIR);
  closedir DIR;

  die("ERROR: No chromosome dirs found\n") unless scalar keys %dirs_in_cache;

  # connect to DB to find all the slices we should have dirs for
  Bio::EnsEMBL::Registry->clear();
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -group   => 'core',
    -species => $self->param('species'),
    -port    => $self->param('port'),
    -host    => $self->param('host'),
    -user    => $self->param('user'),
    -pass    => $self->param('pass'),
    -dbname  => $self->param('dbname'),
    -MULTISPECIES_DB => $self->param('is_multispecies'),
    -species_id => $self->param('species_id'),
  );

  # get slices
  my $sa = $dba->get_SliceAdaptor;
  my @slices = @{$sa->fetch_all('toplevel')};
  push @slices, map {$_->alternate_slice} map {@{$_->get_all_AssemblyExceptionFeatures}} @slices;
  push @slices, @{$sa->fetch_all('lrg', undef, 1, undef, 1)} if $self->param('lrg');

  # check for dirs that don't appear in this list of slices
  delete($dirs_in_cache{$_->seq_region_name}) for @slices;
  die("ERROR: Found ".(scalar keys %dirs_in_cache)." dirs in cache that have no corresponding slice\n") if scalar keys %dirs_in_cache;

  # now filter out slices with no transcripts
  my $ta = $dba->get_TranscriptAdaptor();
  @slices = grep {$ta->count_all_by_Slice($_)} @slices;

  foreach my $slice(@slices) {
    my $chr = $slice->seq_region_name;
    die("ERROR: Dir $dir/$chr missing\n") unless -d "$dir/$chr";

    opendir SUB, "$dir/$chr" or die $!;
    my @files = grep {!/^\./} readdir(SUB);
    closedir SUB;

    die("ERROR: $dir/$chr is empty\n") unless @files;

    foreach my $file(@files) {
      die("ERROR: File name $file does not look right\n")
        unless $file =~ /\d+\-\d+(\_(var|reg))?\.gz/
        or $file eq 'all_vars.gz'
        or $file eq 'all_vars.gz.tbi';

      if($file =~ /var/) {
        die("ERROR: Cache should not contain var files\n") unless $self->param('variation');
        die("ERROR: Converted cache should not contain per-MB files\n") if $converted && $file =~ /\_var\.gz/;
        die("ERROR: Basic cache should not contain all_vars.gz file\n") if !$converted && $file eq 'all_vars.gz';
      }
      elsif($file =~ /reg/) {
        die("ERROR: Cache should not contain var files\n") unless $self->param('regulation');
      }
    }
  }
}

sub run_test_set {
  my $self = shift;
  my $qc_dir = shift;

  # get a test file
  my $test_file = $self->generate_test_file();

  # create a VEP Runner
  my $runner = Bio::EnsEMBL::VEP::Runner->new({
    species => $self->param('species'),
    assembly => $self->param('assembly'),
    cache_version => $self->param('eg_version') || $self->param('ensembl_release'),
    offline => 1,
    database => 0,
    dir => $qc_dir,
    input_file => $test_file,
    format => 'ensembl',
    output_format => 'json',
    safe => 1,
    quiet => 1,
    no_stats => 1,
    regulatory => $self->param('regulation'),
    check_existing => $self->param('variation'),
  });

  my $json = JSON->new();

  while(my $line = $runner->next_output_line) {
    my $data = $json->decode($line);

    die("ERROR: input field not found in JSON output\n".(Dumper $data)) unless $data->{input};
    my @input = split("\t", $data->{input});

    my $feature_type = $input[-1];
    my $expected_cons = join(",", sort split(",", $input[-2]));
    my $feature_id = $input[-3];

    # check consequence type
    die("ERROR: no data for $feature_type found in JSON output\n".(Dumper $data)) unless $data->{$feature_type.'_consequences'};

    if(my ($blob) = grep {$_->{$feature_type.'_id'} eq $feature_id} @{$data->{$feature_type.'_consequences'}}) {
      die("ERROR: no consequence_terms field in blob\n".(Dumper $data)) unless $blob->{consequence_terms};
      my $got_cons = join(",", sort @{$blob->{consequence_terms}});
      die("ERROR: consequence_types don't match, expected: $expected_cons, got: $got_cons\n".(Dumper $data)) unless $expected_cons eq $got_cons;
    }
    else {
      die("ERROR: no data for $feature_id found in JSON output\n".(Dumper $data));
    }

    # check colocated variants
    if($self->param('variation')) {
      my $expected_var_id = $input[-4];
      die("ERROR: no data for colocated_variants found in JSON output\n".(Dumper $data)) unless $data->{colocated_variants};

      die("ERROR: expected var id $expected_var_id not found\n".(Dumper $data))
        unless grep {$_->{id} eq $expected_var_id} @{$data->{colocated_variants}};
    }
  }
}

sub generate_test_file {
  my $self = shift;

  my $file = sprintf(
    '%s/qc/%s_%s_test_input.txt',
    $self->dump_dir(),
    $self->param('species'),
    $self->param('assembly')
  );

  return $file if -e $file;

  # use a lock file in case two processes try writing to the same one
  my $lock = $file.'.lock';

  # wait 3hrs for other process to finish, it will be a long-running query in human
  my $sleep_count = 0;
  if(-e $lock) {
    while(-e $lock) {
      sleep 1;
      die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > (60*60*3);
    }

    die("ERROR: Finished waiting but $file not found\n") unless -e $file;

    return $file;
  }

  open OUT, ">$lock" or die $!;
  print OUT "$$\n";
  close OUT;
  $self->{_lock_file} = $lock;

  my $dbname = $self->param('dbname');
  $dbname =~ s/core|otherfeatures/variation/;

  my $dba = Bio::EnsEMBL::Variation::DBSQL::DBAdaptor->new(
    -group   => 'variation',
    -species => $self->param('species'),
    -port    => $self->param('port'),
    -host    => $self->param('host'),
    -user    => $self->param('user'),
    -pass    => $self->param('pass'),
    -dbname  => $dbname,
    -MULTISPECIES_DB => $self->param('is_multispecies'),
    -species_id => $self->param('species_id'),
  );

  my @tables = ('transcript', ($self->param('regulation') ? qw(regulatory_feature motif_feature) : ()));

  my $fh;

  foreach my $table(@tables) {
    my $sth = $dba->dbc->prepare(qq{
      SELECT s.name, vf.seq_region_start, vf.seq_region_end, vf.allele_string, vf.seq_region_strand,
      vf.variation_name, tv.feature_stable_id, tv.consequence_types
      FROM seq_region s, $table\_variation tv, variation_feature vf
      LEFT JOIN failed_variation fv on vf.variation_id = fv.variation_id
      WHERE s.seq_region_id = vf.seq_region_id
      AND vf.variation_feature_id = tv.variation_feature_id
      AND fv.variation_id IS NULL
      GROUP BY vf.seq_region_id, tv.consequence_types
    }, {mysql_use_result => 1});
    $sth->execute;

    $fh ||= FileHandle->new($file, 'r');

    while(my $arrayref = $sth->fetchrow_arrayref()) {
      print $fh join("\t", @$arrayref, $table)."\n";
    }

    $sth->finish();
  }

  $fh->close();

  unlink($lock);

  return $file;
}

sub DESTROY {
  unlink($_[0]->{_lock_file}) if $_[0]->{_lock_file};
}

1;
