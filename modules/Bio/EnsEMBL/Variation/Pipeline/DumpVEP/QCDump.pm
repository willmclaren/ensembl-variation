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

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::CacheDir;

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
}

sub qc {
  my ($self, $mod) = @_;

  my $type      = $self->param('type');
  my $has_var   = $self->param('variation');
  my $has_reg   = $self->param('regulation');
  my $converted = $mod && $mod =~ /tabix/;
  my $species   = $self->required_param('species');
  my $assembly  = $self->required_param('assembly');

  my $tar_file = $self->get_tar_file_name(
    $dump_dir,
    $species,
    $assembly,
    $mod
  );

  die("ERROR: Tar file $tar_file not found\n") unless -e $tar_file;

  my $dump_dir = $self->dump_dir();
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
  my $config_obj = Bio::EnsEMBL::VEP::Config->new({dir => $extracted_dir, offline => 1, species => $species});
  my $cache_dir_obj = Bio::EnsEMBL::VEP::CacheDir->new({dir => $extracted_dir, config => $config_obj});

  # these subs check different aspects of the cache
  $self->check_info($cache_dir_obj, $converted);

  $self->check_annotation_sources($cache_dir_obj, $converted);
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
      unless $info->{variation_cols} =~ /(\w+\,?)+/;

    if($converted) {
      die("ERROR: variation_cols value does not include chr\n") unless $info->{variation_cols} =~ /^chr\,/;
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

1;
