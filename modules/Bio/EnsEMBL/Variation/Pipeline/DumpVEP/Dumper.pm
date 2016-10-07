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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::Dumper;

use strict;
use warnings;

use Storable qw(nstore_fd);
use File::Path qw(mkpath);

use base qw(Bio::EnsEMBL::Variation::Pipeline::BaseVariationProcess);

sub param_defaults {
  return {
    'sift'     => 0,
    'polyphen' => 0,
    'eg'       => 0,
  };
}

sub get_vep_params {
  my $self = shift;

  my $params = {};

  # basic params
  $params->{eg}      = $self->param('eg');
  $params->{debug}   = $self->param('debug');
  $params->{species} = $self->required_param('species');
  $params->{dir}     = $self->required_param('pipeline_dir');

  $params->{is_multispecies} = 0;

  if($params->{eg}){
     my $meta_container = Bio::EnsEMBL::Registry->get_adaptor($params->{species},'core','MetaContainer');

     if($meta_container->is_multispecies()==1){
        my $collection_db=$1 if($meta_container->dbc->dbname()=~/(.+)\_core/);
        $params->{dir} .= "/".$collection_db;
        make_path($params->{dir});

        $params->{is_multispecies} = 1;
     }

     $params->{assembly}   = $meta_container->single_value_by_key('assembly.default');
     $params->{version}    = $meta_container->schema_version();
     $params->{eg_version} = $self->param('eg_version');
     $params->{host}       = $meta_container->dbc->host();
     $params->{port}       = $meta_container->dbc->port();
     $params->{user}       = $meta_container->dbc->username();
     $params->{pass}       = $meta_container->dbc->password() ? $meta_container->dbc->password() : '';
     $meta_container->dbc()->disconnect_if_idle();
     
     $self->param('assembly', $params->{assembly});
     $self->param('ensembl_release', $params->{version});
  }
  else {
     $params->{assembly} = $self->required_param('assembly');
     $params->{version}  = $self->required_param('ensembl_release');
     $params->{host}     = $self->required_param('host');
     $params->{port}     = $self->required_param('port');
     $params->{user}     = $self->required_param('user');
     $params->{pass}     = $self->required_param('pass') ? $self->required_param('pass') : '';
  }

  # sift, polyphen
  $params->{$_} = 'b' for grep {$self->param($_)} qw(sift polyphen);

  # species-specific
  my $species_flags = $self->param('species_flags');
  
  if(my $flags = $species_flags->{$params->{species}}) {
    
    # assembly-specific
    if(my $as = $flags->{assembly_specific}) {
      delete $flags->{assembly_specific};

      my $assembly = $params->{assembly};
      
      if(my $as_flags = $as->{$assembly}) {
        $params->{$_} = $as_flags->{$_} for keys %$as_flags;
      }
    }

    $params->{$_} = $flags->{$_} for keys %$flags;
  }

  return $params;
}

sub get_cache_dir {
  my ($self, $vep_params) = @_;

  return sprintf(
    '%s/%s/%s_%s',
    $vep_params->{dir},
    $vep_params->{species}.($vep_params->{refseq} ? '_refseq' : ''),
    $vep_params->{eg_version} || $vep_params->{version},
    $vep_params->{assembly}
  );
}

sub dump_chrs {
  my ($self, $db_as, $cache_as) = @_;

  my $sa = $db_as->get_adaptor('core', 'slice');

  my @regions = @{$self->param('regions')};

  my $region_size = $self->param('region_size');

  while(my $region = shift @regions) {
    my ($chr, $sr, $slice_start, $slice_end) = (
      $region->{chr},
      $region->{seq_region_id},
      $region->{start},
      $region->{end}
    );

    my $slice = $sa->fetch_by_seq_region_id($sr);

    my $s = int($slice_start / $region_size);
    my $l = int($slice_end / $region_size);
    my $first = 1;

    while($s <= $l) {
      my $obj = $self->get_dumpable_object($db_as, $sr, $chr, $s);
      my $file = $cache_as->get_dump_file_name($chr, ($s  * $region_size) + 1, ($s + 1) * $region_size);

      if($first) {
        my $filedir = $file;
        $filedir =~ s/\/[^\/]+$//;
        mkpath($filedir) unless -d $filedir;
        $first = 0;
      }

      $self->dump_obj($obj, $file, $chr);

      $self->post_dump($obj, $db_as, $chr);

      $db_as->clean_cache();
      $s++;
    }
  }
}

sub dump_obj {
  my $self = shift;
  my $obj = shift;
  my $file = shift;

  open my $fh, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file";
  nstore_fd($obj, $fh);
  close $fh;
}

sub post_dump {
  return 1;
}

1;
