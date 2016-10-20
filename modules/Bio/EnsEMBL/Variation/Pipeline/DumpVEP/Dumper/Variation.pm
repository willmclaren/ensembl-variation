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
package Bio::EnsEMBL::Variation::Pipeline::DumpVEP::Dumper::Variation;

use strict;
use warnings;

use FileHandle;

use Bio::EnsEMBL::VEP::Config;
use Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation;
use Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation;

our $CAN_USE_TABIX_PM;

BEGIN {
  if (eval { require Bio::EnsEMBL::IO::Parser::VCF4Tabix; 1 }) {
    $CAN_USE_TABIX_PM = 1;
  }
  else {
    $CAN_USE_TABIX_PM = 0;
  }
}

use base qw(Bio::EnsEMBL::Variation::Pipeline::DumpVEP::Dumper);

sub run {
  my $self = shift;

  my $vep_params = $self->get_vep_params();

  # make sure to include failed variants!
  $vep_params->{failed} = 1;

  my $config = Bio::EnsEMBL::VEP::Config->new($vep_params);

  my $region_size = $self->param('region_size');

  my $as = Bio::EnsEMBL::VEP::AnnotationSource::Database::Variation->new({
    config => $config,
    cache_region_size => $region_size,
  });

  my $cache = Bio::EnsEMBL::VEP::AnnotationSource::Cache::Variation->new({
    config => $config,
    cache_region_size => $region_size,
    dir => $self->get_cache_dir($vep_params)
  });

  # create prefixed names here to use later
  if($vep_params->{freq_vcf} && !$self->{freq_vcf}) {
    foreach my $vcf_conf(@{$vep_params->{freq_vcf}}) {
      if(my $prefix = $vcf_conf->{prefix}) {
        $vcf_conf->{prefixed_pops} = [
          map {$_ ? $prefix.'_'.$_ : $prefix}
          @{$vcf_conf->{pops}}
        ];
      }
    }
    $self->{freq_vcf} = $vep_params->{freq_vcf};
  }

  # precache pubmed data
  $self->pubmed($as);

  $self->dump_chrs($as, $cache);

  # bgzip and tabix-index all_vars files
  if($self->param('convert')) {
    my $root = $self->get_cache_dir($vep_params);

    # need to find which column is start
    my @cols = @{$as->get_cache_columns()};
    my %indexes = map {$cols[$_] => $_} 0..$#cols;
    die("Cannot determine start column for indexing\n") unless defined($indexes{start});

    # add 2, 1 to correct for 0->1 indexing, 1 because we've added a chr column
    my $start_i = $indexes{start} + 2;

    foreach my $chr(keys %{{map {$_->{chr} => 1} @{$self->param('regions')}}}) {

      next unless -e "$root/$chr/all_vars";

      $self->run_cmd(
        sprintf(
          "sort -k%i,%in %s/%s/all_vars | bgzip -c > %s/%s/all_vars.gz",
          $start_i, $start_i, $root, $chr, $root, $chr
        )
      );

      $self->run_cmd("tabix -s 1 -b $start_i -e $start_i $root/$chr/all_vars.gz");

      unlink("$root/$chr/all_vars");
    }
  }

  $self->dump_info($as, $self->get_cache_dir($vep_params));
  
  return;
}

sub dump_info {
  my ($self, $as, $dir) = @_;

  my $info_file = $dir.'/info.txt_variation';
  return if -e $info_file;

  open OUT, ">$info_file";

  # var cache cols
  print OUT
    "variation_cols\t".
    join(",",
      @{$as->get_cache_columns()},
      'pubmed',
      map {@{$_->{prefixed_pops} || $_->{pops}}} @{$self->{freq_vcf} || []}
    )."\n";

  my $info = $as->info;
  print OUT "source_$_\t".$info->{$_}."\n" for keys %$info;

  close OUT;

  $self->dump_info_converted($as, $dir) if $self->param('convert');
}

sub dump_info_converted {
  my ($self, $as, $dir) = @_;

  my $info_file = $dir.'/info.txt_variation_converted';
  return if -e $info_file;

  open OUT, ">$info_file";

  # var cache cols
  print OUT
    "variation_cols\t".
    join(",",
      'chr',
      @{$as->get_cache_columns()},
      'pubmed',
      map {@{$_->{prefixed_pops} || $_->{pops}}} @{$self->{freq_vcf} || []}
    )."\n";

  print OUT "var_type\ttabix\n";

  my $info = $as->info;
  print OUT "source_$_\t".$info->{$_}."\n" for keys %$info;

  close OUT;
}

sub get_dumpable_object {
  my ($self, $as, $sr, $chr, $s) = @_;
  return $as->get_features_by_regions_uncached([[$sr, $s]], 1);
}

sub dump_obj {
  my ($self, $obj, $file, $chr) = @_;

  open DUMP, "| gzip -9 -c > ".$file or die "ERROR: Could not write to dump file $file\n";

  my $all_vars_fh;
  if($self->param('convert')) {
    my $all_vars_file = $file;
    $all_vars_file =~ s/[^\/]+?$//;
    $all_vars_file .= 'all_vars';

    $all_vars_fh = FileHandle->new();
    $all_vars_fh->open(">>$all_vars_file") or die $!;
  }

  # get freqs from VCFs?
  $self->freqs_from_vcf($obj, $chr) if $self->{freq_vcf};

  my $pubmed = $self->pubmed;

  foreach my $v(@$obj) {
    my @tmp = (
      $v->{variation_name},
      $v->{failed} == 0 ? '' : $v->{failed},
      $v->{somatic} == 0 ? '' : $v->{somatic},
      $v->{start},
      $v->{end} == $v->{start} ? '' : $v->{end},
      $v->{allele_string},
      $v->{strand} == 1 ? '' : $v->{strand},
      $v->{minor_allele} || '',
      defined($v->{minor_allele_freq}) && $v->{minor_allele_freq} =~ /^[0-9\.]+$/ ? sprintf("%.4f", $v->{minor_allele_freq}) : '',
      $v->{clin_sig} || '',
      $v->{phenotype_or_disease} == 0 ? '' : $v->{phenotype_or_disease},
    );

    push @tmp, $pubmed->{$v->{variation_name}} || '';

    if($self->{freq_vcf}) {
      foreach my $pop(map {@{$_->{prefixed_pops} || $_->{pops}}} @{$self->{freq_vcf}}) {
        push @tmp, $v->{$pop} || '';
      }
    }

    print DUMP join(" ", @tmp);
    print DUMP "\n";

    if($all_vars_fh) {
      print $all_vars_fh join("\t", $chr, map {$_ eq '' ? '.' : $_} @tmp);
      print $all_vars_fh "\n";
    }
  }

  close DUMP;
}

sub freqs_from_vcf {
  my $self = shift;
  my $obj = shift;
  my $chr = shift;

  die("ERROR: Cannot use freqs_from_vcf without Bio::EnsEMBL::IO::Parser::VCF4Tabix module\n") unless $CAN_USE_TABIX_PM;

  # sort by pos and exclude somatic vars
  my @list =
    grep {!$_->{somatic}}
    sort {$a->{start} <=> $b->{start} || $a->{end} <=> $b->{end}}
    @$obj;
  return unless scalar @list;

  # put into a hash so we can lookup by pos
  my %by_pos;
  for(@list) {
    # add both start and start - 1 to match indels with preceding base attached
    push @{$by_pos{$_->{start}}}, $_;
    push @{$by_pos{$_->{start} - 1}}, $_;  
  }
  
  # iterate over each VCF file in the config
  foreach my $vcf_conf(@{$self->{freq_vcf}}) {
    my $file = $vcf_conf->{file};
    next unless -e $file;

    my $prefix = $vcf_conf->{prefix} || '';
    $prefix .= '_' if $prefix && $prefix !~ /\_$/;

    my $parser = $self->{_vcf_parsers}->{$file} ||= Bio::EnsEMBL::IO::Parser::VCF4Tabix->open($file);
    next unless $parser;

    $parser->seek($chr, $list[0]->{start} - 1, $list[-1]->{end} + 1);

    while($parser->next) {
      my $match = 0;

      my $ref = $parser->get_reference;
      my $alts = $parser->get_alternatives;
      my $start = $parser->get_raw_start;

      next unless $by_pos{$start};

      # find out if all the alts start with the same base
      # ignore "*"-types
      my %first_bases = map {substr($_, 0, 1) => 1} grep {!/\*/} ($ref, @$alts);

      if(scalar keys %first_bases == 1) {
        $ref = substr($ref, 1) || '-';
        $start++;

        my @new_alts;

        foreach my $alt_allele(@$alts) {
          $alt_allele = substr($alt_allele, 1) unless $alt_allele =~ /\*/;
          $alt_allele = '-' if $alt_allele eq '';
          push @new_alts, $alt_allele;
        }

        $alts = \@new_alts;
      }

      foreach my $v(@{$by_pos{$start}}) {

        if(grep {$v->{variation_name} eq $_} @{$parser->get_IDs || []}) {
          $match = 1;
        }

        # allele string and coords match
        elsif($v->{start} == $start) {
          if($v->{allele_string} eq join('/', $ref, @$alts)) {
            $match = 1;
          }
          else {
            my @v_alleles = split('/', $v->{allele_string});
            my @vcf_alleles = ($ref, @$alts);

            # ref allele must match
            if($v_alleles[0] eq $vcf_alleles[0]) {

              # check if the VCF alleles exist in the allele string
              my %h = map {$_ => 1} @v_alleles;
              my @m = grep {$h{$_}} @vcf_alleles;
              if(scalar @m == scalar @vcf_alleles) {
                $match = 1;
              }
            }
          }
        }

        if($match) {

          # get the alleles from the VCF entry
          my @vcf_alleles = @$alts;

          # get INFO field data from VCF
          my $info = $parser->get_info;

          # have to process ExAC differently from 1KG and ESP
          if($prefix =~ /exac/i) {
            foreach my $pop(@{$vcf_conf->{pops}}) {
              my $suffix = $pop ? '_'.$pop : '';

              my $store_name = $prefix.$pop;
              $store_name =~ s/\_$//;

              if(my $f = $info->{'AF'.$suffix}) {
                my @split = split(',', $f);
                $v->{$store_name} = join(',',
                  map {$alts->[$_].':'.$split[$_]}
                  0..$#split
                );
              }
              elsif(my ($c, $n) = ($info->{'AC'.$suffix}, $info->{'AN'.$suffix})) {
                
                # safety check against div by zero
                next unless $n;

                my @split = split(',', $c);
                $v->{$store_name} = join(',',
                  map {$alts->[$_].':'.sprintf('%.4g', $split[$_] / $n)}
                  0..$#split
                );
              }
            }
          }

          else {
            foreach my $pop(@{$vcf_conf->{pops}}) {
              my $tmp_prefix = $pop ? $pop.'_' : '';

              my $store_name = $prefix.$pop;
              $store_name =~ s/\_$//;

              if(my $f = $info->{$tmp_prefix.'AF'}) {
                my @split = split(',', $f);
                $v->{$store_name} = join(',',
                  map {$alts->[$_].':'.$split[$_]}
                  0..$#split
                );
              }
              elsif(my $c = $info->{$tmp_prefix.'AC'}) {
                
                my @split = split(',', $c);
                my $n;
                $n += $_ for @split;
                next unless $n;

                # ESP VCFs include REF as last allele, just to annoy everyone
                pop @split;

                $v->{$store_name} = join(',',
                  map {$alts->[$_].':'.sprintf('%.4g', $split[$_] / $n)}
                  0..$#split
                );
              }
            }
          }
        }
      }
    }
  }
}

# we want pubmed IDs by rsID
# this routine gets them from the DB in the first instance
# then caches this fetch to disk so future processes can read it
# as the query takes a little while to run
sub pubmed {
  my $self = shift;
  my $as = shift;

  if(!exists($self->{_pubmed})) {

    my %pm;
    my $file = sprintf(
      '%s/pubmed_%s_%s_%s.txt',
      $self->required_param('pipeline_dir'),
      $self->required_param('species'),
      $self->required_param('ensembl_release'),
      $self->required_param('assembly')
    );
    my $lock = $file.'.lock';

    my $sleep_count = 0;
    if(-e $lock) {
      while(-e $lock) {
        sleep 1;
        die("I've been waiting for $lock to be removed for $sleep_count seconds, something may have gone wrong\n") if ++$sleep_count > 900;
      }
    }

    if(-e $file) {
      open IN, $file;
      while(<IN>) {
        chomp;
        my @split = split;
        $pm{$split[0]} = $split[1];
      }
      close IN;
    }

    elsif($as) {
      open OUT, ">$lock";
      print OUT "$$\n";
      close OUT;
      $self->{_lock_file} = $lock;

      my $sth = $as->get_adaptor('variation', 'variation')->dbc->prepare(qq{
        SELECT v.name, GROUP_CONCAT(p.pmid)
        FROM variation v, variation_citation c, publication p
        WHERE v.variation_id = c.variation_id
        AND c.publication_id = p.publication_id
        AND p.pmid IS NOT NULL
        GROUP BY v.variation_id
      });
      $sth->execute;

      my ($v, $p);
      $sth->bind_columns(\$v, \$p);

      open OUT, ">$file";

      while($sth->fetch()) {
        $pm{$v} = $p;
        print OUT "$v\t$p\n";
      }

      close OUT;

      unlink($lock);

      $sth->finish();
    }

    $self->{_pubmed} = \%pm;
  }

  return $self->{_pubmed};
}

sub DESTROY {
  unlink($_[0]->{_lock_file}) if $_[0]->{_lock_file};
}

1;
