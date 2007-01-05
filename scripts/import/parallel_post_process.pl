use strict;
use warnings;

use Getopt::Long;
use ImportUtils qw(debug);
use Bio::EnsEMBL::Utils::Exception qw(warning throw verbose);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use DBH;
use DBI qw(:sql_types);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use FindBin qw( $Bin );
use Data::Dumper;
use constant MAX_GENOTYPES => 6_000_000; #max number of genotypes per file. When more, split the file into regions
use constant REGIONS => 15; #number of chunks the file will be splited when more than MAX_GENOTYPES

my ($TMP_DIR, $TMP_FILE, $LIMIT);


my ($vhost, $vport, $vdbname, $vuser, $vpass,
    $chost, $cport, $cdbname, $cuser, $cpass,
    $limit, $num_processes, $top_level, $species,
    $variation_feature, $flanking_sequence, $variation_group_feature,
    $transcript_variation, $ld_populations, $reverse_things);

$variation_feature = $flanking_sequence = $variation_group_feature = $transcript_variation = $ld_populations = $reverse_things = '';

GetOptions('tmpdir=s'  => \$ImportUtils::TMP_DIR,
	   'tmpfile=s' => \$ImportUtils::TMP_FILE,
	   'species=s' => \$species,
	   'limit=i'   => \$limit,
	   'top_level=i' => \$top_level,
	   'num_processes=i' => \$num_processes,
	   'variation_feature' => \$variation_feature,
	   'flanking_sequence' => \$flanking_sequence,
	   'variation_group_feature' => \$variation_group_feature,
	   'transcript_variation' => \$transcript_variation,
	   'reverse_things'       => \$reverse_things,
	   'ld_populations' => \$ld_populations );

$num_processes ||= 1;

$LIMIT = ($limit) ? " $limit " : ''; #will refer to position in a slice

usage('-num_processes must at least be 1') if ($num_processes == 0);

warn("Make sure you have a updated ensembl.registry file!\n");

my $registry_file ||= $Bin . "/ensembl.registry";

Bio::EnsEMBL::Registry->load_all( $registry_file );

my $cdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
my $vdba = Bio::EnsEMBL::Registry->get_DBAdaptor($species,'variation');


my $dbVar = $vdba->dbc->db_handle;
my $dbCore = $cdba;

#added default options
$chost = $cdba->dbc->host; 
$cuser    ||= 'ensro';
$cport = $cdba->dbc->port ;
$cdbname = $cdba->dbc->dbname;

$vhost = $vdba->dbc->host;
$vport = $vdba->dbc->port;
$vuser    ||= 'ensadmin';
$vdbname = $vdba->dbc->dbname;
$vpass = $vdba->dbc->password;

$TMP_DIR  = $ImportUtils::TMP_DIR;
$TMP_FILE = $ImportUtils::TMP_FILE;


##Apart from human and mouse, we directly import top_level coordinates from dbSNP
if (! defined $top_level && $species !~ /hum|homo|mouse|mus|mosquito|anoph|dog|can/i) {
  $top_level=1;
}

#we need to create a tmp table to store variations that we filter out
create_failed_variation_table($dbVar);
parallel_variation_feature($dbVar, $top_level) if ($variation_feature);
parallel_flanking_sequence($dbVar) if ($flanking_sequence);
parallel_variation_group_feature($dbVar) if ($variation_group_feature);
parallel_transcript_variation($dbVar) if ($transcript_variation);
parallel_ld_populations($dbVar) if ($ld_populations);
reverse_things($dbVar) if ($reverse_things);

#will take the number of processes, and divide the total number of entries in the variation_feature table by the number of processes
sub parallel_variation_feature{
    my $dbVar = shift;
    my $top_level = shift;

    my $min_variation; #minim variation_feature_id
    my $max_variation; #maximum variation_feature_id
    my $call;
    my $variation_status_file = "status_file_variation_feature_$$\.log";
    #first, create the log file for the variation_feature
    open STATUS, ">$TMP_DIR/$variation_status_file"
	or throw("Could not open tmp file: $TMP_DIR/$variation_status_file\n"); 
    close STATUS;
    #then, calculate the rows for each subprocess
    my $sth = $dbVar->prepare(qq{SELECT min(variation_feature_id),max(variation_feature_id)
			   FROM variation_feature
    });
    $sth->execute();
    ($min_variation, $max_variation) = $sth->fetchrow_array();
    $sth->finish();
    my $dbname = $vdbname; #get the name of the database to create the file    
    
    #to get haplotype seq_region_id
    my $sth1 = $dbCore->dbc->db_handle->prepare(qq{SELECT seq_region_id 
                                                   FROM seq_region_attrib sa, attrib_type at 
                                                   WHERE sa.attrib_type_id=at.attrib_type_id 
                                                   AND at.code='non_ref'});
    $sth1->execute();
    my @hap_seq_region_ids;
    while (my ($seq_region_id) = $sth1->fetchrow_array()) {
      push @hap_seq_region_ids, $seq_region_id;
    }

    my $in_string = "WHERE  seq_region_id not IN (" . join (",",@hap_seq_region_ids) . ")";
    if (scalar @hap_seq_region_ids <1) {
      $in_string = "";
    }
    #create a temporary table to store the map_weight, that will be deleted by the last process
    $dbVar->do(qq{CREATE TABLE tmp_map_weight
                SELECT variation_id, count(*) as count
                FROM   variation_feature
                $in_string
                GROUP BY variation_id}
	       );
    $dbVar->do(qq{ALTER TABLE tmp_map_weight 
		      ADD INDEX variation_idx(variation_id,count)});

    my $sub_variation = int(($max_variation - $min_variation)/ $num_processes);

    #the limit will be (AND variation_feature_id > min and variation_feature_id < max)
    for (my $i = 0; $i < $num_processes ; $i++){
	$limit = "AND variation_feature_id <= " . (($i+1) * $sub_variation + $min_variation-1) . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i+1 < $num_processes);
	$limit =  "AND variation_feature_id <= " .  $max_variation . " AND variation_feature_id >= " . ($i*$sub_variation + $min_variation) if ($i + 1 == $num_processes); #the last one takes the left rows
	$call = "bsub -q normal -J $dbname\_variation_job_$i -m 'bc_hosts' -o $TMP_DIR/output_variation_feature_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_variation_feature.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit '$limit' -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $variation_status_file ";
	$call .= "-cpass $cpass " if ($cpass);
	$call .= "-cport $cport " if ($cport);
	$call .= "-vpass $vpass " if ($vpass);
	$call .= "-toplevel $top_level " if ($top_level);
	#print $call,"\n";
	system($call);      
    }
    $call = "bsub -q normal -K -w 'done($dbname\_variation_job*)' -J waiting_process sleep 1"; #waits until all variation features have finished to continue
    system($call);

    #After post processing, re-calculate map_weight
    $dbVar->do(qq{CREATE TABLE tmp_map_weight
                SELECT variation_id, count(*) as count
                FROM   variation_feature
                $in_string
                GROUP BY variation_id}
	       );
    $dbVar->do(qq{ALTER TABLE tmp_map_weight 
		      ADD INDEX variation_idx(variation_id,count)});

    $dbVar->do(qq{UPDATE variation_feature vf, tmp_map_weight tmw set vf.map_weight = tmw.count
                  WHERE vf.variation_id=tmw.variation_id});
    $dbVar->do(qq{DROP TABLE tmp_map_weight});
}

#will take the number of processes, and divide the total number of entries in the flanking_sequence table by the number of processes
#has to wait until the variation_feature table has been filled
sub parallel_flanking_sequence{
  my $dbVar = shift;  
  my $call;
  my $flanking_status_file = "flanking_status_file_$$\.log";
  #first, create the log file for the variation_feature
  open STATUS, ">$TMP_DIR/$flanking_status_file"
    or throw("Could not open tmp file: $TMP_DIR/$flanking_status_file\n"); 
  close STATUS;

  my $dbname = $vdbname; #get the name of the database to create the file

  #find out the total number of variations to split the into the files
  my $sequences; #total number of flanking sequences in table
  my $sth_variations = $dbVar->prepare(qq{SELECT COUNT(*) from flanking_sequence
 					  });
  $sth_variations->execute();
  ($sequences) = $sth_variations->fetch();
  $sth_variations->finish();

  my $sth = $dbVar->prepare(qq{SELECT fs.variation_id, fs.up_seq, fs.down_seq,
			       vf.seq_region_id, vf.seq_region_start,
			       vf.seq_region_end, vf.seq_region_strand
			       FROM flanking_sequence fs FORCE INDEX (PRIMARY) LEFT JOIN variation_feature vf
			       ON vf.variation_id = fs.variation_id
 			       ORDER BY fs.variation_id
			       $LIMIT},{mysql_use_result => 1});

  $sth->execute();
  my $count = 0; #to know the number of variations
  my $process = 1; #number of file to write the process to
  my $sub_sequences = int($sequences->[0] / $num_processes);
  my $previous_variation_id = 0;
  my $curr_variation_id;
  my $buffer = {}; #hash containing all the files to be parallelized
  #create the files to send to parallelize
  while (my $row = $sth->fetch()){
      $count++; #counting the total number of entries in the table
      $curr_variation_id = $row->[0]; #get the current variation_id
      if ($previous_variation_id == 0){ #initialize the first time
	  $previous_variation_id = $curr_variation_id;
      }
      if ($curr_variation_id ne $previous_variation_id){
	  if((($sub_sequences * $process) <= $count) && ($process < $num_processes)){ #need to write in a new file
	      $process++;
	  }
      }
      my @a = map {defined($_) ? $_ : '\N'} @$row;
      
      &print_buffered($buffer,"$TMP_DIR/$dbname.flanking_sequence_$process\.txt",join("\t",@a) . "\n");
      $previous_variation_id = $curr_variation_id;
      
  }
  $sth->finish();  
  &print_buffered($buffer);
  for (my $i = 1;$i<=$num_processes;$i++){
      $call = "bsub -q normal -m 'bc_hosts' -o $TMP_DIR/output_flanking_$i\_$$.txt /usr/local/ensembl/bin/perl parallel_flanking_sequence.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $flanking_status_file -file $i ";
      $call .= "-cpass $cpass " if ($cpass);
      $call .= "-cport $cport " if ($cport);
      $call .= "-vpass $vpass " if ($vpass);
      system($call);

  }
}

#when the variation_feature table has been filled up, run the variation_group_feature. Not necessary to parallelize as fas as I know....
sub parallel_variation_group_feature{
    my $dbVar = shift;

    my $total_process = 0;
    my $call = "bsub -q normal  -o $TMP_DIR/output_group_feature_$$\.txt /usr/local/ensembl/bin/perl parallel_variation_group_feature.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile $TMP_FILE ";
    $call .= "-cpass $cpass " if ($cpass);
    $call .= "-cport $cport " if ($cport);
    $call .= "-vpass $vpass " if ($vpass);
    $call .= "-limit $limit" if ($limit);
    system($call);
}

#will fill in the transcript variation table. Has to wait until the variation feature table has been filled up. Then, divide the number of entries
# by the number of processes to make the subprocesses
sub parallel_transcript_variation{
    my $dbVar = shift;

    my $total_process = 0;

    my $length_slices = 0; #number of entries in the variation_feature table
    my $call;
    my $transcript_status_file = "transcript_status_file_$$\.log";
    #first, create the log file for the transcript_variation table
    open STATUS, ">$TMP_DIR/$transcript_status_file"
	or throw("Could not open tmp file: $TMP_DIR/$transcript_status_file\n"); 
    close STATUS;
    #then, calculate the rows for each subprocess
    my $sa = $dbCore->get_SliceAdaptor();
    my $inc_non_ref = 1;
    my $slices = $sa->fetch_all('toplevel', undef, $inc_non_ref);
    #order the slices by name
    my @slices_ordered = sort {$a->seq_region_name cmp $b->seq_region_name }  @{$slices};
    # assumes that variation features have already been pushed to toplevel
    foreach my $slice (@slices_ordered) {
      $length_slices += $slice->length;
    }
    
    #I must add up the length of all the slices to find the limit for each 
    my $sub_slice = int($length_slices / $num_processes);
    my $slice_max; #the number of slices in the chunk
    my $slice_min = 0; #first slice in the chunk
    for (my $i = 0; $i < $num_processes ; $i++){
	$length_slices = 0;
	$slice_max = 0;
	foreach my $slice (@slices_ordered){
	    if (($length_slices + $slice->length )<= ($i+1)*$sub_slice){
		$length_slices += $slice->length;
		$slice_max++;	       
	    }
	    else{last;}
	}
	$limit = $slice_min . "," . ($slice_max-$slice_min) if ($i+1 < $num_processes or $num_processes==1);
	$limit = $slice_min . "," . (scalar(@slices_ordered)-$slice_min) 
	  if ($i+1 == $num_processes and $num_processes != 1); #the last slice, get the left slices

	#$call = "bsub -a normal -o $TMP_DIR/output_transcript_$i\_$$.txt  -q normal -m 'bc_hosts' /usr/local/ensembl/bin/perl parallel_transcript_variation.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -limit $limit -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $transcript_status_file ";
	#$call .= "-cpass $cpass " if ($cpass);
	#$call .= "-cport $cport " if ($cport);
	#$call .= "-vpass $vpass " if ($vpass);
	$call = "bsub -a normal -o $TMP_DIR/output_transcript_$i\_$$.txt  -q normal -m 'bc_hosts' /usr/local/ensembl/bin/perl parallel_transcript_variation.pl -species $species -limit $limit -tmpdir $TMP_DIR -tmpfile $TMP_FILE -num_processes $num_processes -status_file $transcript_status_file ";
	$slice_min = $slice_max;
	system($call);
    }
}

##use genotype method to change alleles in genotype table to forward strand, this is for mouse only in order to merge with sanger called SNPs

sub reverse_things {

  my $dbVar = shift;
  #reverse_genotype($dbVar);
  reverse_variation_feature($dbVar);
  #reverse_flanking_sequence($dbVar);
}

sub reverse_genotype{
  my $dbVar = shift;

  my ($population_genotype_id,$variation_id,$allele_id,$allele,$allele_1,$allele_2,$frequency,$sample_id,$allele_string);

  #foreach my $table ("tmp_individual_genotype_single_bp","individual_genotype_multiple_bp","population_genotype","allele") {
    #we only change things that have map_weight=1
    foreach my $table ("test_allele") {
    debug("processling table $table");
    open OUT, ">$TMP_DIR/$table\_out" or die "can't open file $table\_out : $!";
    
    my $sth = $dbVar->prepare(qq{select tg.*,vf.allele_string from $table tg, test_variation_feature vf 
                                 where vf.variation_id=tg.variation_id
                                 #and vf.variation_id in (13,14)
                                 and vf.seq_region_strand = -1 and map_weight=1}, {mysql_use_result=>1} );
    $sth->execute();
    if ($table =~ /pop/i) {
      $sth->bind_columns(\$population_genotype_id,\$variation_id,\$allele_1,\$allele_2,\$frequency,\$sample_id,\$allele_string);
    }
    elsif($table =~ /allele/i) {
      $sth->bind_columns(\$allele_id,\$variation_id,\$allele,\$frequency,\$sample_id,\$allele_string);
    }
    else {
      $sth->bind_columns(\$variation_id,\$allele_1,\$allele_2,\$sample_id,\$allele_string);
    }
    while ($sth->fetch()){
      my ($old_allele,$old_allele_1,$old_allele_2);
      my @alleles = split /\//, $allele_string;
      foreach my $a (@alleles) {
        reverse_comp(\$a);
      }  
      my $new_allele_string = join "/", @alleles;
      
      if ($table =~ /allele/) {
	next if ($allele =~ /\-|\+/);
	next if ($allele =~ /(INDETERMINATE)/i);
        $old_allele = $allele;
        if ($new_allele_string !~ /$allele/ or ($new_allele_string =~ /$allele/ and $allele_string =~ /$allele/)) {
	  if ($allele =~ /\(|\)/g) {
            $allele =~ /\((\S+)\)(\d+)/; print "match is $1 and num is $2\n";
            $allele = $1;
            reverse_comp(\$allele);
            $allele = "($allele)$2";
          }  
          else {
            reverse_comp(\$allele);
          }  
          print OUT "update $table set allele = \"$allele\" where allele_id=$allele_id;\n";
	}
      }
      else {
	next if ($allele_1 =~ /\-|\+/ and $allele_2 =~ /\-|\+/);
	next if ($allele_1 =~/\+/ or $allele_2 =~/\+/);
	next if ($allele_1 =~/homozygous|indeterminate|SEQUENCE/ig or $allele_2 =~/homozygous|indeterminate|SEQUENCE/ig);
        $old_allele_1 = $allele_1;
	$old_allele_2 = $allele_2;
	if (($new_allele_string !~ /$allele_1/ and $new_allele_string !~ /$allele_2/) or ($new_allele_string =~ /$allele_1|$allele_2/ and $allele_string =~ /$allele_1|$allele_2/)) {
          if ($allele_1 =~ /\(|\)/g) {
            $allele_1 =~ /\((\S+)\)(\d+)/; 
            $allele_1 = $1;
            reverse_comp(\$allele_1);
            $allele_1 = "($allele_1)$2";
          }
          else {
            reverse_comp(\$allele_1);
          }
          if ($allele_2 =~ /\(|\)/g) {
            $allele_2 =~ /\((\S+)\)(\d+)/; 
            $allele_2 = $1;
            reverse_comp(\$allele_2);   
            $allele_2 = "($allele_2)$2";      
          }                                               
         else {
           reverse_comp(\$allele_1);   
         }
	 print OUT "update $table set allele_1 = \"$allele_1\", allele_2 = \"$allele_2\" where variation_id=$variation_id and sample_id=$sample_id and allele_1 = \"$old_allele_1\" and allele_2 = \"$old_allele_2\";\n";
        }
      }   
    }
    close OUT;
    #system("mysql -uensadmin -pensembl -h $vhost $vdbname <$TMP_DIR/$table\_out");
  }
}

sub reverse_variation_feature{
  my $dbVar = shift;

  debug("processing reverse variation_feature strand");

  my ($variation_feature_id,$allele_string);

  my $sth = $dbVar->prepare(qq{select variation_feature_id,allele_string from test_variation_feature where seq_region_strand=-1 and map_weight=1});
  $sth->execute();
  $sth->bind_columns(\$variation_feature_id,\$allele_string);

  while ($sth->fetch()){
    my @alleles = split /\//,$allele_string;
    foreach my $allele (@alleles) {
      if ($allele =~ /\((\S+)\)(\d+)/) {
        $allele = $1;                    
        reverse_comp(\$allele);                    
        $allele = "($allele)$2";
      }
      else {
        reverse_comp(\$allele);
      }  
    }
    my $new_allele_string = join "/",@alleles;
    $dbVar->do(qq{update variation_feature set allele_string = "$new_allele_string", seq_region_strand =1 where variation_feature_id=$variation_feature_id});
  }
}

sub reverse_flanking_sequence {

  my $dbVar = shift;

  debug("Processing flanking_sequence reverse strand");

  my $sth=$dbVar->prepare(qq{select * from flanking_sequence where seq_region_strand=-1});
  $sth->execute();
  my ($variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand);

  $sth->bind_columns(\$variation_id,\$up_seq,\$down_seq,\$up_seq_start,\$up_seq_end,\$down_seq_start,\$down_seq_end,\$seq_region_id,\$seq_region_strand);
  while($sth->fetch()) {
    #print "$variation_id,$up_seq,$down_seq,$up_seq_start,$up_seq_end,$down_seq_start,$down_seq_end,$seq_region_id,$seq_region_strand\n";
    if ($up_seq and $down_seq) {
      ($up_seq, $down_seq) = ($down_seq, $up_seq);
      reverse_comp(\$up_seq);
      reverse_comp(\$down_seq);
      $up_seq_start = $up_seq_end = $down_seq_start = $down_seq_end = '\N';
    }
    elsif (! $up_seq and ! $down_seq) {
      my $tmp_seq_start = $up_seq_start;
      my $tmp_seq_end = $up_seq_end;
      ($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
      ($down_seq_start, $down_seq_end) = ($tmp_seq_start, $tmp_seq_end);
      $up_seq = $down_seq = '\N';
    }
    elsif ($up_seq and ! $down_seq) {
      $down_seq = $up_seq;
      reverse_comp(\$down_seq);
      $up_seq = '\N';
      ($up_seq_start, $up_seq_end) = ($down_seq_start, $down_seq_end);
      $down_seq_start = '\N';
      $down_seq_end = '\N';
    }
    elsif (! $up_seq and $down_seq) {
      $up_seq = $down_seq;
      reverse_comp(\$up_seq);
      $down_seq = '\N';
      ($down_seq_start, $down_seq_end) = ($up_seq_start, $up_seq_end);
      $up_seq_start = '\N';
      $up_seq_end = '\N';
    }
    $seq_region_strand = 1 if ($seq_region_strand == -1);
    #my $sth1 = $dbVar->prepare(qq{update flanking_sequence set up_seq = ?, down_seq = ?, up_seq_region_start = $up_seq_start, up_seq_region_end = $up_seq_end, down_seq_region_start=$down_seq_start, down_seq_region_end=$down_seq_end, seq_region_strand=1  where variation_id = $variation_id and seq_region_id = $seq_region_id and seq_region_strand = -1});
    #print "up_seq is $up_seq and down_seq is $down_seq\n";
    #$sth1->bind_param(1,$up_seq,SQL_VARCHAR) if ($up_seq ne '\N');
    #$sth1->bind_param(1,'NULL',SQL_VARCHAR) if ($up_seq eq '\N');
    #$sth1->bind_param(2,$down_seq,SQL_VARCHAR) if ($down_seq ne '\N');
    #$sth1->bind_param(2,'NULL',SQL_VARCHAR) if ($down_seq eq '\N');
    #$sth1->execute();
    $dbVar->do(qq{update flanking_sequence set up_seq = "$up_seq", down_seq = "$down_seq", up_seq_region_start = $up_seq_start, up_seq_region_end = $up_seq_end, down_seq_region_start=$down_seq_start, down_seq_region_end=$down_seq_end, seq_region_strand=1  where variation_id = $variation_id and seq_region_id = $seq_region_id and seq_region_strand = -1});
  }
  $dbVar->do(qq{update flanking_sequence set up_seq = null where up_seq = 'N'});
  $dbVar->do(qq{update flanking_sequence set down_seq = null where down_seq = 'N'});
}

#will have to wait until the variation_feature has finished. Then, select all the genotype, and split the data into files (1 per population)
sub parallel_ld_populations {
    my $dbVar = shift;    
    my $call;

    my %seq_region; #hash containing the mapping between seq_region_id->name region
    my %alleles_variation = (); #will contain a record of the alleles in the variation. A will be the major, and a the minor. When more than 2 alleles
    my %genotype_information; #will contain all the genotype information to write in the file, if necessary
    #, the genotypes for that variation will be discarded
    my %regions; #will contain all the regions in the population and the number of genotypes in each one
    my $previous_variation_id = ''; #to know if it is a new variation and we can get the new alleles
    my $buffer = {}; #will contain a buffer where will be written all the LD information
    my %populations;
    my %genotypes_file; #foreach file, will contain the number of genotypes, so we can split it later

    #going to get the population_id for the HapMap and PerlEgene populations and a hash with the individuals that shouldn't
    #be present in the LD calculation
    my $pop_id;
    my $population_name;
    #get all populations to be tagged (HapMap and PerlEgen)
    my $sth = $dbVar->prepare(qq{SELECT s.sample_id, s.name
				     FROM population p, sample s
				     WHERE (s.name like 'PERLEGEN:AFD%'
				     OR s.name like 'CSHL-HAPMAP%')
				     AND s.sample_id = p.sample_id
				 });
    

    $sth->execute();
    $sth->bind_columns(\$pop_id,\$population_name);
    #get all the children that we do not want in the genotypes
    my $siblings = {}; # hash {$individual_id} ,where the individual is sibling of another one
    my @pops;
    while($sth->fetch){
	if($population_name =~ /CEU|YRI/){
	    &get_siblings($dbVar,$pop_id,$siblings);
	}
	push @pops, $pop_id;
    }
 
    my $in_str = " IN (" . join(',', @pops). ")";

    #necessary the order to know when we change variation. Not get genotypes with a NULL variation or map_weight > 1

    $sth = $dbVar->prepare
	(qq{
	    SELECT  STRAIGHT_JOIN ig.variation_id, vf.variation_feature_id, vf.seq_region_id, vf.seq_region_start, 
                          ig.sample_id, ig.allele_1, ig.allele_2, vf.seq_region_end, ip.population_sample_id
		    FROM  variation_feature vf FORCE INDEX(pos_idx), individual_genotype_single_bp ig, individual_population ip
		   WHERE  ig.variation_id = vf.variation_id

		    AND   ig.allele_2 IS NOT NULL
		    AND   vf.map_weight = 1
		    AND   ip.individual_sample_id = ig.sample_id
		    AND   ip.population_sample_id $in_str
		    ORDER BY  vf.seq_region_id,vf.seq_region_start}, {mysql_use_result => 1} );


    print "Time starting to dump data from database: ",scalar(localtime(time)),"\n";
    $sth->execute();

    my $dbname = $vdbname; #get the name of the database to create the file
    my ($variation_id, $variation_feature_id, $seq_region_id, $seq_region_start, 
	$individual_id, $allele_1,$allele_2,$seq_region_end,$population_id);


    $sth->bind_columns(\$variation_id, \$variation_feature_id, \$seq_region_id, \$seq_region_start, 
		       \$individual_id, \$allele_1,\$allele_2,\$seq_region_end,\$population_id); 
    while ($sth->fetch()){
	if ($previous_variation_id eq ''){
	    $previous_variation_id = $variation_id;
	}
	#only print genotypes without parents genotyped
	if (!exists $siblings->{$population_id . '-' . $individual_id}){ #necessary to use the population_id
	    #if it is a new variation, write to the file (if necessary) and empty the hash
	    if ($previous_variation_id ne $variation_id){
		foreach my $population (keys %alleles_variation){
		    #if the variation has 2 alleles, print all the genotypes to the file
		    if (keys %{$alleles_variation{$population}} == 2){		
			&convert_genotype($alleles_variation{$population},$genotype_information{$population});
			foreach my $individual_id (keys %{$genotype_information{$population}}){
			    &print_individual_file($buffer,$population, 
						   $previous_variation_id, $individual_id,
						   \%genotype_information,$dbname,\%genotypes_file,\%regions);
			}
		    }
		}
		$previous_variation_id = $variation_id;
		%alleles_variation = (); #new variation, flush the hash
		%genotype_information = (); #new variation, flush the hash
	    }
	    #we store the genotype information for the variation
	    if ($allele_1 ne 'N' and $allele_2 ne 'N'){
	      $genotype_information{$population_id}{$individual_id}{variation_feature_id} = $variation_feature_id;
	      $genotype_information{$population_id}{$individual_id}{seq_region_start} = $seq_region_start;
	      $genotype_information{$population_id}{$individual_id}{allele_1} = $allele_1;
	      $genotype_information{$population_id}{$individual_id}{allele_2} = $allele_2;
	      $genotype_information{$population_id}{$individual_id}{seq_region_end} = $seq_region_end;
	      $genotype_information{$population_id}{$individual_id}{seq_region_id} = $seq_region_id;
	      
	      #and the alleles
	      $alleles_variation{$population_id}{$allele_1}++;
	      $alleles_variation{$population_id}{$allele_2}++;
	    
	      $populations{$population_id}++;
	    }
	}
    }
    $sth->finish();
    #we have to print the last variation
    foreach my $population (keys %alleles_variation){
	#if the variation has 2 alleles, print all the genotypes to the file
	if (keys %{$alleles_variation{$population}} == 2){		
	    &convert_genotype($alleles_variation{$population},$genotype_information{$population});
	    foreach my $individual_id (keys %{$genotype_information{$population}}){
		&print_individual_file($buffer,$population, 
				       $previous_variation_id, $individual_id,
				       \%genotype_information,$dbname,\%genotypes_file,\%regions);
	    }
	}
    }

    print_buffered( $buffer );
    print "Time starting to submit jobs to queues: ",scalar(localtime(time)),"\n";
    #let's run a job array
    foreach my $file (keys %genotypes_file){
#	if ($genotypes_file{$file} > MAX_GENOTYPES()){
#	    &split_file($file,\%regions,$genotypes_file{$file},$dbname);
#	    $call = "bsub -q normal -J '$dbname.pairwise_ld[1-".REGIONS()."]' -m 'bc_hosts' -o $TMP_DIR/output_pairwise_ld.txt /usr/local/ensembl/bin/perl calc_genotypes.pl $file\_chunk.\\\$LSB_JOBINDEX $file\_chunk_out.\\\$LSB_JOBINDEX"; #create a job array
#	    $call = "bsub -q normal -J '$dbname.pairwise_ld[1-".REGIONS()."]' -m 'bc_hosts' ./ld_wrapper.sh $file\_chunk.\\\$LSB_JOBINDEX $file\_chunk_out.\\\$LSB_JOBINDEX"; #create a job array
#	}
#	else{
#	    $call = "bsub -q normal -J '$dbname.pairwise_ld' -m 'bc_hosts' -o $TMP_DIR/output_ld_populations\_$$.txt /usr/local/ensembl/bin/perl calc_genotypes.pl $file $file\_out ";	    
	    $call = "bsub -q normal -J '$dbname.pairwise_ld' -m 'bc_hosts' ./ld_wrapper.sh $file $file\_out";	    
#	}
	system($call);
    }
    $call = "bsub -q normal -w 'done($dbname.pairwise_ld)' -m 'ecs4_hosts' -o $TMP_DIR/output_ld_populations_import.txt /usr/local/ensembl/bin/perl parallel_ld_populations.pl -chost $chost -cuser $cuser -cdbname $cdbname -vhost $vhost -vuser $vuser -vport $vport -vdbname $vdbname -tmpdir $TMP_DIR -tmpfile $TMP_FILE ";
    $call .= "-cpass $cpass " if ($cpass);
    $call .= "-cport $cport " if ($cport);
    $call .= "-vpass $vpass " if ($vpass);
#    system($call); #send the last job, that will wait for the rest to finish and upload everything to the table 

}


#prints to the buffer the different individuals
sub print_individual_file{
    my $buffer = shift;
    my $population = shift;
    my $previous_variation_id = shift;
    my $individual_id = shift;
    my $genotype_information = shift;
    my $dbname = shift;
    my $genotypes_file = shift;
    my $regions = shift;


    $regions->{"$TMP_DIR/$dbname.pairwise_ld_$population\.txt"}->{$genotype_information->{$population}->{$individual_id}->{seq_region_id}}++; #add one more genotype in the region for the population

    print_buffered( $buffer, "$TMP_DIR/$dbname.pairwise_ld_$population\.txt", 
		    join("\t",$previous_variation_id,$genotype_information->{$population}->{$individual_id}->{seq_region_start},
			 $individual_id,
			 $genotype_information->{$population}->{$individual_id}->{genotype},
			 $genotype_information->{$population}->{$individual_id}->{variation_feature_id},
			 $genotype_information->{$population}->{$individual_id}->{seq_region_id},
			 $genotype_information->{$population}->{$individual_id}->{seq_region_end},$population)."\n" );
    $genotypes_file->{"$TMP_DIR/$dbname.pairwise_ld_$population\.txt"}++;

}

sub print_buffered {
    my $buffer = shift;
    my $filename = shift;
    my $text = shift;

    local *FH;

    if( ! $filename ) {
	# flush the buffer
	foreach my $file (keys %{$buffer}){
	    open( FH, ">>$file" ) or die "Could not print to file $! \n";
	    print FH $buffer->{ $file };
	    close FH;
	}

    } else {
	$buffer->{ $filename } .= $text;
	if( length( $buffer->{ $filename } ) > 10_000 ) {
	    open( FH, ">>$filename" ) or die;
	    print FH $buffer->{ $filename };
	    close FH;
	    $buffer->{ $filename } = '';
	}
    }
}


#
# Converts the genotype into the required format for the calculation of the pairwise_ld value: AA, Aa or aa
# From the Allele table, will select the alleles and compare to the alleles in the genotype
#
sub convert_genotype{
    my $alleles_variation = shift; #reference to the hash containing the alleles for the variation present in the genotypes
    my $genotype_information = shift; #reference to a hash containing the values to be written to the file
    my @alleles_ordered; #the array will contain the alleles ordered by apparitions in the genotypes (only 2 values possible)
    
    @alleles_ordered = sort({$alleles_variation->{$b} <=> $alleles_variation->{$a}} keys %{$alleles_variation});

    #let's convert the allele_1 allele_2 to a genotype in the AA, Aa or aa format, where A corresponds to the major allele and a to the minor
    foreach my $individual_id (keys %{$genotype_information}){
	#if both alleles are different, this is the Aa genotype
	if ($genotype_information->{$individual_id}{allele_1} ne $genotype_information->{$individual_id}{allele_2}){
	    $genotype_information->{$individual_id}{genotype} = 'Aa';
	}
	#when they are the same, must find out which is the major
	else{	    
	    if ($alleles_ordered[0] eq $genotype_information->{$individual_id}{allele_1}){
		#it is the major allele
		$genotype_information->{$individual_id}{genotype} = 'AA';
	    }
	    else{
		$genotype_information->{$individual_id}{genotype} = 'aa';
	    }
	    
	}
    }
}

#given a file with more than MAX_GENOTYPES, splits it in REGIONS different files
sub split_file{
    my $file = shift;
    my $regions = shift; #hash with all the regions in the population and the number of genotypes in each region
    my $genotypes = shift;
    my $dbname = shift;
    my @regions_ordered = sort {$a<=>$b} keys %{$regions->{$file}};
    my $sub_genotypes = int($genotypes/REGIONS()); #find minimum number of genotypes in each file
    #create the groups of regions for each file: the array will contain the number of genotypes (lines) that the group must contain
    my @groups;
    my $lines = 0;
    my $index = 1; #position in the group. The first position will contain n lines, the second $i*$n,...
    foreach my $region (@regions_ordered){
	$lines += $regions->{$file}->{$region};
	if ($lines > $sub_genotypes * $index){
	    push @groups,$lines;	    
	    $index++;
	}
    }
    open INFILE, "< $file" or die "Could not open file $file: $!";
    my $chunk = 1;
    until(eof INFILE) {
	open OUTFILE, "> $file\_chunk.$chunk"
	    or die "Could not open file $TMP_DIR/$dbname.pairwise_ld_chunk_$$\_$chunk\.txt: $!";	
	while(<INFILE>) {
	    print OUTFILE;
	    if ($chunk != REGIONS()){
		last unless $. % $groups[$chunk-1]; #last chunk will contain all remaining lines in the file
	    }
	}
	++$chunk;
	close OUTFILE;
    }
    close INFILE;
    unlink $file; #remove the original file after splitting it
}

#for a given population, gets all individuals that are children (have father or mother)
sub get_siblings{
    my $dbVariation = shift;
    my $population_id = shift;
    my $siblings = shift;

    my $sth_individual = $dbVariation->prepare(qq{SELECT i.sample_id
							     FROM individual i, individual_population ip
							     WHERE ip.individual_sample_id = i.sample_id
							     AND ip.population_sample_id = ? 
							     AND i.father_individual_sample_id IS NOT NULL
							     AND i.mother_individual_sample_id IS NOT NULL
							 });
    my ($individual_id);
    $sth_individual->execute($population_id);
    $sth_individual->bind_columns(\$individual_id);
    while ($sth_individual->fetch){
	$siblings->{$population_id.'-'.$individual_id}++; #necessary to have in the key the population, since some individuals are shared between
	                                                   #populations
    }
    return $siblings;
}


#method to crete a tmp table to store failed variations
sub create_failed_variation_table{
    my $dbVar = shift;

    $dbVar->do(qq{CREATE TABLE IF NOT EXISTS failed_variation(
			  variation_id int(10) unsigned not null,
			  failed_description_id int(10) unsigned not null,

			   PRIMARY KEY(variation_id))
		  }
	       );
}


sub usage {
  my $msg = shift;

  print STDERR <<EOF;

usage: perl parallel_post_process.pl <options>

options:
    -chost <hostname>    hostname of core Ensembl MySQL database (default = ecs2)
    -cuser <user>        username of core Ensembl MySQL database (default = ensro)
    -cpass <pass>        password of core Ensembl MySQL database
    -cport <port>        TCP port of core Ensembl MySQL database (default = 3364)
    -cdbname <dbname>    dbname of core Ensembl MySQL database
    -vhost <hostname>    hostname of variation MySQL database to write to
    -vuser <user>        username of variation MySQL database to write to (default = ensadmin)
    -vpass <pass>        password of variation MySQL database to write to
    -vport <port>        TCP port of variation MySQL database to write to (default = 3306)
    -vdbname <dbname>    dbname of variation MySQL database to write to
    -limit <number>      limit the number of rows for testing
    -tmpdir <dir>        temp directory to use (with lots of space!)
    -tmpfile <filename>   name of temp file to use
    -num_processes <number> number of processes that are running (default = disabled)
    -variation_feature  fill in the Variation_feature table (default = disabled)
    -flanking_sequence  fill in the flanking sequence tables (default = disabled)
    -variation_group_feature fill in the Variation_group_feature table (default = disabled)
    -transcript_variation  fill in the Transcript_variation table (default = disabled)
    -ld_populations  fill in the Pairwise_ld table (default = disabled)
EOF

  die("\n$msg\n\n");
}
