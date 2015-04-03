#!/usr/bin/env perl
$|=1;
use strict;
use Storable;
use lib "/proj/famgen/resources/FAVA";
use FAVA;
use POSIX qw(strftime);

my $favaVersion = "/proj/famgen/resources/FAVA/latest";
my $frz = 'hg19';
my $fava = new FAVA $frz, $favaVersion;

my $tabix = `which tabix`;
chomp($tabix);
die "Couldn't find tabix\n" unless $tabix;

my $bgzip = `which bgzip`;
chomp($bgzip);
die "Could not find bgzip\n" unless ($bgzip);

my $outdir = "fulloutdir_CGI";
mkdir $outdir, 0755;
die "Couldn't create outdir: $!\n" unless -e $outdir;

my %todo;
my %chr_lengths = (1 => 2.5e8, 2 => 2.4e8, 3 => 2.0e8, 4 =>1.9e8, 5 => 1.8e8, 6 => 1.7e8, 7 => 1.6e8, 8 => 1.5e8, 9 => 1.4e8, 10 =>1.4e8,
                         11 => 1.4e8, 12=> 1.3e8, 13 => 1.2e8, 14 => 1.1e8, 15 => 1.0e8, 16 => 9.1e7, 17 => 8.2e7, 18 => 7.9e7, 19 => 6.0e7,
                   20 => 6.4e7, 21 => 4.9e7, 22 => 5.2e7, 'X' => 1.6e8, 'Y' =>6.0e7, 'M' => 1.0e6);

#grabs list of sample names
my $wanted_samples = grab_list($ARGV[0]); 

my @CHROM = split(/,/, $ARGV[1]);
 
#grabs only the sample ids that have a genome.rcf.gz file
my $fgs = "/40TB_3/gestalt/private/IndividualGenomeAssembly/";
my($sp, $ind_dir) = identifyGenomes_list($wanted_samples, $fgs);

my @ind = sort keys %{$sp}; #individual ids 
my $count = scalar(@ind);
my $Logfile = create_log();

open my $LOG, '>', $Logfile || die ("Could not open $Logfile\n");


print $LOG "
##################################################################################
Number of Individual Files: $count
##################################################################################\n";

print $LOG $0." ".join(" ", @ARGV)."\n"; 

   

foreach my $chr (@CHROM){
    my $chrom = "chr".$chr;
    #break the chromosome up into                                                                                                                                                                                                       

    my $start = 0;
    my $end   = 1e9;
    my @interval = ($end);

    if (($chr<= 20) and ($chr!~m/X|Y/i)){
        undef (@interval);
        $end   = $chr_lengths{$chr};
        my $interval = $end/20;
        for my $int (1..20){
            push (@interval, int($int*$interval));
        }
    }

    print $LOG my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "\nCurrent Chr: $chr\n";
    print $LOG "Number of intervals:".scalar(@interval)."\n";

    my @interval_files;
    my $n=0;

    foreach my $int (@interval){
        $n++;
        my $new_end = $int;
        my $seen_samples = run_main($chr, $start, $new_end);
        my $interval_file_name = join(".","simple.chr$chr", $n);
	print $LOG "Currently running chr: $start - $end\n";
        print_interval($chrom, $interval_file_name, $seen_samples);
        $start = $new_end + 1;
	my $dir_file = join("/", $outdir,$interval_file_name);
        push (@interval_files, $dir_file); 
    }
    my $final_file = $outdir."/simple.chr$chr";
    print $LOG "cat @interval_files > $final_file\n";
    system("cat @interval_files > $final_file\n");
    system("rm @interval_files");
    system("nice \-5 $bgzip $final_file \-f");
    system ("echo \"Finished running $0 on chr=> $chr\" | mail -s \"$chr is completed\" crystal.humphries\@systemsbiology.org");
}


sub run_main {
    my ($chr, $start, $end) = @_;
    my %todo;
    my $chrom = "chr".$chr;
    my $done;
    my $n = 0;
    print $LOG my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "\nCurrent Chr: $chr\n";
  LOOP:  foreach my $ind (@ind){
      my $file = "$ind_dir->{$ind}Genotype.rcf.gz";
      next if (check_file_and_index_and_version($file) ne "TRUE");
      open TBX,"$tabix $ind_dir->{$ind}Genotype.rcf.gz $chrom:$start-$end |  awk '{ if (\$6!~/\\./) print}' |" || die "Could not open $ind_dir->{$ind}Genotype.rcf.gz\n";
      my $dind = $ind;
      $n++;
      while(my $line = <TBX>){
	  chomp ($line);
	  my (undef, $pos, undef, undef, $ref, $var1) = split(/\t/, $line);
	  next if (length($var1)>1);
	  my @alleles = split(/,/, $var1);
	  foreach my $allele (@alleles){
	      if (($allele ne $ref) and ($allele!~m/\./) and ($allele)){
		  $todo{$pos}{$allele}{$dind}++;
	      }
	  }
      }
      close TBX;
      $done++;
      print "." unless $done % 10;    
  }
    print $LOG "number of opened files: $n\n";
    my $snvs = keys %todo;
    #print $LOG "$chrom\t$snvs <=== NUMBER OF SNVS IN THIS CHROMOSOME *******************\n";
    print $LOG "$chrom\t$snvs <=== Number of SNVS in the following porition of the chromosome: $start\t$end\n";
    return (\%todo);
}

sub print_interval{
    my ($chrom, $output_file, $hash) = @_; 
    my $done;
    my $complete_file = join("/", $outdir, $output_file);
    open BIGOUT, '>', $complete_file;
    my @snvs = sort {$a<=>$b} keys %{$hash};
    my $snvs = scalar(@snvs);
      foreach my $pos (@snvs) {
	foreach my $var (keys %{$hash->{$pos}}) {

	     my @ind = sort {$a<=>$b} keys %{$hash->{$pos}{$var}};
	     my $inds = scalar @ind;
	     print BIGOUT join("\t", $chrom, $pos, $var, $inds, join(",", @ind))."\n";
	     
	     $done++;
	     print "." unless $done % 1000;
	     unless ($done % 100000) {
		 print " ", sprintf("%.1f", 100*$done/$snvs), "%\n";
	     }
	}
    }
    print "done\n";
    close BIGOUT;
}

sub create_log{
    use POSIX qw(strftime);

    my $now_string = strftime "%e%b%Y", localtime;
    my $current = $0;
    $current=~s/\.pl//;
    $now_string=~s/\s+//g;
    my $log_F_b = my $temp = join('.', $now_string,$current,"OUT");
    my $log_F = join("/", "script_logs",$log_F_b);
    my $n = 0;
    until (!-e $temp){
	$n++;
	$temp = join(".",$log_F, $n); 
    }
    return ($temp);
}
    

sub check_file_and_index_and_version{
    my $file     = shift;
    my @location = split(/\//, $file);
    my $rcf      = pop (@location);
    my $tbi      = $rcf.".tbi";
    my $flag     = 0;
    my $index    = join("/", @location, $tbi);
    
    #verify existence of RCF
    if ((-e $file) and (-e $index) and (-r $file) and (-r $index))
    {
	#verify that the CGAPipeline used was version 2.* and higher
	my $line =  `zcat $file | head \-20 | grep 'CGAPipeline_'`;
	$line=~m/_(\d{1})/;
        $flag = "TRUE" if ($1>=2);
    }
    return($flag);
}

sub grab_list{
    # grabs file containing list of sample_ids and returns list as referenced array
    my $file = shift;
    open (my $list, $file) || die "Could not open $file\n";
    my @samples;

    while(my $line = <$list>){
	chomp ($line);
	push (@samples, $line);
    }
    return (\@samples);
}

sub identifyGenomes_list {
    my $list = shift;
    my($fgs) = @_;
    my(%sp, %dir);

    #grabs samples from directory, looks further into samples from list, and 
    #returns list of samples that have a genome.vcf.gz file
    foreach my $sample (@{$list}){
	chomp ($sample);
	my $new_dir = $fgs.$sample.'/analyses/';
        my $file_loc = $new_dir."Genotype.rcf.gz";
	if (-e $file_loc){
	    $sp{$sample}=+1;
	    $dir{$sample}=$new_dir;
	}
    }
    return (\%sp, \%dir);
}

sub testEffect {
	my($ind, $chrom, $pos, $var, $quals) = @_;
	my @ind = sort {$a<=>$b} keys %{$ind};
	my $inds = scalar @ind;
	
	my %res = $fava->testSite($chrom, $pos-1, $var);
	
	if ( my @siteEffects = sort keys %{$res{'site'}}){
	    print BIGOUT join("\t", $chrom, $pos, $var, $inds, join(",", @ind), $quals);
	    foreach my $eff (@siteEffects) {
		print BIGOUT "\t".join("\t", $eff, $res{'site'}{$eff}{$var});
	    }
	}else{
	    print BIGOUT join("\t", $chrom, $pos, $var, $inds, join(",", @ind), $quals, "no_FAVA_impact");
	}   
	print BIGOUT "\n";
}

sub get_CAAD{
    my ($chrom, $pos, $var) = @_;
    my $CADD_genome = '/proj/famgen/resources/CADD/whole_genome_SNVs.tsv.gz';
    $chrom=~s/chr//;
    open TBX, "$tabix $CADD_genome $chrom:$pos-$pos |";
    while(my $line = <TBX>){
	chomp ($line);
	my (undef, undef, $ref, $alt, $CADD, $Phred)=split(/\t/,$line);
	if ($var eq $alt){
	    return ($CADD, $Phred);
	}
    }
}

my %chr_lengths = (1 => 2.5e8, 2 => 2.4e8, 3 => 2.0e8, 4 =>1.9e8, 5 => 1.8e8, 6 => 1.7e8, 7 => 1.6e8, 8 => 1.5e8, 9 => 1.4e8, 10 =>1.4e8,
		      11 => 1.4e8, 12=> 1.3e8, 13 => 1.2e8, 14 => 1.1e8, 15 => 1.0e8, 16 => 9.1e7, 17 => 8.2e7, 18 => 7.9e7, 19 => 6.0e7, 
		   20 => 6.4e7, 21 => 4.9e7, 22 => 5.2e7, 'X' => 1.6e8, 'Y' =>6.0e7, 'M' => 1.0e6);
########
sub fulldirlist {
	my($dir) = @_;
	opendir (DIR, $dir);
	my @files = grep /^[^.]/, readdir DIR;
	closedir DIR;
	return @files;
}

sub slicedirlist {
	my($dir, $pat) = @_;
	opendir (DIR, $dir);
	my @files = grep /$pat/, readdir DIR;
	closedir DIR;
	return @files;
}

