#!/usr/bin/env perl
$|=1;
use strict;

use lib "/bigdata0/famgen/FAVA";
use FAVA;

my $favaVersion = "/bigdata0/famgen/FAVA/latest";
my $frz = 'hg19';
my $fava = new FAVA $frz, $favaVersion;

my $tabix = `which tabix`;
chomp($tabix);
die "Couldn't find tabix\n" unless ($tabix);

my $bgzip = `which bgzip`;
chomp($bgzip);
die "Couldn't find bgzip\n" unless ($bgzip);

my $outdir = "fulloutdir1";
mkdir $outdir, 0755;
die "Couldn't create outdir: $!\n" unless -e $outdir;

#grabs list of sample names
my $wanted_samples = grab_list($ARGV[0]); 

my @CHROM = split(/,/, $ARGV[1]);
 
#grabs only the sample ids that have a genome.rcf.gz file
my $fgs = "/bigdata0/users/dmauldin/CGIRCF/";
my($sp, $ind_dir) = identifyGenomes_list($wanted_samples, $fgs);

my @ind = sort keys %{$sp}; #individual ids 

my $count = scalar(@ind);

foreach my $chr (@CHROM){
    my %todo;
    my $chrom = "chr".$chr;
    my $start = 1;
    my $end   = 1e9;
    my $done;
    foreach my $ind (@ind){
	my $file = "$ind_dir->{$ind}$ind.rcf.gz";
	check_file_and_index($file);
	open TBX,"$tabix $ind_dir->{$ind}$ind.rcf.gz $chrom:$start-$end |  awk '{ if (\$6!~/\\./) print}' |" || die "Could not open $ind_dir->{$ind}$ind.rcf.gz\n";
	my $dind = $ind;
	while(my $line = <TBX>){
	    chomp ($line);
	    my (undef, $pos, undef, undef, $ref, $var1) = split(/\t/, $line);
	    my @alleles = split(/,/, $var1);
	    foreach my $allele (@alleles){
		if ($allele ne $ref){
		        #print join("\t",$pos, $allele, $dind)."\n";
		    $todo{$pos}{$allele}{$dind}++;
		}
	    }
	}
	close TBX;
	$done++;
	print "." unless $done % 10;    
    }
    print "\n";
    
    my @snvs = sort {$a<=>$b} keys %todo;
    my $snvs = scalar @snvs;
    print "$chrom\t$snvs <=== NUMBER OF SNVS IN THIS CHROMOSOME *******************\n";
    my $done;
    open BIGOUT, ">$outdir/raw.$chrom";
    foreach my $pos (@snvs) {
	foreach my $var (keys %{$todo{$pos}}) {
	    testEffect($todo{$pos}{$var}, $chrom, $pos, $var);
	}
	$done++;
	print "." unless $done % 1000;
	unless ($done % 100000) {
	    print " ", sprintf("%.1f", 100*$done/$snvs), "%\n";
	}
    }
    print "done\n";
    system ("echo \"Finished running $ARGV[1]\" | mail -s \"$ARGV[1] is completed\" crystal.humphries\@systemsbiology.org");

    close BIGOUT;
    `nice -5 sed -i 's/vcfBeta-//g' $outdir/raw.$chrom`;
    `nice -5 $bgzip -f $outdir/raw.$chrom`;
    `nice -5 $tabix $outdir/raw.$chrom.gz -s 1 -b 2 -e 2 -f`;
}

sub check_file_and_index{
    my $file     = shift;
    my @location = split(/\//, $file);
    my $rcf      = pop (@location);
    my $tbi      = $rcf.".tbi";
    my $flag     = 0;
    
    #verify existence of RCF
    $flag       = "true" if (-e $file); 
    my $idx_loc  = join("/", @location, $tbi);
 
    #create index file
    if ( (!-e $idx_loc) or (-z $idx_loc) ){
	my $wd       = `pwd`;
	my $dir      = join("/", @location);
	chdir $dir; 
	
	if (0 == system("tabix -p vcf $rcf")) {
	    print "tabix created index\n";
	} 
	else {
	    print "Something went wrong with calling tabix\n";
	}
	chdir $wd;
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
	$sample = "vcfBeta-".$sample; 
	my $new_dir = $fgs;
        my $file_loc = $new_dir.$sample.".rcf.gz";
	if (-e $file_loc){
	    $sp{$sample}=+1;
	    $dir{$sample}=$new_dir;
	}
    }
    return (\%sp, \%dir);
}

sub testEffect {
    my($ind, $chrom, $pos, $var) = @_;
    my @ind = sort {$a<=>$b} keys %{$ind};
    my $inds = scalar @ind;
    
    my %res = $fava->testSite($chrom, $pos-1, $var);

    
    if ( my @siteEffects = sort keys %{$res{'site'}}){
	print BIGOUT join("\t", $chrom, $pos, $var, $inds, join(",", @ind));
	foreach my $eff (@siteEffects) {
	    print BIGOUT "\t".join("\t", $eff, $res{'site'}{$eff}{$var}); #,$CADD_raw, $PHRED);
	}
    }else{
	print BIGOUT join("\t", $chrom, $pos, $var, $inds, join(",", @ind), "no_FAVA_impact");
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
