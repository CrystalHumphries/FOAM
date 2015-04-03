#!/usr/bin/env perl
$|=1;
use strict;

#grabs list of sample names
my $wanted_samples = grab_list($ARGV[0]); 

my @CHROM = split(/,/, $ARGV[1]);
 
#grabs only the sample ids that have a genome.rcf.gz file
my $fgs = "/40TB_3/gestalt/private/IndividualGenomeAssembly/";
my($sp, $ind_dir) = identifyGenomes_list($wanted_samples, $fgs);

my @ind = sort keys %{$sp}; #individual ids 

my $count = scalar(@ind);
print "Number of Individual Files: $count\n";


foreach my $chr (@CHROM){
    my $chrom = "chr".$chr;
    my $n = 0;
    foreach my $ind (@ind){
	my $file  = "$ind_dir->{$ind}Genotype.rcf.gz";
	my $e_gen = my $r_gen = my $e_index = my $r_index = "noAccess";
	$e_gen = "FileExists"   if (-e $file);
	$r_gen =  "CanRead"      if (-r $file);
#	print join("\t", $file, $e_gen, $r_gen)."\n";
	my $index = check_index($file);
	$e_index = "FileExists" if (-e $index);
	$r_index = "CanRead"    if (-r $index);
#	my $V    = index_version($index);
	print join("\t", $index, $e_index, $r_index)."\n";
    }
}

sub check_index{
    my $file = shift;
    my @location = split(/\//, $file);
    my $rcf      = pop (@location);
    my $tbi      = $rcf.".tbi";
    my $flag     = 0;
    my $index    = join("/", @location, $tbi);
    return($index);
}

sub index_version{
    my $index     = shift;
    my $flag = "no_version";
    #verify existence of RCF
    if ((-e $index) and (-r $index))
    {
	#verify that the CGAPipeline used was version 2.* and higher
	my $line =  `zcat $index | head \-10 | grep 'CGAPipeline_'`;
	$line=~m/_(\d{1})/;
        $flag = $1;
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

