#!/usr/bin/env perl
$|=1;
use strict;

print_warning() and die  unless ((scalar(@ARGV)==2) and (-e $ARGV[0]) and ($ARGV[1]=~m/wgs|index/i));

#grabs list of sample names
my $wanted_samples = grab_list($ARGV[0]); 

#grabs only the sample ids that have a genome.rcf.gz file
my $fgs = "/40TB_3/gestalt/private/IndividualGenomeAssembly/";
my($sp, $ind_dir) = identifyGenomes_list($wanted_samples, $fgs);

my @ind = sort keys %{$sp}; #individual ids 

my $count = scalar(@ind);
print "Number of Individual Files: $count\n";

foreach my $ind (@ind){
    my $file  = "$ind_dir->{$ind}Genotype.rcf.gz";
    my $e_gen = my $r_gen = my $e_index = my $r_index = "noAccess";
    $e_gen = "FileExists"   if (-e $file);
    $r_gen =  "CanRead"      if (-r $file);
    my $index = check_index($file);
    $e_index = "FileExists" if (-e $index);
    $r_index = "CanRead"    if (-r $index);
    print join("\t", $file, $e_gen, $r_gen)."\n"  if ($ARGV[1]=~m/wgs/i);
    print join("\t", $index, $e_index, $r_index)."\n"  if ($ARGV[1]=~m/index/i);
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

sub print_warning{
    print<<EOF;

    print WGS_file_perm_check.pl <file_with genomes> <wgs|index|wgsindex>

EOF
}

sub grab_list{
    # grabs pipeline version of CGI genome
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

    #grabs samples from directory and checks whether file exists
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

