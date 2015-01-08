#!/usr/bin/env perl
$|=1;
use strict;

my $indCutoff = 15;
my $distCutoff = 200;
my $minLen = 10;
my $base = 1.5;
my $padding = 50;
my $logbase = log($base);
my @hist;
my($prevChrom, $prevPos, $start, $end, %sites, $events, %cadd);
my $version = "1";
my $counter;


while (<>) {
	chomp;
	my($chrom, $pos, $var, $ind, $inds, $FAVA, $cadd_raw, $cadd_scale) = split (/\t/);
	#my ($cadd_raw, $cadd_scale) = check_CADD($cadd_1, $cadd_2); 
	$cadd_raw=~s/NA/0/;
	$cadd_scale=~s/NA/0/;
	if ($chrom ne $prevChrom || $pos-$prevPos>=$distCutoff) {
		output($prevChrom, $start, $end, $events, \%sites, \%cadd) if $end-$start>=$minLen;
		$prevChrom = $chrom;
		$start = $end = $pos;
		$events = 0;
		%sites = ();
		%cadd = ();
	}

	$prevPos = $end = $pos;
	$sites{$pos}++;
	$cadd{$pos}+=$cadd_scale;
	$events += $ind; #assuming independent individuals, which isn't the case...
}
output($prevChrom, $start, $end, $events, \%sites) if $end-$start>=$minLen;


sub output {
    my $sum = 0;
    my($chrom, $start, $end, $events, $sitesref, $cadd) = @_;
    my $sites = scalar keys %{$sitesref};
    my $len = $end-$start+1+$distCutoff;
    my $score = $sites/($len);
    my @cadd = values %{$cadd};
    $sum += $_ foreach @cadd;
    my $avg_cadd = $sum/$sites;
    return unless ($score>=0.02);
    $counter++;
    my $id = join(".", "cms", $version, $chrom, $counter);
    print join("\t", "chr$prevChrom", $start-$padding, $end+$padding, "cms:$id", sprintf("%.2f", $score), $len, $events, $sites, join(",", sort {$a<=>$b} keys %{$sitesref}), $sum, sprintf("%.2f", $avg_cadd)), "\n";
}


sub check_CADD{
    my (@new);
    foreach my $score (@_){
	if ($score == "NA"){
	    push(@new, 0);
	}else{
	    push(@new, $score);
	}	    
    }
    return(@new);
}
    
