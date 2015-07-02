#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);

my $tabix = `which tabix`;
chomp($tabix);

sub get_average_dist{
    my ($file, $query) = @_;
    open TBX,"$tabix $file $query |  awk '{ if (\$5!~/\\./) print}' |" || die "Could not open $file\n";
    my ($oldPos);
    my @distance;
    my $n = 0;

    while(my $line = <TBX>){
	chomp ($line);
	my (undef, $pos) = split(/\t/, $line);
	$n++;
	if ($n >1){
	    my $dist = $pos - $oldPos;
	    push (@distance, $dist);
	}
	$oldPos = $pos;
    }

    my $sum = sum (@distance);
    my $avg_dist = $sum / $n;
    
    return($sum, $avg_dist);
}


open my $fh, $ARGV[0] || die "cannot open $ARGV[0]";
my $chr = $ARGV[1];
my $query = $chr.":0-9000000000";

#create file to write results to
my $new_file = $chr."_distance_avg.txt";
open NFH, '>', $new_file;

while (my $WGS = <$fh>) {
    chomp $WGS;
    my @dir_struct = split(/\//, $WGS); 
    my $main_file  = $dir_struct[-1];
    $main_file     =~s/\.genome\.vcf\.gz//g;
    my ($sum, $avg_d) = get_average_dist($WGS, $query);
    print NFH join("\t", $main_file, $sum, $avg_d)."\n";
}

close(NFH);
