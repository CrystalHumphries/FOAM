#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(sum);


my $file = $ARGV[0];
print $file."\n";

my $tabix = `which tabix`;
chomp($tabix);

open TBX,"$tabix $file chr22:0-90000000 |  awk '{ if (\$5!~/\\./) print}' |" || die "Could not open $file\n";
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
print "Sum: $sum\n";
my $avg_dist = $sum / $n;
print "Distance: $avg_dist \n";