#!/bin/env perl
use POSIX;
use Tie::File; 
use Time::HiRes qw( gettimeofday tv_interval );


#written by Crystal Humphries and adapted library MSS
#the library MSS was written by Max Robinson and Nigel Clegg
my $file = create_temp($ARGV[0]);

open my $fh, $file;
my @vector;
chomp(my @vector = <$fh>);

print_header( scalar(@vector));
my $start_pos = $result[0];

$outLR = &mss(\@vector);
undef(@vector);

my $loc_file = create_locations($ARGV[0]);
open my $fh, $loc_file;
my @locations;
chomp(my @locations = <$fh>);

for my $lr (@{$outLR}) {
    my $start     = $locations[ $lr->[0]];
    my $end_pos   = $lr->[1] + 1;
    my $end       = $locations[ $end_pos] - 1 ;
    my $bin_size  = $end - $start;
    print join("\t", $start, $end, $bin_size, $lr->[2])."\n";
}
my $count = scalar (@{$outLR});
print $count."\n";

system ("echo \"Finished running $0\" | mail -s \"$0 is completed\" crystal.humphries\@systemsbiology.org");
system("rm -f $file $loc_file");

sub create_temp{
    my $file = shift;
    my $temp = my $orig = $file.".trunc";
    my $n=0;
    unless (!-e $temp){
	$n++;
	$temp = join(".",$orig, $n);
    }
    system ("cut \-f1 $ARGV[0] > $temp");
    return ($temp);
}

sub print_header{
    my $vector_size = shift;
    my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
    my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
    print "## $now_string\n";
    print "## $username\n";
    print "## $ARGV[0]: $vector_size length\n";
    print join("\t", "Start", "Stop", "SegScore", "BinSize")."\n";
}

sub create_locations{
    my $file = shift;
    my $temp = my $orig = $file.".loc";
    my $n=0;
    unless (!-e $temp){
        $n++;
        $temp = join(".",$orig, $n);
    }
    system ("cut \-f2 $ARGV[0] > $temp");
    return ($temp);
}

sub mss {
    my $t0      = [gettimeofday()];
    my ($scoreLR) = @_;
    my $scores  = scalar @{$scoreLR};
    my($blocks, $sLeft, @left, @right, @sLeft, @sRight);
    for (my $s = 0; $s < $scores; $s++) {
        my $score = $scoreLR->[$s];
        if (0 < $score) {
        	$left[$blocks]	= $s;
        	$right[$blocks]	= $s;
        	$sLeft[$blocks]	= $sLeft;
        	$sRight[$blocks]	= my $sRight = $sLeft + $score;
        	my $flag = 1;
        	while ($flag) { # until we stop finding blocks to extend
        		my $j = $blocks-1;
        		until (($sLeft[$j] < $sLeft[$blocks]) || ($j < 0)) {
        			$j--;
        		}
        		if ($j < 0) {
        			$flag = 0;
        		} else {
        			my $j1 = $j;
        			if ($sRight[$j1] >= $sRight[$blocks]) {
        				$flag = 0;
        			} else {
        				$right[$j1]  = $s;
        				$sRight[$j1] = $sRight;
        				$blocks = $j1;
        			}
        		}
        	}
        	$blocks++;
        }
        $sLeft += $score;
    }
    my $blockLR;
    for (my $b = 0; $b < $blocks; ++ $b) {
        push @{$blockLR}, [$left[$b], $right[$b], $sRight[$b]-$sLeft[$b]];
    }
    return $blockLR;
}


