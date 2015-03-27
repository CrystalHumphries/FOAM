#!/usr/bin/perl
use strict;
use POSIX;

open my $fh, $ARGV[0];
my $x = "T";
chomp(my @array = <$fh>);

my $start_pos = $array [0];
my @result =  create_vector(@array);

if ($x eq "T"){
    print_header();
    my $outLR = mss(\@result);
    for my $lr (@{$outLR}) {
	my ($blk_start, $blk_end) = get_positions($lr->[0], $lr->[1]);
	print join("\t", $blk_start, $blk_end, $lr->[2])."\n";
    }
}


sub get_positions{
    my ($begin, $end) = @_;
    my $block_start   = $start_pos + $begin;
    my $block_end     = $start_pos + $end;
    return($block_start, $block_end);
}

sub print_header{
    my $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime;
    my $username = $ENV{LOGNAME} || $ENV{USER} || getpwuid($<);
    print "## $now_string\n";
    print "## $username\n";
    print join("\t", "Start", "Stop", "SegmentScore")."\n";
}

sub create_vector{
    my @positions = @_;
    my $start     = my $begin =  shift(@positions);
    my @new_vector;

    push (@new_vector, 2);
    foreach my $pos (@positions){
	my $diff = $pos - $begin;
	$diff--;
	push (@new_vector, 2);  # pushes a value of 2 for existing variant; this number can be changed

	#pushes a value of -5 for all the missing values
	push (@new_vector, ('-5') x $diff);
	$begin = $pos;
    }
    my $array_distance = $begin - $start;
    $array_distance++;
    print "Error in distance\n" unless ($array_distance==scalar(@new_vector));
    return (@new_vector);
}
		

######################################################################
### $blockLR = &mss($scoreLR);
###
### Ruzzo & Tompa's maximal scoring subsequences algorithm (whence mss),
### adapted from the maximal_seg() subroutine from BroadPeak (Wang,
### Luniak, and Jordan 2013, Bioinformatics 29(4):492-3.
###
### Given a list of scores, finds the non-overlapping blocks of
### positive total score that are "maximal" in the following sense:
###
### 1.  All blocks within a maximal block have lower score, and
### 
### 2.  Any block containing a maximal block contains some block
###with a score at least as high as its own.
###
### The input is an array of scores; additional information
### associated with the score should reference the score by its array
### index, and is not passed to this subroutine.
###
### The output is a list of blocks, which are "records"
### (array references) containing the range of (0-based) indices
### [left..right] (length 1+right-left) and the block's
### score.
###
### Example input: [0,0, 1, -1, 1,1, -1,-1,-1,-1,-1,-1,-1,-1, 1,1,1,0,0,1]
### Example output: [ [2,2,1], [4,5,2], [14,19,4] ]
###
### Notes:
### 1.  This algorithm finds clusters of _positive_ scores.  Any region with
###an average score of just greater than zero will cause the cluster to be
###extended; so an average score of zero is the "clustering threshold."
###
###The theory of scoring systems with this property was presented
###in "Methods for assessing the statistical significance of molecular
###sequence features by using general scoring schemes."
###Proc Natl Acad Sci U S A. 1990 Mar;87(6):2264-8.  Karlin S, Altschul SF.
###
###The Karlin-Altschul theory is one of the foundations of biological
###sequence analysis, and IMHO every bioinformatician should understand it.
###
### 2.  For use in clustering genomic features, it is likely that the
###features to be clustered will have been identified and assigned
###positive scores by some other method, which generates a sparse list
###of the sites (or regions) by their genomic coordinates along with
###a score for each.  Each item on such a list should be represented
###in the input to this algorithm as a single positive score entry,
###and the uninteresting "background" genomic sequence between two
###consecutive features should be represented as a single total negative
###score for the region.
###
###Even for "dense" scores like CADD, where every nucleotide is assigned
###a score, it is useful to preprocess the scores into runs of scores <= 0
###and runs of scores > 0 (every cluster found by this algorithm will
###ALWAYS start and end with a positive score, never a zero).  If the
###dense scores are all positive, you will need to make some of them
###negative by subtracting some score threshold T.
###
### 3.The essential question for using this algorithm is, what should be
###the value of the negative background scores for the intervening
###regions?  Unless the degree of clustering should depend on features of
###the surrounding genomic sequence (e.g. GC content, proximity to a gene,
###etc.), it is sufficient to choose a single per-nucleotide negative
###"background" score -B, and score a background region of n nucleotides
###as -B*n.
###
### 4.  One of the points made by Karlin and Altschul is that unless
###the expected score of a randomly-chosen region containing both
###positive and negative segments is negative, any cluster of
###positive scores is expected to be extended even at random, and
###the clustering will be pretty meaningless.  Thus, if the total
###score of the pre-identified features is P and the total
###number of bases in the genome outside the pre-identified
###features is N, -B should be chosen to at least be less than -P/N.
###Karlin and Altschul give additional information about what it
###means to choose -B to be progressively more and more negative
###than -P/N; basically, the larger the magnitude of B (the more
###negative -B is), the more stringent the clustering and the
###smaller the clusters of features will be.
######################################################################
sub mss {
    my ($scoreLR) = @_;
    my $scores  = scalar @{$scoreLR};
    my($blocks, $sLeft, @left, @right, @sLeft, @sRight);
    $blocks = 0;
    for (my $s = 0; $s < $scores; $s++) {
        my $score = $scoreLR->[$s];
        if (0 < $score) {
	    $left[$blocks]= $s;
	    $right[$blocks]= $s;
	    $sLeft[$blocks]= $sLeft;
	    $sRight[$blocks]= my $sRight = $sLeft + $score;
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

    
    
	    
