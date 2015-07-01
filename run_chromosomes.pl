#!/usr/bin/perl
use strict;
use warnings;

my $loop_no = 0;
my $old_file = "NULL";
my $n = 0;


LOOP: for my $chr (@ARGV){

    my $chromosome="chr".$chr;
    
    system("nice \-5 perl findFrequentSimplev2_CGI.pl /isb/chumphri/FOAM/Scripts/CGI_founder_list_20150120.txt $chr");

    my $file="fulloutdir1/raw.".$chromosome.".gz";
    
    if ($old_file eq $file){
	$loop_no++;
    }else{
	$loop_no=1;
    }

    $old_file = $file;
    if ( ( ! -e $file) or (!-z $file) ){

	my $var= `zcat $file |  awk 'BEGIN {max = 0} {if (\$4>max) max=\$4} END {print max}' `;

	if ( $var == 1669 ){
	    print "you're done with ".$chromosome."\n";
	}else{
	    system ("echo \"Running $chromosome for the $loop_no time\" | mail -s \"rerun $chromosome\" crystal.humphries\@systemsbiology.org");
	    print "finished with loop 3 of $chromosome\n" and next LOOP if ($loop_no >= 3);
	    redo LOOP;
	}
    }
    else{
	print "$file is empty.\n";
    }
}

LOOP2: for my $chr (@ARGV){
    my $chromosome="chr".$chr;
    my $file="fulloutdir1/raw.".$chromosome.".gz";
    my $new_file="fulloutdir1/".$chromosome."CADD.gz";
    system("nice -5 zcat $file | perl OverlapwithCADD.pl | bgzip -c > $new_file");
}
