#!/usr/bin/perl
use strict;
use warnings;


my $tabix = `which tabix`;
chomp($tabix);
die unless ($tabix);
open my $command_line, '>', "Command_file.txt";
my @chr = split(/,/, $ARGV[1]);

foreach my $chr (@chr) {
    print "chr$chr\n";
    my $chrom = "chr".$chr;
    my $new_file = join(".", "singletons", $chr, "temp");
    open my $nf, ">", $new_file;
    #grab regions
    my $old_pos = 0;
    my $pos; 
    open my $fh, "$tabix $ARGV[0] $chrom|" or die ("Could not open $ARGV[0]");
    my $n=0;
    while(<$fh>){
	next if (m/\#\#/);
        (undef, $pos) = split(/\t/);
	print $nf "$pos\n" unless (defined($pos) and ($pos eq $old_pos));
	$old_pos = $pos;
    }
    close($nf);

    my $log_file        = "random_log_file".$chr.".txt";
    my $chr_vector_file = $chrom.".vector_file";
    my $screen_name     = $chrom."mss";
    my $final_file      = "test.".$chrom.".file";

    system("nice \-5 perl create_vector.pl $new_file $chr_vector_file > $log_file");
    system("nice \-5 perl mssv2.pl $chr_vector_file > $final_file\; rm $chr_vector_file");
}
