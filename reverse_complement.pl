#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_fasta\"    \"output_reverse_complemented_fasta\"   " if (@ARGV < 1);
my $filein=$ARGV[0];
my $fileout=$filein.".rc";
if (scalar(@ARGV) > 1) { $fileout=$ARGV[1]; }
open(IN, $filein) or die "Cannot open input_fasta file : $filein";
open(OUT,">".$fileout);
while (<IN>) {
    chomp;
    if(m/^>/||m/^@/){
        s/^>//;
        print OUT ">rc".$_,"\n";
    }
    else{
        my $t=$_;
        $t=~tr/[ATCG]/[TAGC]/;
        my $r=scalar reverse $t;
        print OUT $r,"\n";
    }
}
