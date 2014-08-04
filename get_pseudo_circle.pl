#!/usr/bin/perl -w
use strict;
# concatenate the first N-bases to the end of the sequence, thus creating a pseudo-circle
die "Usage: $0  \"input_fasta\"    \"output\"  \"\(optional\) extend_N_bases\"  " if (@ARGV < 2);
my $filein=$ARGV[0];    
my $fileout=$ARGV[1];
my $extend=150;
if (scalar(@ARGV) > 2) { $extend=$ARGV[2]; }
if ($extend < 0) { die "extend_N_bases must be larger than 0 and optimal to be the sequencing length"; }
open(IN, $filein) or die "Cannot open input_fasta : $filein";
open(OUT, ">".$fileout) or die "Cannot open output file : $fileout";
while(<IN>){
    chomp;
    my $id=$_;
    my $seq=<IN>;
    chomp $seq;
    print OUT $id,"\n";
    print OUT $seq.substr($seq,0,$extend),"\n";
}



