#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_raw_fasta\"   \"output\"  \"info_to_add_to_header\"  \"Sequencing Mode \(SE==0 or PE==2\)\"  \"\(optional\) remove_sequences_with_N == 1\"" if (@ARGV < 3);
my $filein1=$ARGV[0];
my $fileout=$ARGV[1];
my $info=$ARGV[2];
my $SM=0;
if (scalar(@ARGV) > 3) { $SM=$ARGV[3]; }
if (($SM ne 0) and ($SM ne 2)) { die "Sequencing Mode can only be : 0 or 2!"; }
my $sanity=1;
if (scalar(@ARGV) > 3) { $sanity=$ARGV[4]; }
if (($sanity ne 0) and ($sanity ne 1)) { die "remove_sequences_with_N==1 or keep_sequences_with_N==0 !"; }
open(IN, $filein1) or die "Cannot open input_raw_fasta $filein1 ";
open(OUT, ">".$fileout) or die "Cannot open output file $fileout ";
while (<IN>) {
    chomp;
    if (m/^>/) {
        s/^>//;
        my @a=split(" ",$_);
        my $seq=<IN>;
        if (($sanity eq 1) and ($seq=~m/N/)) { next; }
        if ($SM eq 0) {
            print OUT ">".$info."_".$a[0],"\n",$seq;
        }
        else {
            print OUT ">".$info."_".$a[0]."/".$SM,"\n",$seq;
        }
    }
    elsif (m/^@/) {
        s/^@//;
        my @a=split(" ",$_);
        my $seq=<IN>;
        my $qua=<IN>;
        $qua=<IN>;
        if (($sanity eq 1) and ($seq=~m/N/)) { next; }
        if ($SM eq 0) {
            print OUT "@".$info."_".$a[0],"\n",$seq,"+\n",$qua;
        }
        else {
            print OUT "@".$info."_".$a[0]."/".$SM,"\n",$seq,"+\n",$qua;
        }
    }
    
}
close IN;
close OUT;
