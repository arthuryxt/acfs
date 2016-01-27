#!/usr/bin/perl -w
use strict;
my $filein=$ARGV[0];    # list
my $DB=$ARGV[1];        # big.fa
my $fileout=$ARGV[2];   # output
open IN,$filein;
my %uniq;
while(<IN>) {
    chomp $_;
    #my @a=split(/\_/,$_);
    $uniq{$_}=1;
}
close IN;
open IN1,$DB;
open OUT1,">".$fileout;
while(<IN1>) {
    chomp $_;
    next unless m/^>/;
    s/>//;
    my @a=split("\t",$_);
    if (exists $uniq{$a[0]}) {
        print OUT1 ">".$_,"\n";
        my $seq=<IN1>;
        print OUT1 $seq;
    }
}
