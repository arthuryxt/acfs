#!/usr/bin/perl -w
use strict;
my $filein=$ARGV[0];    # interesting list
my $DB=$ARGV[1];        # database file, all entry
my $fileout=$ARGV[2];   # resulting interesting entry
my $pos=0;              # 0-based column position
if ($ARGV[3]) {$pos=$ARGV[3]}
open IN,$filein;
my %uniq;
while(<IN>) {
    chomp $_;
    $uniq{$_}=1;
}
close IN;
open IN1,$DB;
open OUT1,">".$fileout;
while(<IN1>) {
    chomp $_;
    my @a=split("\t",$_);
    if (exists $uniq{$a[$pos]}) {
        print OUT1 $_,"\n";
    }
}
