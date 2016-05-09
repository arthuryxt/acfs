#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"List\"  \"Database_to_screen_from\"  \"Output\"  \"zero-based_column\"  \(select==1_or_-1\)" if (@ARGV < 3);
my $filein=$ARGV[0];    # interesting list
my $DB=$ARGV[1];        # database file, all entry
my $fileout=$ARGV[2];   # resulting interesting entry
my $pos=0;              # 0-based column position
if (scalar(@ARGV) > 3) {$pos=$ARGV[3];}
my $select=1;
if (scalar(@ARGV) > 4) {$select=$ARGV[4];}
if (($select ne 1) and ($select ne -1)) {
    die "select can be either 1 (select) or -1 (unselect) \n";
}

open IN,$filein;
my %uniq;
while(<IN>) {
    chomp $_;
    $uniq{$_}=1;
}
close IN;
open IN1,$DB;
open OUT1,">".$fileout;
if ($select eq 1) {
    while(<IN1>) {
        chomp $_;
        my @a=split("\t",$_);
        if (exists $uniq{$a[$pos]}) {
            print OUT1 $_,"\n";
        }
    }
}
else {
    while(<IN1>) {
        chomp $_;
        my @a=split("\t",$_);
        if (exists $uniq{$a[$pos]}) {}
        else {
            print OUT1 $_,"\n";
        }
    }
}

