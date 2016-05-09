#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"List\"  \"Database_to_screen_from\"  \"Output\"  \(select==1_or_-1\)" if (@ARGV < 3);
my $filein=$ARGV[0];    # empty1
my $DB=$ARGV[1];        # 454Isotigs.fna.fa
my $fileout=$ARGV[2];   # redo_empty
my $select=1;
if (scalar(@ARGV) > 3) {$select=$ARGV[3];}
if (($select ne 1) and ($select ne -1)) {
    die "select can be either 1 (select) or -1 (unselect) \n";
}

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
    if ($select eq 1) {
        if (exists $uniq{$a[0]}) {
            print OUT1 ">".$_,"\n";
            my $seq=<IN1>;
            print OUT1 $seq;
        }
    }
    else {
        if (!exists $uniq{$a[0]}) {
            print OUT1 ">".$_,"\n";
            my $seq=<IN1>;
            print OUT1 $seq;
        }
    }
    
}
