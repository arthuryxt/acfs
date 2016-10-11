#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"output_basename\"   \"input_sam\"   \"\(optional\)remove_N==1\"  \"\(optional\)separate_output_every_line==10000\"  " if (@ARGV < 2);
my $fileout=$ARGV[0];
my $filein=$ARGV[1];    # sam must be ordered by id.
my $removeN=1;
if (scalar(@ARGV) > 2) { $removeN=$ARGV[2]; }
my $sep=10000;
if (scalar(@ARGV) > 3) { $sep=$ARGV[3]; }

my $rcnt=0;
my $lastid="";
my $fcnt=0;

open(IN, $filein) or die "Cannot open input sam : $filein\n";
open(OUT, ">".$fileout."_".$fcnt);
while (<IN>) {
    chomp;
    if (m/^@/) { next; }
    my @a=split("\t",$_);
    if (($rcnt > $sep) and ($lastid ne $a[0])) {
        $fcnt++;
        $rcnt=0;
        close OUT;
        open(OUT, ">".$fileout."_".$fcnt);
        if (($removeN eq 1) and ($a[9] !~m/N/)){
            $lastid=$a[0];
            $rcnt++;
            print OUT join("\t",@a),"\n";
        }
    }
    else {
        if (($removeN eq 1) and ($a[9] !~m/N/)){
            $lastid=$a[0];
            $rcnt++;
            print OUT join("\t",@a),"\n";
        }
    }
}




