#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_refFlat\"   \"input_expression\"    \"output_bed\"    \"\(optional\)track_name\"   \"\(optional\)min_Splicing_Score_sum==10\"    \"\(optional\)min_expression==2\"    \"\(optional\)erase_CDS == 0\"  \"\(optional\)add_\"chr\" > 0\" " if (@ARGV < 2);
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $fileout=$ARGV[2];
my $track_name=$filein1;
if (scalar(@ARGV) > 3) { $track_name=$ARGV[3]; }
my $SSum=10;
if (scalar(@ARGV) > 4) { $SSum=$ARGV[4]; }
my $count=2;
if (scalar(@ARGV) > 5) { $count=$ARGV[5]; }
my $erase=1;
if (scalar(@ARGV) > 6) { $erase=$ARGV[6]; }
my $addchr=0;
if (scalar(@ARGV) > 7) { $addchr=$ARGV[7]; }
open(IN, $filein1) || die "Can't open $filein1 for reading!\n";
open(OUT, ">".$fileout) || die "Can't open $fileout for writing!\n";
#
# format examples: 
# I. refFlat   half-open zero-based. This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99.
# XLOC_005566 TCONS_00012033  chr6  + 171030441 171044890 171030441 171030441 3 171030441,171037018,171044838,  171030573,171037146,171044890,
# XLOC_005566 TCONS_00011428  chr6  + 171044838 171045633 171044838 171044838 2 171044838,171045435,  171044973,171045633,
# XLOC_006319 TCONS_00014242  chr7  + 158743206 158750045 158743206 158743206 3 158743206,158749273,158749907,  158743609,158749374,158750045,
#
# II. Bed:     half-open 1-based.
# chr1  11873 14409 uc001aaa.3  0 + 11873 11873 0 3 354,109,1189, 0,739,1347,
# chr1  11873 14409 uc010nxr.1  0 + 11873 11873 0 3 354,52,1189,  0,772,1347,
# chr1  11873 14409 uc010nxq.1  0 + 12189 13639 0 3 354,127,1007, 0,721,1529,

print OUT "track name=\"".$track_name."\" description=\"".$track_name."\" visibility=2 itemRgb=\"On\"\n";
my $score=0;
my $col=0;
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    my $name=$a[1];
    $name = $a[0].'|'.$a[1] unless($a[0] eq $a[1]);
    my $lengths="";
    my $starts="";
    my @Start=split(/\,/,$a[9]);
    my @End=split(/\,/,$a[10]);
    $a[4]--;
    $a[5]--;
    for (my $i=0; $i<$a[8]; $i++) {
        my $cur_s = $Start[$i] -1 - $a[4];
        my $cur_l = $End[$i] - $Start[$i];
        $lengths = $lengths.$cur_l.',';
        $starts = $starts.$cur_s.',';
    }
    if ($erase eq 0) {
        $a[6]=$a[4];
        $a[7]=$a[4];
    }
    else {
        $a[6]=$a[4];
        $a[7]=$a[5];
    }
    if ($addchr > 0) {
        $a[2]="chr".$a[2];
    }
    if ($a[2] eq "chrMT") { $a[2]="chrM"; }
    if (scalar(@a) eq 12) {
        my @b=split(/\,/,$a[11]);
        print OUT join("\t",$a[2],$a[4],$a[5],$name."|".scalar(@b),scalar(@b),$a[3],$a[6],$a[7],$col,$a[8],$lengths,$starts),"\n";
    }
    else {
        print OUT join("\t",$a[2],$a[4],$a[5],$a[1],$score,$a[3],$a[6],$a[7],$col,$a[8],$lengths,$starts),"\n";
    }
}
