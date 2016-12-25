#!/usr/bin/perl -w
use strict;
die "Usage: $0   \"input_refFlat\"   \"input_expression\"    \"output_bed\"    \"\(optional\)track_name\"   \"\(optional\)min_Splicing_Score_sum==10\"   \"\(optional\)min_sample==1\"   \"\(optional\)min_expression==2\"    \"\(optional\)erase_CDS == 0\"  \"\(optional\)add_\"chr\" > 0\" " if (@ARGV < 2);
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $fileout=$ARGV[2];
my $track_name=$filein1;
if (scalar(@ARGV) > 3) { $track_name=$ARGV[3]; }
my $SSum=10;
if (scalar(@ARGV) > 4) { $SSum=$ARGV[4]; }
my $min_sample=1;
if (scalar(@ARGV) > 5) { $min_sample=$ARGV[5]; }
my $min_count=2;
if (scalar(@ARGV) > 6) { $min_count=$ARGV[6]; }
my $erase=1;
if (scalar(@ARGV) > 7) { $erase=$ARGV[7]; }
my $addchr=0;
if (scalar(@ARGV) > 8) { $addchr=$ARGV[8]; }
open(IN1, $filein1) || die "Can't open $filein1 for reading!\n";

open(OUT, ">".$fileout) || die "Can't open $fileout for writing!\n";

open(IN2, $filein2) || die "Can't open $filein2 for reading!\n";
my %Anno;
my %AnnoS;
while (<IN2>) {
    next if (m/^@/);
    next if (m/^#/);
    chomp;
    my @a=split("\t",$_);
    my $Nr=scalar(@a);
    if ($Nr <= 2) { next; }
    my @b=split(/\_\_\_/,$a[0]);
    my $skip=0;
    for(my $i=2; $i<$Nr; $i++) { if ($a[$i]=~m/\D/) { $skip++ }  }
    if ($skip > 0) { next; }
    #if ($a[2]=~m/\D/) { next; }
    my $cnt1=0;
    my $cnt2=0;
    for(my $i=2; $i<$Nr; $i++) {
        $cnt1+=$a[$i];
        if ($a[$i] > 0) { $cnt2++; }
    }
    $Anno{$b[0]}=$cnt1;
    $AnnoS{$b[0]}=$cnt2;
}
close IN2;

print OUT "track name=\"".$track_name."\" description=\"".$track_name."\" visibility=2 itemRgb=\"On\"\n";
my $score=0;
my $col=0;
while (<IN1>) {
    chomp;
    my @a=split("\t",$_);
    $score=0;
    if (($a[11] >= $SSum) and (exists $Anno{$a[1]}) and ($Anno{$a[1]} >= $min_count) and (exists $AnnoS{$a[1]}) and ($AnnoS{$a[1]} >= $min_sample)) {
        if ($a[11] <= 1) { $score=0; }
        else { $score=log($a[11]); }
        $score+=log($AnnoS{$a[1]});
        $score+=log($Anno{$a[1]});
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
        if ($addchr > 0) { $a[2]="chr".$a[2]; }
        if ($a[2] eq "chrMT") { $a[2]="chrM"; }
        my $name=$a[1]."|".$a[11]."|".$AnnoS{$a[1]}."|".$Anno{$a[1]};
        print OUT join("\t",$a[2],$a[4],$a[5],$name,sprintf("%.3f",$score),$a[3],$a[6],$a[7],$col,$a[8],$lengths,$starts),"\n";
    }
}
