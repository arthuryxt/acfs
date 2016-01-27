#!/usr/bin/perl -w
use strict;
# convert circRNA structure file into raw_gtf format 
die "Usage: $0  \"circRNA\"  \"split_exon_gtf\"   \"output_basename\"  \"\(optional\) extend N bases\"  " if (@ARGV < 3);
my $filein=$ARGV[0];    
my $gtf=$ARGV[1];       
my $fileout=$ARGV[2];   
my $Extend=50;         # 100nt by default.
if (scalar(@ARGV) > 3) {$Extend=$ARGV[3];}
my $command="rm -f Step4_MEA_finished";
system($command);

my %anno;
my %ExonNr;
open IN,$gtf;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[2] eq "exon") {
        if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
        if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
        my @b=split(/\_\_\_/,$a[8]);
        #     ENSG   exon_number
        $anno{$b[0]}{$b[1]}=join("\t",@a);
        $ExonNr{$b[0]}=$b[2];
    }
}
close IN;

open IN,$filein;
open OUT,">".$fileout.".gtf";               # close zero-based. This means that the first 100 bases of a chromosome are represented as [0,99]
open OUTrefFlat,">".$fileout.".refFlat";    # half-open zero-based. This means that the first 100 bases of a chromosome are represented as [0,100), i.e. 0-99. The second 100 bases are represented as [100,200), i.e. 100-199.
open OUT1,">".$filein.".ext".$Extend;
open OUT11,">".$fileout.".gtf.ext".$Extend;
open OUT2,">".$fileout.".gtf_2G";
open OUT3,">".$fileout.".err";
while(<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
    if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
    my @left=split(/\_\_\_/,$a[5]);
    my @right=split(/\_\_\_/,$a[11]);
    if ($left[0] ne $right[0]) { print OUT2 join("\t",@a),"\n"; next;}
    my $start=$left[1] < $right[1] ? $left[1] : $right[1];
    my $end=$left[1] > $right[1] ? $left[1] : $right[1];
    my $Nr=$end - $start +1;
    my @ExonL;
    my @ExonR;
    my $strand=$a[7];
    if ($a[7] ne $a[13]) { print OUT3 "discrepancy in two-part-annotation\t",join("\t",@a),"\n"; next; }
    if ($a[7] ne $a[20]) { print OUT3 "discrepancy in annotation and read strandness\t",join("\t",@a),"\n"; next; }
    if ($a[2] ne $a[10]) { print OUT3 "discrepancy in left-border\t",join("\t",@a),"\n"; }
    if ($a[3] ne $a[15]) { print OUT3 "discrepancy in right-border\t",join("\t",@a),"\n"; }
    for(my $i=1; $i<=$Nr; $i++) {
        my @b=split("\t",$anno{$left[0]}{$i+$start-1});
        my $info=join(" ","gene_id","\"$a[0]\"\;","transcript_id","\"$a[0]\"\;","gene_name","\"$a[6]\"\;","gene_name2","\"$left[0]\"\;","exon_number","\"$i\"");
        my $tmp_s=$b[3];
        my $tmp_e=$b[4];
        if (($tmp_s <= $a[3]) and ($a[3] <= $tmp_e)) { $tmp_s=$a[3]; }
        elsif (($tmp_s <= $a[2]) and ($a[2] <= $tmp_e)) { $tmp_e=$a[2]; }
        print OUT join("\t",$b[0],$b[1],$b[2],$tmp_s,$tmp_e,$b[5],$b[6],$b[7],$info),"\n";
        if ($strand eq "+") {
            $ExonL[$i-1]=$tmp_s;
            $ExonR[$i-1]=$tmp_e+1;
        }
        else {
            $ExonL[$Nr-$i]=$tmp_s;
            $ExonR[$Nr-$i]=$tmp_e+1;
        }
    }
    if (($ExonL[0] ne $a[3]) or ($ExonR[$Nr-1] ne ($a[2]+1))) {
        print OUT3 "fix-border\t",join("\t",@a),"\n";
        $ExonL[0]=$a[3];
        $ExonR[$Nr-1]=$a[2]+1;
    }
    print OUTrefFlat join("\t",$a[6],$a[0],"chr".$a[1],$strand,$a[3],$a[2]+1,$a[3],$a[3],$Nr,join(",",@ExonL),join(",",@ExonR),$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
    my $cumlen=0;
    # check for extension on left side
    my $start_move=0;
    for (my $i=$start-1; $i>0; $i--) {
        my @tmp=split("\t",$anno{$left[0]}{$i});
        my $tmplen=$tmp[4] - $tmp[3];
        $cumlen+=$tmplen;
        if ($cumlen < $Extend) {$start_move++;}
        else {last;}
    }
    # check for extension on right side
    my $end_move=0;
    for (my $i=$end+1; $i<=$ExonNr{$left[0]}; $i++) {
        my @tmp=split("\t",$anno{$left[0]}{$i});
        my $tmplen=$tmp[4] - $tmp[3];
        $cumlen+=$tmplen;
        if ($cumlen < $Extend) {$end_move++;}
        else {last;}
    }
    for (my $i=0; $i<=$start_move; $i++) {
        for (my $j=0; $j<=$end_move; $j++) {
            # output all possible combinations, except the original ones
            if (($i eq 0) and ($j eq 0)) {next;}
            my $nstart=$start-$i;
            my $nend=$end+$j;
            my @Anstart=split("\t",$anno{$left[0]}{$nstart});
            my @Anend=split("\t",$anno{$left[0]}{$nend});
            my $nleftmost=$Anstart[3] < $Anend[3] ? $Anstart[3] : $Anend[3];
            my $nrightmost=$Anstart[4] > $Anend[4] ? $Anstart[4] : $Anend[4];
            my $nid="C".$a[0]."_".$i."_".$j;
            print OUT1 join("\t",$nid,$a[1],$nleftmost,$nrightmost,($nrightmost - $nleftmost),$Anstart[8],$Anstart[7],$Anstart[6],$Anstart[5],$Anstart[3],$Anstart[4],$Anend[8],$Anend[7],$Anend[6],$Anend[5],$Anend[3],$Anend[4]),"\n";
            
            my $Nr=$nend - $nstart +1;
            for(my $k=1; $k<=$Nr; $k++) {
                my @b=split("\t",$anno{$left[0]}{$k+$nstart-1});
                my $info=join(" ","gene_id","\"$a[0]\"\;","transcript_id","\"$nid\"\;","gene_name","\"$a[6]\"\;","gene_name2","\"$left[0]\"\;","exon_number","\"$k\"");
                print OUT11 join("\t",$b[0],$b[1],$b[2],$b[3],$b[4],$b[5],$b[6],$b[7],$info),"\n";
            }
        }
    }
    
}



open(OUTFLAG,">Step4_MEA_finished");
print OUTFLAG "Step4_MEA_finished\n";
close OUTFLAG;

