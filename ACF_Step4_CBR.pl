#!/usr/bin/perl -w
use strict;
# convert CBR structure file into agtf format 
die "Usage: $0  \"circRNA\"  \"split_exon_gtf\"   \"output\" " if (@ARGV < 3);
my $filein=$ARGV[0];    
my $gtf=$ARGV[1];       
my $fileout=$ARGV[2];   
my $command="rm -f Step4_CBR_finished";
system($command);

my %anno;
if ($gtf ne "no") {
    open IN,$gtf;
    while(<IN>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[2] eq "exon") {
            if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
            if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
            my @b=split(/\_\_\_/,$a[8]);
            $anno{$a[0]}{$a[6]}{$a[3]."\t".$a[4]}=join("\t",@a);
        }
    }
    close IN;
}

sub getname ($) {
    my $Rank=shift;
    my $fill=6 - length($Rank);
    my $lead="";
    for(my $i=1; $i<=$fill; $i++) { $lead=$lead."0"; }
    $fill=$lead.$Rank;
    return $fill;
}

my $rank=0;
open IN,$filein;
open OUT,">".$fileout.".agtf";
open OUTrefFlat,">".$fileout.".refFlat";
open OUTg,">".$fileout.".agtf_gene";
while(<IN>) {
    chomp;
    if (m/^#/) { next; }
    my @a=split("\t",$_);
    if ($a[1]=~m/^chromosome/i) {$a[1]=~s/chromosome//i;}
    if ($a[1]=~m/^chr/i) {$a[1]=~s/chr//i;}
    my $left=$a[2] < $a[3] ? $a[2] : $a[3];
    my $right=$a[2] > $a[3] ? $a[2] : $a[3];
    $rank++;
    if (exists $anno{$a[1]}{$a[20]}) {
        my %record;
        my %Gname;
        foreach my $info (keys %{$anno{$a[1]}{$a[20]}}) {
            my @b=split("\t",$anno{$a[1]}{$a[20]}{$info});
            if (($left <= $b[3]) and ($b[4] <= $right)) {
                $Gname{$b[7]}++;
                $record{$b[3]}=$b[4];
            }
        }
        foreach my $pos (sort{$a <=> $b} keys %record) {
            my $tmp=$record{$pos};
            delete $record{$pos};
            $record{$left}=$tmp;
            last;
        }
        foreach my $pos (sort{$b <=> $a} keys %record) {
            $record{$pos}=$right;
            last;
        }
        my $gname="na";
        foreach my $id (sort{$Gname{$b} <=> $Gname{$a}} keys %Gname) { $gname=$id; last; }
        my @ExonL;
        my @ExonR;
        my $Nr=0;
        foreach my $pos (keys %record) { $Nr++; }
        if ($Nr > 0) {
            if ($a[20] eq "+") {
                my $count=1;
                foreach my $pos (sort {$a <=> $b} keys %record) {
                    print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gname,$a[0]."___".$count."___".$Nr),"\n";
                    #print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gname,"CBR_".getname($rank)."___".$count."___".$Nr),"\n";
                    $ExonL[$count-1]=$pos;
                    $ExonR[$count-1]=$record{$pos}+1;
                    $count++;
                }
            }
            else {
                my $count=$Nr;
                foreach my $pos (sort {$a <=> $b} keys %record) {
                    print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gname,$a[0]."___".$count."___".$Nr),"\n";
                    #print OUT join("\t",$a[1],"split","exon",$pos,$record{$pos},"na",$a[20],$gname,"CBR_".getname($rank)."___".$count."___".$Nr),"\n";
                    $ExonL[$Nr - $count]=$pos;
                    $ExonR[$Nr - $count]=$record{$pos}+1;
                    $count--;
                }
            }
            print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],$gname,$a[0]),"\n";
            #print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],$gname,"CBR_".getname($rank)),"\n";
            print OUTrefFlat join("\t",$gname,$a[0],"chr".$a[1],$a[20],$left,$right+1,$left,$left,($Nr),join(",",@ExonL),join(",",@ExonR),$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
        }
        else {
            print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],"na",$a[0]."___1___1"),"\n";
            print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],"na",$a[0]),"\n";
            #print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],"na","CBR_".getname($rank)."___1___1"),"\n";
            #print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],"na","CBR_".getname($rank)),"\n";
            print OUTrefFlat join("\t","na",$a[0],"chr".$a[1],$a[20],$left,$right+1,$left,$left,1,$left,$right+1,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
        }
    }
    else {
        print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],"na",$a[0]."___1___1"),"\n";
        print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],"na",$a[0]),"\n";
        #print OUT join("\t",$a[1],"split","exon",$left,$right,"na",$a[20],"na","CBR_".getname($rank)."___1___1"),"\n";
        #print OUTg join("\t",$a[1],"split","gene",$left,$right,"na",$a[20],"na","CBR_".getname($rank)),"\n";
        print OUTrefFlat join("\t","na",$a[0],"chr".$a[1],$a[20],$left,$right+1,$left,$left,1,$left,$right+1,$a[17],$a[18],$a[19],$a[21],$a[22],$a[23]),"\n";
    }
}



open(OUTFLAG,">Step4_CBR_finished");
print OUTFLAG "Step4_CBR_finished\n";
close OUTFLAG;

