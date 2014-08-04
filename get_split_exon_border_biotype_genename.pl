#!/usr/bin/perl -w
use strict;
die "Usage: $0 \"input\"  \"output\"  \"\(optional\)debug\" " if (@ARGV < 2);
my $filein=$ARGV[0];    
my $fileout=$ARGV[1];   
my $debug=0;
if (scalar(@ARGV) > 2) {$debug=$ARGV[2];}
open IN,$filein;
open OUT,">".$fileout;
open OUT2,">".$fileout."_gene";
my %uniq;
my %biotype;
my %Gname;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    # MT      protein_coding  exon    3307    4262    .       +       .        gene_id "ENSG00000198888"; transcript_id "ENST00000361390"; exon_number "1"; gene_name "MT-ND1"; gene_biotype "protein_coding"; transcript_name "MT-ND1-201";
    if (($a[2] eq "exon") and ($a[0]!~m/NT/)){
        my @b=split(/\"/,$a[8]);
        my $Nr=scalar(@b);
        my $gene_id="";
        my $transcript_id="";
        #my $exon_number="";
        my $gene_name="";
        my $gene_biotype=$a[5];
        for(my $i=0; $i<$Nr; $i=$i+2){
            $b[$i]=~s/\s//g; $b[$i]=~s/\;//g;
            if ($b[$i] eq "gene_id") {$gene_id=$b[$i+1];}
            elsif ($b[$i] eq "transcript_id") {$transcript_id=$b[$i+1];}
            elsif ($b[$i] eq "gene_name") {$gene_name=$b[$i+1];}
            elsif ($b[$i] eq "gene_biotype") {$gene_biotype=$b[$i+1];}
        }
        if ($gene_name eq "") {$gene_name=$gene_id;}
        my $id=$gene_id."\t".$a[0]."\t".$a[6];
        if (!exists $uniq{$id}) {
            $biotype{$id}=$gene_biotype;
            $Gname{$id}=$gene_name;
        }
        $uniq{$id}{$a[3]}{$a[4]}=1;
    }
}

foreach my $id (keys %uniq) {
    my %tmp;
    my %End;
    foreach my $left (keys %{$uniq{$id}}) {
        foreach my $right (keys %{$uniq{$id}{$left}}) {
            $tmp{$right}=1;
            $End{$right}=1;
        }
        $tmp{$left}=2;
    }
    
    foreach my $left (keys %{$uniq{$id}}) {
        foreach my $right (keys %{$uniq{$id}{$left}}) {
            foreach my $test_end (keys %End) {
                if (($left < $test_end) and ($test_end < $right)) {$End{$test_end}++;}
            }
        }
    }
    
    my @border=();
    my $Nr=0;
    my $last=-1;
    my %used;
    foreach my $pos (sort{$a <=> $b} keys %tmp) {
        if ($debug) { print join("\t",$Nr,$last,$pos),"\n"; }
        if ($last eq -1) { $last=$pos; }
        else {
            if (exists $End{$last}) {
                if ($End{$last} > 1) {
                    if ($tmp{$pos} eq 2) {
                        my $p5=($last+1);
                        my $p3=($pos-1);
                        if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; }
                    }
                    else {
                        my $p5=($last+1);
                        my $p3=($pos);
                        if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; }
                    }
                }
            }
            else {
                if ($tmp{$pos} eq 2) {
                    my $p5=($last);
                    my $p3=($pos-1);
                    if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; }
                } 
                else {
                    my $p5=($last);
                    my $p3=($pos);
                    if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; }
                }
            }
            $last=$pos;
        }
    }
    
    my @a=split("\t",$id);
    my $start_g=-1;
    my $end_g=-1;
    my $biotype_g="";
    my $gene_id="";
    my $chr="";
    my $strand="";
    if ($a[2] eq "+") {
        for (my $i=0; $i<$Nr; $i++) {
            print OUT join("\t",$a[1],"split","exon",$border[$i],$biotype{$id},$a[2],$Gname{$id},$a[0]."___".($i+1)."___".$Nr),"\n";
            my @b=split("\t",$border[$i]);
            if ($start_g eq -1) {$start_g=$b[0]; $chr=$a[1]; $biotype_g=$biotype{$id}; $gene_id=$a[0]; $strand=$a[2];}
            if ($end_g < $b[1]) {$end_g=$b[1];}
        }
    }
    else {
        for (my $i=0; $i<$Nr; $i++) {
            print OUT join("\t",$a[1],"split","exon",$border[$i],$biotype{$id},$a[2],$Gname{$id},$a[0]."___".($Nr - $i)."___".$Nr),"\n";
            my @b=split("\t",$border[$i]);
            if ($start_g eq -1) {$start_g=$b[0]; $chr=$a[1]; $biotype_g=$biotype{$id}; $gene_id=$a[0]; $strand=$a[2];}
            if ($end_g < $b[1]) {$end_g=$b[1];}
        }
    }
    print OUT2 join("\t",$chr,".","gene",$start_g,$end_g,$biotype_g,$strand,$Gname{$id},$gene_id),"\n";
}