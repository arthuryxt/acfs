#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"output\"    \"fusion.gtf\"   \"\(optional\)mimic_stranded_Seq_reads==1(default) othersize==0\" " if (@ARGV ne 2);
my $fileout=$ARGV[0];	# simu_fusion_circRNA_id
my $filein1=$ARGV[1];	# simu_fusion.gtf
my $mimic=1;
if (scalar(@ARGV) > 2) { $mimic=$ARGV[2]; }
if (($mimic ne 0) and ($mimic ne 1)) {
    die "mimic_stranded_Seq_reads can only be 1 or 0\n";
}

open(IN, $filein1) or die "Cannot open fusion gtf file : $filein1";
my %uniq;
my %anno;
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    # 7       fusion  exon    76143264        76143393        .       +       .       gene_id "NM_030570__2__2__NM_001297713__3__5__NM_030570";
    my @b=split(/\"/,$a[8]);
    my $Nr=scalar(@b);
    my $gene_id="";
    my $transcript_id="";
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
    if ($transcript_id eq "") {$transcript_id=$gene_id;}
    my @c=split(/\_\_/,$gene_id);
    $anno{join("__",$c[0],$c[1],$c[2],$c[3],$c[4],$c[5])}=1;
    if (exists $uniq{$c[-1]}) {
        my @c1=split("\t",$uniq{$c[-1]});
        if ($a[3] < $c1[1]) { $c1[1]=$a[3]; }
        if ($a[4] > $c1[2]) { $c1[2]=$a[4]; }
        $uniq{$c[-1]}=join("\t",@c1);
    }
    else { $uniq{$c[-1]}=join("\t",$a[0],$a[3],$a[4],$a[6]); }    
}
close IN;

open(OUT, ">".$fileout);
foreach my $id (sort keys %anno) {
    my @a=split(/\_\_/,$id);
    my @L=split("\t",$uniq{$a[0]});
    my @R=split("\t",$uniq{$a[3]});
    if ($mimic eq 1) {
        if ($L[3] eq "+") { $L[3]="-"; }
        else { $L[3]="+"; }
        if ($R[3] eq "+") { $R[3]="-"; }
        else { $R[3]="+"; }
    }
    if (($L[3] eq "-") and ($R[3] eq "-")) {
        print OUT join("_",$L[0],$L[1],$L[3],$R[0],$R[2],$R[3]),"\t",$id,"\n";
        print OUT join("_",$R[0],$R[1],$R[3],$L[0],$L[2],$L[3]),"\t",$id,"\n";
    }
    elsif (($L[3] eq "+") and ($R[3] eq "+")) {
        print OUT join("_",$L[0],$L[2],$L[3],$R[0],$R[1],$R[3]),"\t",$id,"\n";
        print OUT join("_",$R[0],$R[2],$R[3],$L[0],$L[1],$L[3]),"\t",$id,"\n";
    }
    elsif (($L[3] eq "+") and ($R[3] eq "-")) {
        print OUT join("_",$L[0],$L[2],$L[3],$R[0],$R[2],$R[3]),"\t",$id,"\n";
        print OUT join("_",$R[0],$R[1],$R[3],$L[0],$L[1],$L[3]),"\t",$id,"\n";
    }
    elsif (($L[3] eq "-") and ($R[3] eq "+")) {
        print OUT join("_",$L[0],$L[1],$L[3],$R[0],$R[1],$R[3]),"\t",$id,"\n";
        print OUT join("_",$R[0],$R[2],$R[3],$L[0],$L[2],$L[3]),"\t",$id,"\n";
    }
}

