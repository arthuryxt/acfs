#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_gtf\"   \"output\"  \"\(optional\)Nr_fusion_transcripts==10\"   \"\(optional\)min_Exons==5\"  \"\(optional\)Nr_circRNA_per_transcript==2\"    " if (@ARGV < 2);
my $gtf=$ARGV[0];
my $fileout=$ARGV[1];
my $NrFT=10;
if (scalar(@ARGV) > 2) { $NrFT=$ARGV[2]; }
my $min_Exs=5;
if (scalar(@ARGV) > 3) { $min_Exs=$ARGV[3]; }
my $GenerateNr=2;
if (scalar(@ARGV) > 4) { $GenerateNr=$ARGV[4]; }
my %uniq;
my %biotype;
my %Gname;
my %ExonCnt;
open IN,$gtf;
my @Gid;
my $gcnt=0;
while(<IN>) {
	chomp;
	my @a=split("\t",$_);
    # MT      protein_coding  exon    3307    4262    .       +       .        gene_id "ENSG00000198888"; transcript_id "ENST00000361390"; exon_number "1"; gene_name "MT-ND1"; gene_biotype "protein_coding"; transcript_name "MT-ND1-201";
    if (($a[2] eq "exon") and ($a[0]!~m/NT/) and ($a[0]!~m/Y/)){
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
		if ($gene_id=~m/^NR/) { next; }
        if ($gene_name eq "") {$gene_name=$gene_id;}
		if ($transcript_id eq "") {$transcript_id=$gene_id;}
        my $id=$gene_id."\t".$transcript_id."\t".$a[0]."\t".$a[6];
        if (!exists $uniq{$id}) {
            $biotype{$id}=$gene_biotype;
            $Gname{$id}=$gene_name;
			$Gid[$gcnt]=$id;
			$gcnt++;
        }
		$uniq{$id}{$a[3]}=$a[4];
		$ExonCnt{$id}++;
    }
}
close IN;
open(OUT, ">".$fileout);
my %circ;
for (my $i=1; $i<=$NrFT; $i++) {
	my $left_gid=int(rand($gcnt-1));
	while ($ExonCnt{$Gid[$left_gid]} < $min_Exs) { $left_gid=int(rand($gcnt-1)); }
	my @aL=split("\t",$Gid[$left_gid]);		# Gid	Tid	chr	strand
	my $right_gid=int(rand($gcnt-1));
	my @aR=split("\t",$Gid[$right_gid]);	# Gid	Tid	chr	strand
	while (($left_gid eq $right_gid) or ($ExonCnt{$Gid[$right_gid]} < $min_Exs) or ($aL[2] eq $aR[2])) {
		$right_gid=int(rand($gcnt-1));
		@aR=split("\t",$Gid[$right_gid]);
	}
    my $left_Ex1=2+int(rand($ExonCnt{$Gid[$left_gid]}-3));
	my $left_Ex2=2+int(rand($ExonCnt{$Gid[$left_gid]}-3));
	while ($left_Ex1 > $left_Ex2) {
        $left_Ex1=2+int(rand($ExonCnt{$Gid[$left_gid]}-3));
		$left_Ex2=2+int(rand($ExonCnt{$Gid[$left_gid]}-3));
    }
    my $right_Ex1=2+int(rand($ExonCnt{$Gid[$right_gid]}-3));
	my $right_Ex2=2+int(rand($ExonCnt{$Gid[$right_gid]}-3));
	while ($right_Ex1 > $right_Ex2) {
        $right_Ex1=2+int(rand($ExonCnt{$Gid[$right_gid]}-3));
		$right_Ex2=2+int(rand($ExonCnt{$Gid[$right_gid]}-3));
    }
	
	my $circ_id0=join("_",$aL[1],$left_Ex1,$left_Ex2,$aR[1],$right_Ex1,$right_Ex2);
	my $circ_id1=join("_",$aR[1],$right_Ex1,$right_Ex2,$aL[1],$left_Ex1,$left_Ex2);
	if ((exists $circ{$circ_id0}) or (exists $circ{$circ_id1})) {
        $i--;
		next;
    }
    $circ{$circ_id0}=1;
	$circ{$circ_id1}=1;
	my @StartL;
	my @EndL;
	my $cntL=0;
	foreach my $pos (sort{$a <=> $b} keys %{$uniq{$Gid[$left_gid]}}) {
		$StartL[$cntL]=$pos;
		$EndL[$cntL]=$uniq{$Gid[$left_gid]}{$pos};
		$cntL++;
	}
	my @StartR;
	my @EndR;
	my $cntR=0;
	foreach my $pos (sort{$a <=> $b} keys %{$uniq{$Gid[$right_gid]}}) {
		$StartR[$cntR]=$pos;
		$EndR[$cntR]=$uniq{$Gid[$right_gid]}{$pos};
		$cntR++;
	}
	my $circ_id=join("__",$aL[1],$left_Ex1,$left_Ex2,$aR[1],$right_Ex1,$right_Ex2,$aL[1]);
	my $exoncnt=1;
	my $tmp_gtf="";
	for(my $j=$left_Ex1; $j<=$left_Ex2; $j++) {
		my $tmp_info=join(";","gene_id \"".$circ_id."\"", " transcript_id \"".$circ_id."\""," exon_number \"".$exoncnt."\""," gene_name \"".$circ_id."\"");
		if ($tmp_gtf eq "") {
            $tmp_gtf=join("\t",$aL[2],"fusion","exon",$StartL[$j],$EndL[$j],".",$aL[3],".",$tmp_info);
        }
        else {
			$tmp_gtf=$tmp_gtf."\n".join("\t",$aL[2],"fusion","exon",$StartL[$j],$EndL[$j],".",$aL[3],".",$tmp_info);
		}
		$exoncnt++;
	}
	print OUT $tmp_gtf,"\n";
	$circ_id=join("__",$aL[1],$left_Ex1,$left_Ex2,$aR[1],$right_Ex1,$right_Ex2,$aR[1]);
	$tmp_gtf="";
	for(my $j=$right_Ex1; $j<=$right_Ex2; $j++) {
		my $tmp_info=join(";","gene_id \"".$circ_id."\"", " transcript_id \"".$circ_id."\""," exon_number \"".$exoncnt."\""," gene_name \"".$circ_id."\"");
		if ($tmp_gtf eq "") {
            $tmp_gtf=join("\t",$aR[2],"fusion","exon",$StartR[$j],$EndR[$j],".",$aR[3],".",$tmp_info);
        }
        else {
			$tmp_gtf=$tmp_gtf."\n".join("\t",$aR[2],"fusion","exon",$StartR[$j],$EndR[$j],".",$aR[3],".",$tmp_info);
		}
		$exoncnt++;
	}
	print OUT $tmp_gtf,"\n";
}
