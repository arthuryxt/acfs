#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_gtf\"   \"output\"   \"\(optional\)Nr_circRNA_per_transcript\"    " if (@ARGV < 2);
my $gtf=$ARGV[0];
my $fileout=$ARGV[1];
my $GenerateNr=2;
if (scalar(@ARGV) > 2) { $GenerateNr=$ARGV[2]; }
my %uniq;
my %biotype;
my %Gname;
my %ExonCnt;
open IN,$gtf;
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
		if ($transcript_id eq "") {$transcript_id=$gene_id;}
        my $id=$gene_id."\t".$transcript_id."\t".$a[0]."\t".$a[6];
        if (!exists $uniq{$id}) {
            $biotype{$id}=$gene_biotype;
            $Gname{$id}=$gene_name;
        }
		$uniq{$id}{$a[3]}=$a[4];
		$ExonCnt{$id}++;
    }
}
my %circ;
foreach my $id (keys %uniq){
	if ($ExonCnt{$id} > 2) {
		my @Start;
		my @End;
		my $cnt=0;
		foreach my $pos (sort{$a <=> $b} keys %{$uniq{$id}}) {
			$Start[$cnt]=$pos;
			$End[$cnt]=$uniq{$id}{$pos};
			$cnt++;
		}
		if ($cnt eq 3) {
			# self-circulate the 2nd exon
			my @a=split("\t",$id);
			my $circ_id=join("_",$a[2],$End[1],$Start[1],($Start[1]-$End[1]));
			if (exists $circ{$circ_id}) {}
			else {
				my $tmp_info=join(";","gene_id \"".$circ_id."\"", " transcript_id \"".$circ_id."\""," exon_number \"1\""," gene_name \"".$Gname{$id}."\"");
				$circ{$circ_id}=join("\t",$a[2],$biotype{$id},"exon",$Start[1],$End[1],".",$a[3],".",$tmp_info);
			}
		}
		elsif ($cnt > 3) {
			# generate $GenerateNr circRNAs
			my @a=split("\t",$id);
			for(my $trail=0; $trail<$GenerateNr; $trail++){
				my $pos1=int(rand($cnt-1));
				my $pos2=int(rand($cnt-1));
				my $min=$pos1 < $pos2 ? $pos1 : $pos2;
				my $max=$pos1 > $pos2 ? $pos1 : $pos2;
				if ($min eq 0) {$min++;}
				if ((0 < $min) and ($min < $max) and ($max < ($cnt-1))) {
					my $circ_id=join("_",$a[2],$End[$max],$Start[$min],($Start[$min] - $End[$max]));
					if (exists $circ{$circ_id}) { } #$trail--; }
					else{
						my $exoncnt=1;
						my $tmp_info=join(";","gene_id \"".$circ_id."\"", " transcript_id \"".$circ_id."\""," exon_number \"".$exoncnt."\""," gene_name \"".$Gname{$id}."\"");
						my $tmp_gtf=join("\t",$a[2],$biotype{$id},"exon",$Start[$min],$End[$min],".",$a[3],".",$tmp_info);
						for(my $i=$min+1; $i<=$max; $i++) {
							$exoncnt++;
							$tmp_info=join(";","gene_id \"".$circ_id."\"", " transcript_id \"".$circ_id."\""," exon_number \"".$exoncnt."\""," gene_name \"".$Gname{$id}."\"");
							$tmp_gtf=$tmp_gtf."\n".join("\t",$a[2],$biotype{$id},"exon",$Start[$i],$End[$i],".",$a[3],".",$tmp_info);
							
						}
						$circ{$circ_id}=$tmp_gtf;
					}
				}
				#else { $trail--; }
			}
		}
	}
}
open OUT,">".$fileout;
foreach my $id (sort keys %circ){
	print OUT $circ{$id},"\n";
}
