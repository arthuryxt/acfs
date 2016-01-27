#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_circRNA_fa\"   \"output\"   \"\(optional\)name default=circRNA\"   \"\(optional\)read_length default=100\"   \"\(optional\)min_circRNA_length default=50\"   \"\(optional\)read_per_transcript default\=20\"  \"\(optional\)overlap_junction default\=6\" " if (@ARGV < 2);
my $filein1=$ARGV[0];	# circ.fa
my $fileout=$ARGV[1];	# simulated_read.fa, assuming Illumina Stranded RNA-Seq protocol.
my $newid="circRNA";
if (scalar(@ARGV) > 2) { $newid=$ARGV[2]; }
my $readlen=100;
if (scalar(@ARGV) > 3) { $readlen=$ARGV[3]; }
my $min_len=100;
if (scalar(@ARGV) > 4) { $min_len=$ARGV[4]; }
my $Nread=20;
if (scalar(@ARGV) > 5) { $Nread=$ARGV[5]; }
my $overlap=6;
if (scalar(@ARGV) > 6) { $overlap=$ARGV[6]; }
my $debug=0;
if (scalar(@ARGV) > 7) { $debug=$ARGV[7]; }
open IN, $filein1;
open OUT,">".$fileout;
while(<IN>) {
	chomp;
	s/^>//;
	my $id=$_;
	my $seq=<IN>;
	chomp $seq;
	my $len=length($seq);
	my $pseq=$seq.$seq;
	if ($len > $min_len) {
		for (my $i=0; $i<$Nread; $i++) {
			my $pos=$len - $readlen + $overlap + int(rand($readlen-$overlap));
			my $covered=0;
			if (($pos < ($len - $overlap)) and (($len + $overlap) < ($pos + $readlen))) { $covered = 1;};
			while (($pos < 0) or (($pos + $readlen) >= (2*$len - $overlap)) or ($covered != 1)) {
				$pos=$len - $readlen + $overlap + int(rand($readlen-$overlap));
				#if ((($len-$readlen+$overlap) < $pos) and ($pos < ($len-$overlap))){$covered=1;}
				if (($pos < ($len - $overlap)) and (($len + $overlap) < ($pos + $readlen))) { $covered = 1;}
				if ($debug) {print join("\t",$id,$len,$readlen,$pos),"\n";}
			}
			my $sreadorg=substr($pseq,$pos,$readlen);
			$sreadorg=~tr/[ATCGatcg]/[TAGCTAGC]/;
			my $sread=scalar reverse $sreadorg;
			print OUT ">Truseq_".$newid."_".$id."_".$len."_".$pos."\/1\n";
			print OUT $sread,"\n";
		}
	}
}

