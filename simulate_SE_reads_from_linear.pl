#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_mRNA_fa\"   \"output\"   \"\(optional\)name default=circRNA\"   \"\(optional\)read_length default=100\"   \"\(optional\)min_mRNA_length default=300\"   \"\(optional\)read_per_transcript default\=20\"  " if (@ARGV < 2);
my $filein1=$ARGV[0];	# mRNA.fa
my $fileout=$ARGV[1];	# simulated_read.fa, assuming Illumina Stranded RNA-Seq protocol.
my $newid="circRNA";
if (scalar(@ARGV) > 2) { $newid=$ARGV[2]; }
my $readlen=100;
if (scalar(@ARGV) > 3) { $readlen=$ARGV[3]; }
my $min_len=300;
if (scalar(@ARGV) > 4) { $min_len=$ARGV[4]; }
my $Nread=20;
if (scalar(@ARGV) > 5) { $Nread=$ARGV[5]; }
open IN, $filein1;
open OUT,">".$fileout;
while(<IN>) {
	chomp;
	s/^>//;
	my $header=$_;	# >ENSMUSG00000000544___ENSMUST00000027847___processed_transcript
	my @id=split("___",$header);
	if (scalar(@id) > 1) { $header=$id[1];}
	my $seq=<IN>;
	chomp $seq;
	my $len=length($seq);
	if ($len > $min_len) {
		for (my $i=0; $i<$Nread; $i++) {
			my $pos=int(rand($len - $readlen));
			my $sreadorg=substr($seq,$pos,$readlen);
			$sreadorg=~tr/[ATCGatcg]/[TAGCTAGC]/;
			my $sread=scalar reverse $sreadorg;
			print OUT ">Truseq_".$newid."_".$header."_".$len."_".$pos."_".$i."\/1\n";
			print OUT $sread,"\n";
		}
	}
}

