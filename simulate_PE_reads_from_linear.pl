#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_mRNA_fa\"   \"output\"   \"\(optional\)name default=circRNA\"   \"\(optional\)read_length default=100\"   \"\(optional\)insert_size default=200\"   \"\(optional\)insert_size_SD default=50\"   \"\(optional\)min_mRNA_length default=300\"   \"\(optional\)read_per_transcript default\=20\"  " if (@ARGV < 2);
my $filein1=$ARGV[0];	# mRNA.fa
my $fileout=$ARGV[1];	# simulated_read.fa, assumed Illumina Stranded RNA-Seq protocol.
my $newid="circRNA";
if (scalar(@ARGV) > 2) { $newid=$ARGV[2]; }
my $readlen=100;
if (scalar(@ARGV) > 3) { $readlen=$ARGV[3]; }
my $insertSize=200;
if (scalar(@ARGV) > 4) { $insertSize=$ARGV[4]; }
my $insertSD=50;
if (scalar(@ARGV) > 5) { $insertSD=$ARGV[5]; }
my $min_len=300;
if (scalar(@ARGV) > 6) { $min_len=$ARGV[6]; }
my $Nread=20;
if (scalar(@ARGV) > 7) { $Nread=$ARGV[7]; }
open IN, $filein1;
open OUT,">".$fileout;

sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}

while(<IN>) {
	#my $trailMax=100;
	#my $trail=0;
	chomp;
	my @id=split("___",$_);
	my $seq=<IN>;
	chomp $seq;
	my $len=length($seq);
	if ($len > $min_len) {
		for (my $i=0; $i<$Nread; $i++) {
			my @info=gaussian_rand;
			my $tmp_insertSize=int($info[0]*$insertSD + $insertSize);
			if ($tmp_insertSize < $readlen) {$tmp_insertSize = $readlen};
			my $pos=int(rand($len - $tmp_insertSize));
			while(($pos <= 0)or(($pos + $tmp_insertSize) > $len)) {
				@info=gaussian_rand;
				$tmp_insertSize=int($info[0]*$insertSD + $insertSize);
				if ($tmp_insertSize < $readlen) {$tmp_insertSize = $readlen};
				$pos=int(rand($len - $tmp_insertSize));
				#$trail++;
				#if ($trail >= $trailMax) {last;}
			}
			#if ($trail >= $trailMax) {last;}
			my $sread1org=substr($seq,$pos,$readlen);
			$sread1org=~tr/[ATCGatcg]/[TAGCTAGC]/;
			my $sread1=scalar reverse $sread1org;
			my $sread2=substr($seq,($pos+$tmp_insertSize-$readlen),$readlen);
			print OUT ">Truseq_".$newid."_".$id[1]."_".$len."_".$pos."_".$tmp_insertSize."\/1\n";
			print OUT $sread1,"\n";
			print OUT ">Truseq_".$newid."_".$id[1]."_".$len."_".$pos."_".$tmp_insertSize."\/2\n";
			print OUT $sread2,"\n";
		}
	}
}

