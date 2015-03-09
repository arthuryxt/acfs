#!/usr/bin/perl -w
use strict;
# convert bed_formatted pre_defined_circRNAs into proprietray sum format
die "Usage: $0   \"bed_input\"  \(optional\)output" if (@ARGV < 1);
my $filein=$ARGV[0];
my $fileout="pre_defined_circRNA.sum";
if (scalar(@ARGV) > 1) {$fileout=$ARGV[1];}
open IN,$filein;
open OUT,">".$fileout;

my $count=0;
while(<IN>) {
	chomp;
	if (m/^>/){ next; }
	if (m/^@/){ next; }
	if (m/^#/){ next; }
	if (m/^track/){ next; }
	my @a=split("\t",$_);
	$count++;
	my $id="tpid-".$count."/1__1";
	my $chr=$a[0];
	$chr=~s/chromosome//;
	$chr=~s/chr//;
	my $SS=100;
	my $HS=50;
	my $strand="+";
	if ($a[5] eq "+") {
		$strand="-";
		if (scalar(@a) >= 12) {
			my @b=split(/\|/,$a[3]);
			if (scalar(@b) eq 4) {
				$SS=$b[1];
				$HS=sprintf("%.2f",$SS/2);
			}
			print OUT join("\t",$id,60,$chr,30,30,($a[2]-30+1),$a[2],$strand,0,30,($a[1]+1),($a[1]+30),$strand,($a[1]-$a[2]+1),$SS,$HS,$HS,$a[5],0,0,0,10,0),"\n";
		}
		elsif (scalar(@a) >= 6 ) {
			print OUT join("\t",$id,60,$chr,30,30,($a[2]-30+1),$a[2],$strand,0,30,($a[1]+1),($a[1]+30),$strand,($a[1]-$a[2]+1),$SS,$HS,$HS,$a[5],0,0,0,10,0),"\n";
		}
	}
	else {
		# $a[5] eq "-"
		if (scalar(@a) >= 12) {
			my @b=split(/\|/,$a[3]);
			if (scalar(@b) eq 4) {
				$SS=$b[1];
				$HS=sprintf("%.2f",$SS/2);
			}
			print OUT join("\t",$id,60,$chr,0,30,($a[2]-30+1),($a[2]),$strand,30,30,($a[1]+1),($a[1]+30),$strand,($a[1]-$a[2]+1),$SS,$HS,$HS,$a[5],0,0,0,10,0),"\n";
		}
		elsif (scalar(@a) >= 6 ) {
			print OUT join("\t",$id,60,$chr,0,30,($a[2]-30+1),($a[2]),$strand,30,30,($a[1]+1),($a[1]+30),$strand,($a[1]-$a[2]+1),$SS,$HS,$HS,$a[5],0,0,0,10,0),"\n";
		}
	}
	
}
