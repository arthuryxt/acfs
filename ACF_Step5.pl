#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"BWA_sam\"   \"ref.fa\"    \"output\"    \"\(optional\) expand\"   \"\(optional\) Junc\"  \"\(optional\) NM==5%\"   \"\(optional\)stranded\"  \(debug == 1\)" if (@ARGV < 3);
# to estimate the expresison of circs after remap to pseudo-transcripts
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $fileout=$ARGV[2];
my $expand=150;
if(scalar(@ARGV) > 3) {$expand=$ARGV[3];}
my $JUN=6;
if(scalar(@ARGV) > 4) {$JUN=$ARGV[4];}
my $ER=0.05;
if(scalar(@ARGV) > 5) {$ER=$ARGV[5];}
my $stranded="no";
if(scalar(@ARGV) > 6) {
	$stranded=$ARGV[6];
	if (($stranded ne "no") and ($stranded ne "+") and ($stranded ne "-")) {
		die "stranded parameter accept only + or -, as reads should be mapped to sense or antisense of mRNA transcripts respectively. Leave this option when the sequencing is un-stranded."
	}
}
my $debug=0;
if(scalar(@ARGV) > 7) {$debug=$ARGV[7];}

my %Ref;
open IN2,$filein2;
while(<IN2>) {
    chomp;
    s/^>//;
    my $id=$_;
    my $seq=<IN2>;
    chomp $seq;
    $Ref{$id}=length($seq);
}
close IN2;

my %uniq;   # store processed ID
open IN, $filein1;
open OUT1,">".$fileout.".1";
open OUT2,">".$fileout.".2";
while(<IN>) {
    chomp;
    if ((m/^#/) or (m/^@/)) {next;}
    my @a=split("\t",$_);
	if (($a[2] eq "*") or ($a[1] eq 4)) { next; }
    my $Len=length($a[9]);
    my $strand="+";
    if ($a[1] eq 16) {$strand="-";}
    my @CIGAR_op=($a[5]=~m/[MSID]/g); 
    my @CIGAR_va=($a[5]=~m/\d+/g);
    my $start=0;
    my $length=0;
    my $NM=0;		# aligned edit-distance
    my $Nr=scalar(@CIGAR_op);
    my $readlen=length($a[9]);
	my $mismatch=0;	# full-read edit-distance
	my $exp=$expand;
	if ($Ref{$a[2]} < 2*$expand) { $exp=int($Ref{$a[2]}/2); }
	for(my $i=11; $i<scalar(@a); $i++) {
		if ($debug > 0 ) {print $a[$i],"\n";}
		if ($a[$i]=~/NM/) {
			my @tmp=split(/\:/,$a[$i]);
			if ($debug > 0 ) {print $a[$i],"\n";}
			$mismatch=$tmp[-1];
			$NM=$tmp[-1];
			last;
		}
	}
    if ($strand eq "+") {
        for(my $i=0; $i<$Nr; $i++){
            if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; }
            elsif ($CIGAR_op[$i] eq "I") {$NM+=$CIGAR_va[$i]; }
			elsif ($CIGAR_op[$i] eq "D") {$NM+=$CIGAR_va[$i];}
            elsif ($CIGAR_op[$i] eq "S") {
                if (($start > 0) or ($length > 0)) {last;}
                $start=$CIGAR_va[$i];
            }
        }
		for(my $i=0; $i<$Nr; $i++){
		    if ($CIGAR_op[$i] eq "S") {
				#$NM+=$CIGAR_va[$i];
				$mismatch+=$CIGAR_va[$i];
		    }
		}
    }
    else {
        for(my $i=$Nr-1; $i>=0; $i--){
            if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i];}
            elsif ($CIGAR_op[$i] eq "I") {$NM+=$CIGAR_va[$i];}
			elsif ($CIGAR_op[$i] eq "D") {$NM+=$CIGAR_va[$i];}
            elsif ($CIGAR_op[$i] eq "S") {
                if (($start > 0) or ($length > 0)) {last;}
                $start=$CIGAR_va[$i];		
            }
        }
		for(my $i=$Nr-1; $i>=0; $i--){
		    if ($CIGAR_op[$i] eq "S") {
				#$NM+=$CIGAR_va[$i];
				$mismatch+=$CIGAR_va[$i];
		    }
		}
    }
	if ($length > 0.8*$readlen) {
		my $SSR=$Ref{$a[2]}-$exp-$JUN;
		my $SSL=$Ref{$a[2]}-$exp-$readlen+$JUN;
		print OUT1 join("\t",$a[0],$strand,$a[2],$a[3],$a[5],$length,$NM,$mismatch,$readlen,($Ref{$a[2]}-$exp),$SSL,$SSR,$a[9]),"\n";
		if (($stranded eq "no") or ($stranded eq $strand)) {
			if (($SSL <= $a[3]) and ($a[3] <= $SSR) and ($NM <= 5) and (($mismatch/$readlen) < $ER)) {
				my $dist1=$a[3] - $SSL;
				my $dist2=$SSR - $a[3];
				my $dist=$dist1 < $dist2 ? $dist1 : $dist2;
				my $JNM=0;
				# should take CIGAR+MD signal here,
				$JNM=$mismatch;
				if ($dist < (1+$JUN)) {
					if ($JNM eq 0) {
						print OUT2 join("\t",$a[0],$strand,$a[2],$a[3],$a[5],$length,$NM,$mismatch,$readlen,($Ref{$a[2]}-$exp),$SSL,$SSR),"\n";
					}
				}
				else {
					print OUT2 join("\t",$a[0],$strand,$a[2],$a[3],$a[5],$length,$NM,$mismatch,$readlen,($Ref{$a[2]}-$exp),$SSL,$SSR),"\n";
				}
			}
		}
	}	
}

