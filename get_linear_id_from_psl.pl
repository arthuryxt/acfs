#!/usr/bin/perl -w
use strict;
# report the ID of the sequences which have significant parts that can be aligned to multiple places
die "Usage: $0  \"psl_input\"   \"output\"  \"\(optional\)max_diff_score_ratio==0.02\"   \"\(optional\)overlap==0.2\"   \"\(optional\)largest_percentage=0.66\"  " if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my $max_diff_score_ratio=0.02;
if (scalar(@ARGV) > 2) { $max_diff_score_ratio = $ARGV[2]; }
my $overlap=0.2;
if (scalar(@ARGV) > 3) { $overlap = $ARGV[3]; }
if (($overlap < 0) or ($overlap >1)) { die "overlap MUST be within [0,1]. It is the max percentage two valid alignment can overlap.\n"; }
my $diff=0.66;
if (scalar(@ARGV) > 4) { $diff = $ARGV[4]; }
if (($diff < 0) or ($diff >1)) { die "largest_percentage MUST be within [0,1]. It is the max percentage the largest alignment allow to cover the whole sequence, otherwise considered as linear.\n"; }

sub check_overlap($$$$){
	if (($_[1] <= $_[2]) or ($_[3] <= $_[0])) { return(-1); }		# 01 23, 23 01
	if ($_[0] < $_[2]) {	# 0 2 1 3
		if ($_[3] < $_[1]) { return(abs($_[3] - $_[2] + 1));}	# 0 2 3 1
		else{ return(abs($_[1] - $_[2]) + 1);}					# 0 2 1 3
	}
	else {
		if ($_[1] < $_[3]) { return(abs($_[1] - $_[0] + 1));}	# 2 0 1 3 
		else{ return(abs($_[3] - $_[0]) + 1);}					# 2 0 3 1
	}
}

my $fileins=$filein.".s";
my $command="sort -k10,10 -k1,1nr $filein > $fileins";
system($command);

my %uniq;
my %bad;
open(IN, $fileins) or die "Cannot open psl file : $fileins\n";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
	if ($a[5] > 10) { next; }
    if (exists $bad{$a[0]}) { next;}
    if (exists $uniq{$a[9]}){
        my $f=0;
        foreach my $id (keys %{$uniq{$a[9]}}) {
            if ($f eq -1) { last; }
            my @b=split("\t",$uniq{$a[9]}{$id});
            my $myOverlap=check_overlap($a[11],$a[12],$b[11],$b[12]);
            my $critLen=$overlap * $a[10];
            if ($myOverlap < $critLen) {
                # seems to be a new partial alignment
            }
            else {
                # the almost-same part has another hit. If the score diff is big enough, store the one with higher score; otherwise, report to %bad
                my $diff_score=$a[0]-$b[0];
                if (abs($diff_score)/$a[0] < $max_diff_score_ratio) {
                    $f=-1;
                    $bad{$a[9]}=2;
                }
				else { $f=1; }
            }
        }
        if ($f eq 0) {
            $uniq{$a[9]}{$a[11]."\t".$a[12]}=join("\t",@a);
        }
    }
    else {
        if (($a[0] > ($diff * $a[10])) or (($a[12]-$a[11]) > ($diff * $a[10]) )) {
            $bad{$a[9]}=1;
        }
        else{
            $uniq{$a[9]}{$a[11]."\t".$a[12]}=join("\t",@a);
        }
    }
}
close IN;
$command="rm -rf $fileins";
system($command);
foreach my $id (keys %uniq) {
	my $cnt=0;
	foreach my $pos (keys %{$uniq{$id}}) { $cnt++; }
	if ($cnt ne 2) { $bad{$id}=3; }
}

open(OUT, ">".$fileout);
foreach my $id (sort keys %bad) {
    print OUT $id,"\n";
}


