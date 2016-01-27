#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input1\"   \"input_refFlat\"   \"\(optional\)window\"  \"\(optional\)debug\"" if (@ARGV < 2);
# check the if input_refFlat track is approximately the same (within the given window) by tracks in input1
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $window=10;
if (scalar(@ARGV) > 2){ $window=$ARGV[2]; }
my $debug=0;
if (scalar(@ARGV) > 3){ $debug=$ARGV[3]; }
my $command="rm -f Step3_MuSeg_finished";
system($command);

my %uniq;
open(IN1, $filein1);
while(<IN1>) {
    chomp;
    if (m/^track/) {}
    elsif (m/^#/) {}
    elsif (m/^@/) {}
    else {
        my @a=split("\t",$_);
        $a[1]=~s/chr//;
        $uniq{$a[1]}{$a[3]}{$a[2]}=$a[0];
    }
}
close IN1;

open(IN2, $filein2);
open(OUT1, ">".$filein2.".matched");
open(OUT2, ">".$filein2.".novel");
open(OUT22, ">".$filein2.".novel2");
my %selected;
while(<IN2>) {
    chomp;
    if (m/^track/) {}
    elsif (m/^#/) {}
    elsif (m/^@/) {}
    else {
        my @a=split("\t",$_);
        $a[2]=~s/chr//;
        my $f="";
        if (exists $uniq{$a[2]}) {
            foreach my $left (sort {$a <=> $b} keys %{$uniq{$a[2]}}) {
                if (abs($left - $a[4]) <= $window) {
                    foreach my $right(sort{$a <=> $b} keys %{$uniq{$a[2]}{$left}}) {
                        if (abs($right - $a[5]) <= $window) {
                            $f=$left."\t".$right;
                            last;
                        }
                    }
                    if ($f ne "") { last; }
                }
            }
        }
        if ($f ne "") {
            my @b=split("\t",$f);
            print OUT1 $uniq{$a[2]}{$b[0]}{$b[1]},"\t",join("\t",@a),"\n";
        }
        else {
            print OUT2 join("\t",@a),"\n";
            $selected{$a[1]}=1;
            # output in the same format as input1
            my @tmp=split(/\_/,$a[1]);
            print OUT22 join("\t",$a[1],$a[2],$tmp[1],$tmp[2],$tmp[3],$a[11],$a[12],$a[13],$a[3],$a[14],$a[15],$a[16],0,$a[15],"newid"),"\n";
        }
    }
}
close IN2;

open(IN3, $filein2.".fa");
open(OUT3, ">".$filein2.".novel.fa");
while (<IN3>) {
    chomp;
    if (m/^>/) {
        s/^>//;
        my $id=$_;
        my $seq=<IN3>;
        if (exists $selected{$id}) {
            print OUT3 ">".$id,"\n",$seq;
        }
    }
}

open(OUTFLAG,">Step3_MuSeg_finished");
print OUTFLAG "Step3_MuSeg_finished\n";
close OUTFLAG;

