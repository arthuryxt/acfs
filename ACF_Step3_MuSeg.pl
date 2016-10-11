#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input1\"   \"input_refFlat\"   \"\(optional\)window\"  \"\(optional\)debug\"" if (@ARGV < 2);
# check the if input_refFlat track is approximately the same (within the given window) by tracks in input1
my $filein1=$ARGV[0];
my $filein2=$ARGV[1];
my $window=0;
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
        #$uniq{$a[1]}{$a[3]}{$a[2]}{$a[0]}=join("\t",@a);
        $uniq{$a[0]}=join("\t",@a);
        # 14_36777395_36769991_-7404	14	36777395	36769991	-7404	19.04	8.96	10.08	-	2	0	2	-3	0	newid-203725__1
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
        # 14_36777395_36769991_-7404_|	14_36777395_36769991_-7404	chr14	-	36769991	36777395	36769991	36769991	2	36769991,36777312	36770029,36777395	19.04	10.08	8.96	2	2	0	newid-33948__2,newid-267135__1,newid-334715__1,newid-645333__1,newid-766734__1,newid-836151__1,newid-1019877__1,newid-1243254__1,newid-1366244__1,newid-1628932__1,newid-1819963__1
        my @a=split("\t",$_);
        if (exists $uniq{$a[1]}) {
            my @b=split("\t",$uniq{$a[1]});
            $b[-1]=$b[-1].",".$a[-1];
            $uniq{$a[1]}=join("\t",@b);
        }
        else {
            print OUT2 join("\t",@a),"\n";
            $selected{$a[1]}=1;
            # output in the same format as input1
            my @tmp=split(/\_/,$a[1]);
            print OUT22 join("\t",$a[1],$a[2],$tmp[1],$tmp[2],($tmp[2]-$tmp[1]),$a[11],$a[12],$a[13],$a[3],$a[14],$a[15],$a[16],0,$a[15],$a[-1]),"\n";
        }
    }
}
close IN2;
foreach my $id (sort keys %uniq) {
    print OUT1 $uniq{$id},"\n";
}

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

