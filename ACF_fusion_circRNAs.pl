#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"output\"   \"unmap.trans.splicing\"   \"\(optional\)searching_window==1000000\"    \"\(optional\)expr==unmap.trans.splicing.expr\"   \"\(optional\)debug==0\" " if (@ARGV < 2);
my $fileout=$ARGV[0];	# fusion_circRNAs
my $filein=$ARGV[1];	# unmap.trans.splicing
my $window=1000000;
if (scalar(@ARGV) > 2) { $window=$ARGV[2]; }
my $exprfile="";
if (scalar(@ARGV) > 3) { $exprfile=$ARGV[3]; }
my $debug=0;
if (scalar(@ARGV) > 4) { $debug=$ARGV[4]; }

my $fileinname=$filein.".tsloci.anno";
open(IN, $fileinname) or die "Cannot open inpur file: $fileinname. \n";
my %uniq;
my %TS;
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    # 11_119185973_+_X_107396317_-	yes	MCAM	ATG4A	119185849	119185973	ENSG00000076706___13___55	MCAM	-	protein_coding	107396241	107396317	ENSG00000101844___25___29	ATG4A	+	protein_coding
    # X_107372011_-_11_119181059_+	yes	ATG4A	MCAM	107372011	107372082	ENSG00000101844___10___29	ATG4A	+	protein_coding	119181059	119181176	ENSG00000076706___46___55	MCAM	-	protein_coding
    $uniq{$a[0]}=join("\t",@a);
    my @b=split(/\_/,$a[0]);
    $TS{$b[0]}{$b[2]}{$b[1]}=join("\t",@a);
    $TS{$b[3]}{$b[5]}{$b[4]}=join("\t",@a);
}
close IN;

my %EXPR;
my @Sample;
my $print_exp=0;
if ($exprfile ne "") {
    open(IN, $exprfile) or die "Cannot open expression file : $exprfile\n";
    my $header=<IN>;
    chomp $header;
    my @H=split("\t",$header);
    for(my $i=2; $i<scalar(@H); $i++) { $Sample[$i-2]=$H[$i]; }
    while (<IN>) {
        chomp;
        my @a=split("\t",$_);
        $EXPR{$a[0]}=join("\t",@a);
        $print_exp++;
    }
    close IN;
}


open(OUT, ">".$fileout);
my %used;
my $cnt=0;
print OUT join("\t","Fusion_circRNA_id","Fusion_anno","Junction_1","Junction_2",join("\t",@Sample)),"\n";
foreach my $id (keys %uniq) {
    my @a=split("\t",$uniq{$id});
    my @b=split(/\_/,$a[0]);
    foreach my $pos (keys %{$TS{$b[3]}{$b[5]}}) {
        my @a2=split("\t",$TS{$b[3]}{$b[5]}{$pos});
        my @b2=split(/\_/,$a2[0]);
        if ((!exists $used{$a2[0]}) and ($a2[0] ne $id) and ($b[0] eq $b2[3]) and ($b[2] eq $b2[5]) and (abs($pos - $b[4]) < $window) and (abs($b2[4] - $b[1]) < $window)) {
            #if ($debug > 0) { print join("\t",@a),"\n",join("\t",@a2),"\n\n"; }
            if (($b[2] eq "+") and ($b[5] eq "+")){
                if (($b[1] > $b2[4]) and ($b[4] < $b2[1])) {
                    if ($print_exp > 0) {
                        if ((exists $EXPR{$a[0]}) and (exists $EXPR{$a2[0]})){
                            $cnt++;
                            my @exp1=split("\t",$EXPR{$a[0]});
                            my @exp2=split("\t",$EXPR{$a2[0]});
                            my $expr=$exp1[2] > $exp2[2] ? $exp2[2] : $exp1[2];
                            for (my $i=3; $i<scalar(@exp1); $i++) {
                                $expr=$expr."\t".($exp1[$i] > $exp2[$i] ? $exp2[$i] : $exp1[$i]);
                            }
                            print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\t",$expr,"\n";
                            if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                        }
                    }
                    else{
                        print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\n";
                        if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                    }
                }
            }
            elsif (($b[2] eq "-") and ($b[5] eq "-") ){
                if(($b[1] < $b2[4]) and ($b[4] > $b2[1])) {
                    if ($print_exp > 0) {
                        if ((exists $EXPR{$a[0]}) and (exists $EXPR{$a2[0]})){
                            $cnt++;
                            my @exp1=split("\t",$EXPR{$a[0]});
                            my @exp2=split("\t",$EXPR{$a2[0]});
                            my $expr=$exp1[2] > $exp2[2] ? $exp2[2] : $exp1[2];
                            for (my $i=3; $i<scalar(@exp1); $i++) {
                                $expr=$expr."\t".($exp1[$i] > $exp2[$i] ? $exp2[$i] : $exp1[$i]);
                            }
                            print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\t",$expr,"\n";
                            if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                        }
                    }
                    else{
                        print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\n";
                        if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                    }
                }
            }
            elsif (($b[2] eq "-") and ($b[5] eq "+") ){
                if (($b[1] < $b2[4]) and ($b[4] < $b2[1])) {
                    if ($print_exp > 0) {
                        if ((exists $EXPR{$a[0]}) and (exists $EXPR{$a2[0]})){
                            $cnt++;
                            my @exp1=split("\t",$EXPR{$a[0]});
                            my @exp2=split("\t",$EXPR{$a2[0]});
                            my $expr=$exp1[2] > $exp2[2] ? $exp2[2] : $exp1[2];
                            for (my $i=3; $i<scalar(@exp1); $i++) {
                                $expr=$expr."\t".($exp1[$i] > $exp2[$i] ? $exp2[$i] : $exp1[$i]);
                            }
                            print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\t",$expr,"\n";
                            if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                        }
                    }
                    else{
                        print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\n";
                        if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                    }
                }
            }
            elsif (($b[2] eq "+") and ($b[5] eq "-") ){
                if(($b[1] > $b2[4]) and ($b[4] > $b2[1])) {
                    if ($print_exp > 0) {
                        if ((exists $EXPR{$a[0]}) and (exists $EXPR{$a2[0]})){
                            $cnt++;
                            my @exp1=split("\t",$EXPR{$a[0]});
                            my @exp2=split("\t",$EXPR{$a2[0]});
                            my $expr=$exp1[2] > $exp2[2] ? $exp2[2] : $exp1[2];
                            for (my $i=3; $i<scalar(@exp1); $i++) {
                                $expr=$expr."\t".($exp1[$i] > $exp2[$i] ? $exp2[$i] : $exp1[$i]);
                            }
                            print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\t",$expr,"\n";
                            if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                        }
                    }
                    else{
                        print OUT "fusion-circ-".$cnt,"\t",$a[2]."|".$a[3],"\t",$a[0],"\t",$a2[0],"\n";
                        if ($debug > 0) {print OUT join("\t",@a),"\n",join("\t",@a2),"\n";}
                    }
                }
            }
        }
    }
    $used{$id}=1;
}

