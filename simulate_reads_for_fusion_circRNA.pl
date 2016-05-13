#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"input_fa\"   \"output\"   \"\(optional\)read_length==100\"  \"\(optional\)border==20\"   \"\(optional\)sample_reads==10\"   \"\(optional\)strandness==-\"   \"\(optional\)PE==0\"" if (@ARGV < 2);
my $filein=$ARGV[0];    # Must be the sequences generated for fusion_circRNA: there are two sequences for every translocation
my $fileout=$ARGV[1];
my $readLen=100;
if (scalar(@ARGV) > 2) { $readLen=$ARGV[2]; }
my $border=20;
if (scalar(@ARGV) > 3) { $border=$ARGV[3]; }
my $Nr=10;
if (scalar(@ARGV) > 4) { $Nr=$ARGV[4]; }
my $strandness="-";
if (scalar(@ARGV) > 5) { $strandness=$ARGV[5]; }
my $PE=0;   # set to 0 to simulate only SE data; set to 1 to simulate PE data
if (scalar(@ARGV) > 6) { $PE=$ARGV[6]; }
open(IN, $filein) or die "Cannot open the input fasta file: $filein";
open(OUT, ">".$fileout);
open(OUT1, ">".$fileout.".pseudo");
my %uniq;
while (<IN>) {
    chomp;
    s/^>//;
    #>NM_001242672_1_3_NM_001321384_10_10_NM_001242672___NM_001242672_1_3_NM_001321384_10_10_NM_001242672___split
    my @a=split(/\_\_/,$_);
    my $seq=<IN>;
    chomp $seq;
    my $id=join("__",$a[0],$a[1],$a[2],$a[3],$a[4],$a[5]);
    if ($a[0] eq $a[6]) { $uniq{$id}{1}=$seq; }
    else { $uniq{$id}{2}=$seq; }
}
close IN;
my @Pos;
my $step=int(($readLen-$border-$border-1)/$Nr);
for(my $i=0; $i<$Nr; $i++){ $Pos[$i]=0 + $i * $step; }
foreach my $id (keys %uniq) {
    if ((exists $uniq{$id}{1}) and (exists $uniq{$id}{2})){    
        my $seq1=$uniq{$id}{1};
        my $seq2=$uniq{$id}{2};
        my $tmp1=substr($seq1,(0-$readLen+$border));
        my $tmp2=substr($seq2,0,($readLen-$border));
        my $tmp=$tmp1.$tmp2;
        print OUT1 ">".$id."_1_2","\n",$tmp,"\n";
        for(my $i=0; $i<$Nr; $i++){
            my $tmpseq=substr($tmp,$Pos[$i],$readLen);
            if ($strandness eq "-") {
                my $cseq=$tmpseq;
                $cseq=~tr/[atcgATCG]/[TAGCTAGC]/;
                my $rc=scalar reverse $cseq;
                print OUT ">Truseq_fcirc_".$id."_1_2_".$i."\/1","\n",$rc,"\n";
                if ($PE eq 1) {
                    print OUT ">Truseq_fcirc_".$id."_1_2_".$i."\/2","\n",$tmpseq,"\n";
                }
            }
            else {
                print OUT ">Truseq_fcirc_".$id."_1_2_".$i."\/1","\n",$tmpseq,"\n";
                if ($PE eq 1) {
                    $tmpseq=~tr/[atcgATCG]/[TAGCTAGC]/;
                    my $rc=scalar reverse $tmpseq;
                    print OUT ">Truseq_fcirc_".$id."_1_2_".$i."\/2","\n",$rc,"\n";
                }
            }       
        }
        $tmp2=uc substr($seq2,(0-$readLen+$border));
        $tmp1=uc substr($seq1,0,($readLen-$border));
        $tmp=$tmp2.$tmp1;
        print OUT1 ">".$id."__2__1","\n",$tmp,"\n";
        for(my $i=0; $i<$Nr; $i++){
            my $tmpseq=substr($tmp,$Pos[$i],$readLen);
            if ($strandness eq "-") {
                my $cseq=$tmpseq;
                $cseq=~tr/[atcgATCG]/[TAGCTAGC]/;
                my $rc=scalar reverse $cseq;
                print OUT ">Truseq_fcirc_".$id."_2_1_".$i."\/1","\n",$rc,"\n";
                if ($PE eq 1) {
                    print OUT ">Truseq_fcirc_".$id."_2_1_".$i."\/2","\n",$tmpseq,"\n";
                }
            }
            else {
                print OUT ">Truseq_fcirc_".$id."_2_1_".$i."\/1","\n",$tmpseq,"\n";
                if ($PE eq 1) {
                    $tmpseq=~tr/[atcgATCG]/[TAGCTAGC]/;
                    my $rc=scalar reverse $tmpseq;
                    print OUT ">Truseq_fcirc_".$id."_2_1_".$i."\/2","\n",$rc,"\n";
                }
            }
        }
    }
}
close OUT;
