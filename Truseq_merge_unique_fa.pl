#!/usr/bin/perl -w
use strict;
# unique "Truseq" style fasta/fq reads, output fasta
die "Usage: $0    \"output_name\"    \"header_name\"    \"fasta-1\"    \"fasta-2 ...\"" if (@ARGV < 3);
my $fileout=$ARGV[0];
my $header=$ARGV[1];    # batch1 for example
my $fileNR=scalar(@ARGV);
my %uniq;
my %SEQ;
my %Sample;
open OUT,">".$fileout;
open OUT2,">".$fileout."_expr";
for(my $i=2; $i<$fileNR; $i++) {
    open IN,$ARGV[$i];
    print "reading file : ",$ARGV[$i],"\n";
    while(<IN>) {
        chomp;
        if (m/^>/) {
            my $id=$_;
            #$id=~s/^>//;
            my $seq=<IN>;
            chomp $seq;
            my @a=split(/\_/,$id);
            # >Truseq_F4.RiboZero_ACAGTGA_DJG6PNM1_230_1_1101_18024_2222/1__3828
            my $count=1;
            if ($id=~m/\_\_/) {
                my @b=split(/\_\_/,$id);
                if (scalar(@b) > 1) {$count=$b[-1];}
            }
            if (exists $uniq{$seq}) {
                $uniq{$seq}{$a[1]}+=$count;
                $SEQ{$seq}+=$count;
            }
            else {
                $uniq{$seq}{$a[1]}=$count;
                $SEQ{$seq}=$count;
            }
            $Sample{$a[1]}+=$count;
        }
        elsif (m/^@/) {
            my $id=$_;
            #$id=~s/^@//;
            my $seq=<IN>;
            chomp $seq;
            my @a=split(/\_/,$id);
            # >Truseq_F4.RiboZero_ACAGTGA_DJG6PNM1_230_1_1101_18024_2222/1__3828
            my $count=1;
            if ($id=~m/\_\_/) {
                my @b=split(/\_\_/,$id);
                if (scalar(@b) > 1) {$count=$b[-1];}
            }
            if (exists $uniq{$seq}) {
                $uniq{$seq}{$a[1]}+=$count;
                $SEQ{$seq}+=$count;
            }
            else {
                $uniq{$seq}{$a[1]}=$count;
                $SEQ{$seq}=$count;
            }
            $Sample{$a[1]}+=$count;
            $seq=<IN>;
            $seq=<IN>;
        }
    }
}


my $count=0;
my $info="newid";
foreach my $sample (sort keys %Sample) {
    $info=$info."\t".$sample;
}
print OUT2 $info,"\n";

foreach my $seq (sort{$SEQ{$b} <=> $SEQ{$a}} keys %SEQ) {
    $count++;
    my $sum_exp=0;
    my $expr="";
    foreach my $sample (sort keys %Sample) {
        if (exists $uniq{$seq}{$sample}) {
            $sum_exp+=$uniq{$seq}{$sample};
            if ($expr eq "") {$expr=$uniq{$seq}{$sample};}
            else {$expr=$expr."\t".$uniq{$seq}{$sample};}
        }
        else {
            if ($expr eq "") {$expr=0;}
            else {$expr=$expr."\t0";}
        }
    }
    my $id=$header."-".$count."__".$sum_exp;
    print OUT ">".$id,"\n",$seq,"\n";
    print OUT2 $id,"\t",$expr,"\n";
}