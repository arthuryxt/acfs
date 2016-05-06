#!/usr/bin/perl -w
use strict;
# check splicing sites for 2-segment-diff_chr_strand
die "Usage: $0    \"unmap.parsed.multi.parsed\"  \"genome_location\"  \"output basename\"  \"\(optional\)read_Length==100\"  \"\(optional\)min_Score==10\"  \"\(optional\)agtf\"  \"\(DEBUG set to 1\)\"" if (@ARGV < 3);
my $anno=$ARGV[0];      # unmap.parsed.multi.parsed
my $genome=$ARGV[1];    # /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes/
my $fileout=$ARGV[2];   # unmap.parsed.multi.parsed.tsloci
my $readLen=100;
if (scalar(@ARGV) > 3) {$readLen=$ARGV[3];}
my $minSS=10;
if (scalar(@ARGV) > 4) {$minSS=$ARGV[4];}
my $agtf="";            # /home/xintian/genome/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf
if (scalar(@ARGV) > 5) {$agtf=$ARGV[5];}
my $debug=0;
if (scalar(@ARGV) > 6) {$debug=$ARGV[6];}

my %Anno;
my %Gene;
my $bin=1000;
if ($agtf ne "") {
    open IN,$agtf;
    while(<IN>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[2] eq "exon") {
            if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
            if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
            my $bin1=int($a[3]/$bin);
            my $bin2=int($a[4]/$bin)+1;
            my @b=split(/\_\_\_/,$a[8]);
            for(my $i=$bin1; $i<=$bin2; $i++) { $Anno{$a[0]}{$i}{$a[3]}{$a[4]}=join("\t",$a[8],$a[7],$a[6],$a[5]); $Gene{$a[0]}{$i}{$a[3]}{$a[4]}=$b[0]; }
        }
    }
    close IN;
}

open IN0,$anno.".01";
my %uniq;
my %reads;
my %CHR;
while (<IN0>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[16] < $minSS) { next; }
    # 0             1   2   3   4   5   6   7           8           9   10  11  12  13          14          15  16      17      18      19  20  21  22  23  24
    # newid-12__8   100 100 100 7   0   60  44246987    44247046    -   6   60  40  44200126    44200165    -   20.28   10.24   10.04   +   4   3   1   3   3
    # AGCAGGGGATCCTACTGGCCAGTCTATCCTGTCGACTTGCTTGGAGAATTCATCTAGTACCCACATGAATACAGCTGTGAGGCTCCGGCCCAACCAGTCA
    if (($a[9] eq "-") and ($a[15] eq "-")) {
        my $id=join("_",$a[4],$a[7],$a[9],$a[10],$a[14],$a[15]);
        $CHR{$a[4]}=1;
        $CHR{$a[10]}=1;
        if (exists $uniq{$id}) {
            my @b=split("\t",$uniq{$id});
            if ($a[8] > $b[8]) { $b[8]=$a[8]; }
            if ($a[13] < $b[13]) { $b[13]=$a[13]; }
            $reads{$id}=$reads{$id}.";".$a[0];
            $uniq{$id}=join("\t",@b);
        }
        else { $uniq{$id}=join("\t",@a); $reads{$id}=$a[0];}
    }
    # 0             1   2   3   4   5   6   7           8           9   10  11  12  13          14          15  16      17      18      19  20  21  22  23  24
    # newid-2__9	100	100	100	2	0	65	30379494	30379558	-	12	65	35	10570943	10570977	+	19.90	9.48	10.42	+	2
    # AGACGGGTACCACCGATATGATCAAGGAAAATTCTGCCCATTTTTATGGCTGAAGTTCTAAAAACCATTTCTTCTTCATTATCTATACAAAGCAGACTAG
    elsif (($a[9] eq "-") and ($a[15] eq "+")) {
        my $id=join("_",$a[4],$a[7],$a[9],$a[10],$a[13],$a[15]);
        $CHR{$a[4]}=1;
        $CHR{$a[10]}=1;
        if (exists $uniq{$id}) {
            my @b=split("\t",$uniq{$id});
            if ($a[8] > $b[8]) { $b[8]=$a[8]; }
            if ($a[14] > $b[14]) { $b[14]=$a[14]; }
            $reads{$id}=$reads{$id}.";".$a[0];
            $uniq{$id}=join("\t",@b);
        }
        else { $uniq{$id}=join("\t",@a); $reads{$id}=$a[0];}
    }
    # 0             1   2   3   4   5   6   7           8           9   10  11  12  13          14          15  16      17      18      19  20  21  22  23  24
    # newid-36__6	100	100	100	11	0	50	119185924	119185973	+	X	50	50	107396268	107396317	-	22.56	9.88	12.68	-	2
    # TTCCACCTCCACCAGCTCAGGCGCAGGCTGCTCAGCCTCTCCGGGCACACCTGCCCCAGTGGTTGTCACTTCTGGCTTGGCTGGAGGTACAAAGGGAGGC
    elsif (($a[9] eq "+") and ($a[15] eq "-")) {
        my $id=join("_",$a[4],$a[8],$a[9],$a[10],$a[14],$a[15]);
        $CHR{$a[4]}=1;
        $CHR{$a[10]}=1;
        if (exists $uniq{$id}) {
            my @b=split("\t",$uniq{$id});
            if ($a[7] < $b[7]) { $b[7]=$a[7]; }
            if ($a[13] < $b[13]) { $b[13]=$a[13]; }
            $reads{$id}=$reads{$id}.";".$a[0];
            $uniq{$id}=join("\t",@b);
        }
        else { $uniq{$id}=join("\t",@a); $reads{$id}=$a[0];}
    }
    # 0             1   2   3   4   5   6   7           8           9   10  11  12  13          14          15  16      17      18      19  20  21  22  23  24
    # newid-11__8	100	100	100	17	0	60	40971568	40971627	+	12	60	40	49167252	49167291	+	17.76	7.15	10.61	-	1
     # GGCTGGGGGGATGAATCTGCGAGAGACACCATCCTGGCGAGGAGTTTCAATAAATGGCTCAAGCCATGGACGCCAAGCAGTAGGTCATAGTTGTCAAAGA
    elsif (($a[9] eq "+") and ($a[15] eq "+")) {
        my $id=join("_",$a[4],$a[8],$a[9],$a[10],$a[13],$a[15]);
        $CHR{$a[4]}=1;
        $CHR{$a[10]}=1;
        if (exists $uniq{$id}) {
            my @b=split("\t",$uniq{$id});
            if ($a[7] < $b[7]) { $b[7]=$a[7]; }
            if ($a[14] > $b[14]) { $b[14]=$a[14]; }
            $reads{$id}=$reads{$id}.";".$a[0];
            $uniq{$id}=join("\t",@b);
        }
        else { $uniq{$id}=join("\t",@a); $reads{$id}=$a[0];}
    }
}
close IN0;

my %reported;
my %genome;
my $f=0; 
foreach my $chr (sort keys %CHR) {
    #if (length($chr) > 3) {next}
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    print "searching chromosome : ".$chr.".fa";
    $f=-1;
    if (exists $reported{$chr}) {next;}
    open (IN1, $genome."/".$chr.".fa") or next;
    <IN1>;
    my $SEQ;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $SEQ=$SEQ.$_;}
	$genome{$chr}=$SEQ;
    close IN1;
}
if ($f eq 1) {print "\t=======\tyes\n";}
elsif ($f eq -1) {print "\t=======\tno\n";}
if ($f eq -1) {
    foreach my $chr (sort keys %CHR) {
    #if (length($chr) > 3) {next}
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    print "searching chromosome : chr".$chr.".fa";
    $f=-1;
    if (exists $reported{$chr}) {next;}
    open (IN1, $genome."/chr".$chr.".fa") or next;
    <IN1>;
    my $SEQ;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $SEQ=$SEQ.$_;}
	$genome{$chr}=$SEQ;
    close IN1;
}
if ($f eq 1) {print "\t=======\tyes\n";}
elsif ($f eq -1) {print "\t=======\tno\n";}
}


open(OUT1, ">".$fileout.".tsloci");
open(OUT2, ">".$fileout.".tsloci.fa");
foreach my $id (keys %uniq) {
    my @a=split("\t",$uniq{$id});
    $a[0]=$id;
    my $seq="";
    print OUT1 join("\t",@a),"\t",$reads{$id},"\n";
    if (($a[9] eq "-") and ($a[15] eq "-")) {
        # AGCAGGGGATCCTACTGGCCAGTCTATCCTGTCGACTTGCTTGGAGAATTCATCTAGTACCCACATGAATACAGCTGTGAGGCTCCGGCCCAACCAGTCA
        my $tmp1=substr($genome{$a[4]},$a[7]-1,$readLen);
        my $tmp2=substr($genome{$a[10]},($a[14]-$readLen),$readLen);
        $seq=$tmp2.$tmp1;
    }
    elsif (($a[9] eq "-") and ($a[15] eq "+")) {
        # AGACGGGTACCACCGATATGATCAAGGAAAATTCTGCCCATTTTTATGGCTGAAGTTCTAAAAACCATTTCTTCTTCATTATCTATACAAAGCAGACTAG
        my $tmp1=substr($genome{$a[4]},$a[7]-1,$readLen);
        my $tmp2=substr($genome{$a[10]},$a[13]-1,$readLen);
        $tmp2=~tr/[atcgATCG]/[TAGCTAGC]/;
        my $rc2=scalar reverse $tmp2;
        $seq=$rc2.$tmp1;
    }
    elsif (($a[9] eq "+") and ($a[15] eq "-")) {
        # TTCCACCTCCACCAGCTCAGGCGCAGGCTGCTCAGCCTCTCCGGGCACACCTGCCCCAGTGGTTGTCACTTCTGGCTTGGCTGGAGGTACAAAGGGAGGC
        my $tmp1=substr($genome{$a[4]},($a[8]-$readLen),$readLen);
        my $tmp2=substr($genome{$a[10]},($a[14]-$readLen),$readLen);
        $tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
        my $rc1=scalar reverse $tmp1;
        $seq=$tmp2.$rc1;
    }
    elsif (($a[9] eq "+") and ($a[15] eq "+")) {
        # GGCTGGGGGGATGAATCTGCGAGAGACACCATCCTGGCGAGGAGTTTCAATAAATGGCTCAAGCCATGGACGCCAAGCAGTAGGTCATAGTTGTCAAAGA
        my $tmp1=substr($genome{$a[4]},($a[8]-$readLen),$readLen);
        my $tmp2=substr($genome{$a[10]},($a[13]-1),$readLen);
        $tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
        my $rc1=scalar reverse $tmp1;
        $tmp2=~tr/[atcgATCG]/[TAGCTAGC]/;
        my $rc2=scalar reverse $tmp2;
        $seq=$rc2.$rc1;
    }
    print OUT2 ">".$id,"\n",$seq,"\n";
}
close OUT1;
close OUT2;

# try to add annotation to the trans-splicing events
if ($agtf ne "") {
    open(OUT3, ">".$fileout.".tsloci.anno");
    foreach my $id (keys %uniq) {
        my @a=split("\t",$uniq{$id});
        $a[0]=$id;
        if (($a[9] eq "-") and ($a[15] eq "-")) {
            my $leftBorder=$a[7];
            my $rightBorder=$a[14];
            my $leftChr=$a[4];
            my $rightChr=$a[10];
            my $leftStr=$a[9];
            my $rightStr=$a[15];
            my $bin1=int($leftBorder/$bin);
            my $left="";
            my $bin2=int($rightBorder/$bin);
            my $right="";
            foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                    my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                    if (($b[2] ne $leftStr) and ($leftBorder eq $posS)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                }
                if ($left ne "") { last;}
            }
            if ($left eq ""){
                foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                    foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                        my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                        if (($b[2] ne $leftStr) and ($leftBorder eq $posS)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                    }
                    if ($left ne "") { last;}
                }
            }
            foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                    my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                    if (($b[2] ne $rightStr) and ($rightBorder eq $posE)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                }
                if ($right ne "") { last;}
            }
            if ($right eq ""){
                foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                    foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                        my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                        if (($b[2] ne $rightStr) and ($rightBorder eq $posE)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                    }
                    if ($right ne "") { last;}
                }
            }
            if (($left ne "") and ($right ne "")) {
                my @L=split("\t",$left);
                my @R=split("\t",$right);
                my $gL=$L[3];
                my $gR=$R[3];
                print OUT3 join("\t",$id,"yes",$gL,$gR,$left,$right),"\n";
            }
            elsif (($left ne "") and ($right eq "")) {
                my @L=split("\t",$left);
                my $gL=$L[3];
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no",$gL,"NA",$left,$right),"\n";
            }
            elsif (($left eq "") and ($right ne "")) {
                my @R=split("\t",$right);
                my $gR=$R[3];
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA",$gR,$left,$right),"\n";
            }
            else{
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA","NA",$left,$right),"\n";
            } 
        }
        elsif (($a[9] eq "-") and ($a[15] eq "+")) {
            my $leftBorder=$a[7];
            my $rightBorder=$a[13];
            my $leftChr=$a[4];
            my $rightChr=$a[10];
            my $leftStr=$a[9];
            my $rightStr=$a[15];
            my $bin1=int($leftBorder/$bin);
            my $left="";
            my $bin2=int($rightBorder/$bin);
            my $right="";
            foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                    my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                    if (($b[2] ne $leftStr) and ($leftBorder eq $posS)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                }
                if ($left ne "") { last;}
            }
            if ($left eq ""){
                foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                    foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                        my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                        if (($b[2] ne $leftStr) and ($leftBorder eq $posS)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                    }
                    if ($left ne "") { last;}
                }
            }
            foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                    my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                    if (($b[2] ne $rightStr) and ($rightBorder eq $posS)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                }
                if ($right ne "") { last;}
            }
            if ($right eq ""){
                foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                    foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                        my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                        if (($b[2] ne $rightStr) and ($rightBorder eq $posS)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                    }
                    if ($right ne "") { last;}
                }
            }
            if (($left ne "") and ($right ne "")) {
                my @L=split("\t",$left);
                my @R=split("\t",$right);
                my $gL=$L[3];
                my $gR=$R[3];
                print OUT3 join("\t",$id,"yes",$gL,$gR,$left,$right),"\n";
            }
            elsif (($left ne "") and ($right eq "")) {
                my @L=split("\t",$left);
                my $gL=$L[3];
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no",$gL,"NA",$left,$right),"\n";
            }
            elsif (($left eq "") and ($right ne "")) {
                my @R=split("\t",$right);
                my $gR=$R[3];
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA",$gR,$left,$right),"\n";
            }
            else{
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA","NA",$left,$right),"\n";
            } 
        }
        elsif (($a[9] eq "+") and ($a[15] eq "-")) {
            my $leftBorder=$a[8];
            my $rightBorder=$a[14];
            my $leftChr=$a[4];
            my $rightChr=$a[10];
            my $leftStr=$a[9];
            my $rightStr=$a[15];
            my $bin1=int($leftBorder/$bin);
            my $left="";
            my $bin2=int($rightBorder/$bin);
            my $right="";
            foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                    my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                    if (($b[2] ne $leftStr) and ($leftBorder eq $posE)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                }
                if ($left ne "") { last;}
            }
            if ($left eq ""){
                foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                    foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                        my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                        if (($b[2] ne $leftStr) and ($leftBorder eq $posE)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                    }
                    if ($left ne "") { last;}
                }
            }
            foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                    my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                    if (($b[2] ne $rightStr) and ($rightBorder eq $posE)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                }
                if ($right ne "") { last;}
            }
            if ($right eq ""){
                foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                    foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                        my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                        if (($b[2] ne $rightStr) and ($rightBorder eq $posE)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                    }
                    if ($right ne "") { last;}
                }
            }
            if (($left ne "") and ($right ne "")) {
                my @L=split("\t",$left);
                my @R=split("\t",$right);
                my $gL=$L[3];
                my $gR=$R[3];
                print OUT3 join("\t",$id,"yes",$gL,$gR,$left,$right),"\n";
            }
            elsif (($left ne "") and ($right eq "")) {
                my @L=split("\t",$left);
                my $gL=$L[3];
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no",$gL,"NA",$left,$right),"\n";
            }
            elsif (($left eq "") and ($right ne "")) {
                my @R=split("\t",$right);
                my $gR=$R[3];
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA",$gR,$left,$right),"\n";
            }
            else{
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA","NA",$left,$right),"\n";
            } 
        }
        elsif (($a[9] eq "+") and ($a[15] eq "+")) {
            my $leftBorder=$a[8];
            my $rightBorder=$a[13];
            my $leftChr=$a[4];
            my $rightChr=$a[10];
            my $leftStr=$a[9];
            my $rightStr=$a[15];
            my $bin1=int($leftBorder/$bin);
            my $left="";
            my $bin2=int($rightBorder/$bin);
            my $right="";
            foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                    my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                    if (($b[2] ne $leftStr) and ($leftBorder eq $posE)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                }
                if ($left ne "") { last;}
            }
            if ($left eq ""){
                foreach my $posS (keys %{$Anno{$leftChr}{$bin1}}) {
                    foreach my $posE (keys %{$Anno{$leftChr}{$bin1}{$posS}}) {
                        my @b=split("\t",$Anno{$leftChr}{$bin1}{$posS}{$posE});
                        if (($b[2] ne $leftStr) and ($leftBorder eq $posE)) { $left=join("\t",$posS,$posE,$Anno{$leftChr}{$bin1}{$posS}{$posE}); last;}
                    }
                    if ($left ne "") { last;}
                }
            }
            foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                    my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                    if (($b[2] ne $rightStr) and ($rightBorder eq $posS)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                }
                if ($right ne "") { last;}
            }
            if ($right eq ""){
                foreach my $posS (keys %{$Anno{$rightChr}{$bin2}}) {
                    foreach my $posE (keys %{$Anno{$rightChr}{$bin2}{$posS}}) {
                        my @b=split("\t",$Anno{$rightChr}{$bin2}{$posS}{$posE});
                        if (($b[2] ne $rightStr) and ($rightBorder eq $posS)) { $right=join("\t",$posS,$posE,$Anno{$rightChr}{$bin2}{$posS}{$posE}); last; }
                    }
                    if ($right ne "") { last;}
                }
            }
            if (($left ne "") and ($right ne "")) {
                my @L=split("\t",$left);
                my @R=split("\t",$right);
                my $gL=$L[3];
                my $gR=$R[3];
                print OUT3 join("\t",$id,"yes",$gL,$gR,$left,$right),"\n";
            }
            elsif (($left ne "") and ($right eq "")) {
                my @L=split("\t",$left);
                my $gL=$L[3];
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no",$gL,"NA",$left,$right),"\n";
            }
            elsif (($left eq "") and ($right ne "")) {
                my @R=split("\t",$right);
                my $gR=$R[3];
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA",$gR,$left,$right),"\n";
            }
            else{
                $left=join("\t",$a[7],$a[8],"NA","NA","NA","NA");
                $right=join("\t",$a[13],$a[14],"NA","NA","NA","NA");
                print OUT3 join("\t",$id,"no","NA","NA",$left,$right),"\n";
            } 
        }
    }
}