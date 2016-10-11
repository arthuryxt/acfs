#!/usr/bin/perl -w
use strict;
use strict;
# check overlapping of known exons border, NOTE that borders must NOT overlap!
die "Usage: $0  \"defined_clusters\"  \"split_exon_gtf\"   \"output_basename\"  \"\(optional\) extend_N_bases==0\"   \"\(optional\) minimum_BS_distance==100\"   \"\(optional\) maximum_BS_distance==1000000\"    \"\(optional\)minimum_SS_Score==10\"   \"\(DEBUG==1\)\"" if (@ARGV < 3);
my $filein=$ARGV[0];   
my $gtf=$ARGV[1];      
my $fileout=$ARGV[2];   
my $Extend=0;           # 0nt by default.
if (scalar(@ARGV) > 3) {$Extend=$ARGV[3];}
my $mingap=100;
if (scalar(@ARGV) > 4) {$mingap=$ARGV[4];}
my $maxgap=1000000;
if (scalar(@ARGV) > 5) {$maxgap=$ARGV[5];}
my $minSSScore=10;
if (scalar(@ARGV) > 6) {$minSSScore=$ARGV[6];}
my $debug=0;
if (scalar(@ARGV) > 7) {$debug=$ARGV[7];}
my $command="rm -f Step4_finished";
system($command);

my %anno;
my %Gene;
my %rawanno;
open IN,$gtf;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[2] eq "exon") {
        if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
        if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
        my $bin1=int($a[3]/1000);
        my $bin2=int($a[4]/1000)+1;
        my @b=split(/\_\_\_/,$a[8]);
        for(my $i=$bin1; $i<=$bin2; $i++) { $anno{$a[0]}{$i}{$a[3]}{$a[4]}=join("\t",$a[8],$a[7],$a[6],$a[5]); $Gene{$a[0]}{$i}{$a[3]}{$a[4]}=$b[0]; }
        $rawanno{$b[0]}{$a[8]}=join("\t",@a);
    }
}
close IN;

open IN,$filein;
open OUT,">".$fileout;
while(<IN>) {
    chomp;
    if (m/^#/) {
        my @a=split("\t",$_);
        $a[4]=join("\t",$a[4],"Left_Exon","Left_N","Left_St","Left_Bt","Left_pos1","Left_pos2","Right_Exon","Right_N","Right_St","Right_Bt","Right_pos1","Right_pos2");
        print OUT join("\t",@a),"\n";
        next;
    }
    my @a=split("\t",$_);
    if ($debug) {print join("\t",@a),"\n";}
    my $info_left="";
    my $info_right="";
    my $bin=int($a[2]/1000);
    my $overlap=999999999999;
    my $f=999999999;
    if (exists $anno{$a[1]}{$bin}) {
        foreach my $start (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin}}) {
            if (($start - $Extend) <= $a[2]) {
                foreach my $end (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin}{$start}}) {
                    if (($end + $Extend) >= $a[2]) {
                        if ($f eq 0) {last;}
                        else {
                            my $gap=($start - $a[2])*($end - $a[2]);
                            if ($gap eq 0) {
                                $f=0;
                                $info_left=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                last;
                            }
                            elsif ($gap > 0) {
                                if ($overlap > $gap) {
                                    $info_left=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                            else {
                                if ($overlap > 0) {
                                    $info_left=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                                elsif (($overlap < 0) and ($gap > $overlap)){
                                    $info_left=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                        }
                    }
                }
            }
        } 
    }
    elsif (exists $anno{$a[1]}{$bin+1}) {
        foreach my $start (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin+1}}) {
            if (($start - $Extend) <= $a[2]) {
                foreach my $end (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin+1}{$start}}) {
                    if (($end + $Extend) >= $a[2]) {
                        if ($f eq 0) {last;}
                        else {
                            my $gap=($start - $a[2])*($end - $a[2]);
                            if ($gap eq 0) {
                                $f=0;
                                $info_left=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                last;
                            }
                            elsif ($gap > 0) {
                                if ($overlap > $gap) {
                                    $info_left=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                            else {
                                if ($overlap > 0) {
                                    $info_left=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                                elsif (($overlap < 0) and ($gap > $overlap)){
                                    $info_left=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                        }
                    }
                }
            }
        } 
        
    }
    
    $bin=int($a[3]/1000);
    $overlap=99999999999;
    $f=999999999;
    if (exists $anno{$a[1]}{$bin}) {
        foreach my $start (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin}}) {
            if (($start - $Extend) <= $a[3])  {
                foreach my $end (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin}{$start}}) {
                    if (($end + $Extend) >= $a[3]) {
                        if ($f eq 0) {last;}
                        else {
                            my $gap=($start - $a[3])*($end - $a[3]);
                            if ($gap eq 0) {
                                $f=0;
                                $info_right=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                last;
                            }
                            elsif ($gap > 0) {
                                if ($overlap > $gap) {
                                    $info_right=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                            else {
                                if ($overlap > 0) {
                                    $info_right=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                                elsif (($overlap < 0) and ($gap > $overlap)){
                                    $info_right=$anno{$a[1]}{$bin}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                        }
                    }
                }
            }
        } 
    }
    elsif (exists $anno{$a[1]}{$bin+1}) {
        foreach my $start (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin+1}}) {
            if (($start - $Extend) <= $a[3])  {
                foreach my $end (sort{$a <=> $b} keys %{$anno{$a[1]}{$bin+1}{$start}}) {
                    if (($end + $Extend) >= $a[3]) {
                        if ($f eq 0) {last;}
                        else {
                            my $gap=($start - $a[3])*($end - $a[3]);
                            if ($gap eq 0) {
                                $f=0;
                                $info_right=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                last;
                            }
                            elsif ($gap > 0) {
                                if ($overlap > $gap) {
                                    $info_right=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                            else {
                                if ($overlap > 0) {
                                    $info_right=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                                elsif (($overlap < 0) and ($gap > $overlap)){
                                    $info_right=$anno{$a[1]}{$bin+1}{$start}{$end}."\t".$start."\t".$end;
                                    $overlap=$gap;
                                }
                            }
                        }
                    }
                }
            }
        } 
        
    }
    my $gap=-1;
    my $gap1=-1;
    my $gap2=-1;
    {
        if ($info_left ne "") {
            my @tmp=split("\t",$info_left);
            $gap1=abs($a[2] - $tmp[-1]);
            $gap+=$gap1;
        }
        if ($info_right ne "") {
            my @tmp=split("\t",$info_right);
            $gap2=abs($a[3] - $tmp[-2]);
            $gap+=$gap2;
        }
        if ($gap > $Extend) {
            if ($gap1 eq 0) { $info_right=join("\t","na","na","na","na","0","0"); }
            elsif ($gap2 eq 0) { $info_left=join("\t","na","na","na","na","0","0"); }
            else {
                $info_left=join("\t","na","na","na","na","0","0");
                $info_right=join("\t","na","na","na","na","0","0");
            }
        }
    }

    if ($info_left eq "") {$info_left=join("\t","na","na","na","na","0","0");}
    if ($info_right eq "") {$info_right=join("\t","na","na","na","na","0","0");}
    $a[4]=$a[4]."\t".$info_left."\t".$info_right;
    print OUT join("\t",@a),"\n";
}
close OUT;

open(IN, $fileout);
open(OUTT1, ">".$fileout."_MEA");
open(OUTT2, ">".$fileout."_CBR");
open(OUTT3, ">".$fileout."_discard");
while (<IN>) {
    chomp;
    next if (m/^#/);
    my @a=split("\t",$_);
    if (($a[4] <= (0-$mingap)) and ($a[4] >= (0-$maxgap)) and ($a[17] >= $minSSScore)) {
        if (($a[6] eq $a[12]) and ($a[6] ne "na")) {
            my @left=split(/\_\_\_/,$a[5]);
            my @right=split(/\_\_\_/,$a[11]);
            my $flag=0;
            if ($left[0] ne $right[0]) {
                if ($debug) { print join("\t",$left[0],$right[0]),"\t",join("\t",@a),"\n"; }
                # make some effort to make left and right annotation on the same gene so that the "_2G" is minimized
                # ENSG00000256407 and ENSG00000215883 both point to the same gene : CYB5RL
                # try to fit the right exon borders to the left Gene_anno first
                foreach my $exon (keys %{$rawanno{$right[0]}}) {
                    if ($flag ne 0) { last; }
                    my @tmp=split("\t",$rawanno{$right[0]}{$exon});
                    my $dist=abs($a[9] - $tmp[3]) + abs($a[10] - $tmp[4]);
                    if ($dist < $Extend) {
                        $flag=$exon;
                    }
                }
                if ($flag eq 0) {
                    # try to fit the left exon borders to the right Gene_anno first
                    foreach my $exon (keys %{$rawanno{$left[0]}}) {
                        if ($flag ne 0) { last; }
                        my @tmp=split("\t",$rawanno{$left[0]}{$exon});
                        my $dist=abs($a[15] - $tmp[3]) + abs($a[16] - $tmp[4]);
                        if ($dist < $Extend) {
                            $flag=$exon;
                        }
                    }
                    if ($flag ne 0) { $a[11]=$flag; }
                }
                else {
                    $a[5]=$flag;
                }
                if ($flag eq 0) { print OUTT2 join("\t",@a),"\n";}
                else { print OUTT1 join("\t",@a),"\n"; }
            }
            else { print OUTT1 join("\t",@a),"\n"; }
        }
        else {
            print OUTT2 join("\t",@a),"\n";
        }
    }
    else {
        print OUTT3 join("\t",@a),"\n";
    }
}

open(OUTFLAG,">Step4_finished");
print OUTFLAG "Step4_finished\n";
close OUTFLAG;


