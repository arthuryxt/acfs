#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"BWA_sam\"    \"output\"   \"\(optional\)min_AS\"   \"\(optional\)cutoff\"" if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my $MAS=30;
if (scalar(@ARGV) > 2) {$MAS=$ARGV[2];}
my $cutoff=0.9;
if (scalar(@ARGV) > 3) {$cutoff=$ARGV[3];}
if (($cutoff <= 0) or ($cutoff >=1 )) { die "cutoff must be within (0,1) !" }
my %uniq;   # store processed ID
open IN, $filein;
open OUT,">".$fileout.".tmp";
open OUT1,">".$fileout.".1";    # single hit
open OUT2,">".$fileout.".multi";
open OUT4,">".$fileout.".segs";
open OUT21,">".$fileout.".2pp"; # two hits on the same chr and same strand
open OUT22,">".$fileout.".2pm"; # two hits on the same chr BUT Different strand
open OUT3,">".$fileout.".unmap";
open OUTMT,">".$fileout.".MT";
open OUT2S,">".$fileout.".2pp.S1"; # selected candidates from two hits on the same chr and same strand
open OUTUID,">".$fileout.".UID";
my $command="rm -f Step1_finished";
system($command);

while(<IN>) {
    chomp;
    if ((m/^#/) or (m/^@/)) {next;}
    my @a=split("\t",$_);
    my $Len=length($a[9]);
    if ($a[1] eq 4) {
        print OUT3 $a[0],"\t",$a[9],"\n";
        print OUTUID $a[0],"\n";
        next;
    }
    if (exists $uniq{$a[0]}) {
        next;
    }
    else {
        $uniq{$a[0]}=1;
        my $strand="+";
        if ($a[1] eq 16) {$strand="-";}
        my $info="";
        my @d=split(/\:/,$a[11]);
        if (scalar(@a) > 14) {
            my @b=split(/\:/,$a[14]);
            my @c=split(/\;/,$b[2]);
            $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]),@c);
        }
        else {
            $info=join("\t",join("\,",$a[2],$strand.$a[3],$a[5],$a[4],$d[2]));
        }
        if ($info =~/MT/) {
            print OUTMT join("\t",$a[0],$Len,$info),"\n";
            next;
        }
        if(scalar(@a) <= 14){
            my @CIGAR_op=($a[5]=~m/[MSID]/g); 
            my @CIGAR_va=($a[5]=~m/\d+/g);
            my $start=0;
            my $length=0;
            my $Nr=scalar(@CIGAR_op);
            if ($strand eq "+") {
                for(my $i=0; $i<$Nr; $i++){
                    if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; }
                    elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; }
                    elsif ($CIGAR_op[$i] eq "S") {
                        if (($start > 0) or ($length > 0)) {last;}
                        $start=$CIGAR_va[$i];
                    }
                }
            }
            else {
                for(my $i=$Nr-1; $i>=0; $i--){
                    if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i];}
                    elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i];}
                    elsif ($CIGAR_op[$i] eq "S") {
                        if (($start > 0) or ($length > 0)) {last;}
                        $start=$CIGAR_va[$i];
                    }
                }
            }
            print OUT1 join("\t",$a[0],$Len,$start,$length,$d[2],$a[2],$strand,$a[3],$a[4]),"\n";
            print OUT join("\t",$a[0],$Len,$info),"\n";
            next;
        }
        print OUT join("\t",$a[0],$Len,$info),"\n";
        # process $info
        my %anno;
        my %Chr;
        my @p=split("\t",$info);
        for(@p) {
            my @tmp=split(/\,/,$_);
            # ignore bad-quality alignments
            if ($tmp[3] < $MAS) {next;}
            my $strandy="+";
            if ($tmp[1]=~m/^-/){
                $strandy="-";
                $tmp[1]=~s/^\-//;
            }
            else {
                $tmp[1]=~s/^\+//;
            }
            my @CIGAR_op=($tmp[2]=~m/[MSID]/g); 
            my @CIGAR_va=($tmp[2]=~m/\d+/g);
            my $start=0;
            my $length=0;
            my $Nr=scalar(@CIGAR_op);
            if ($strandy eq "+") {
                for(my $i=0; $i<$Nr; $i++){
                    if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i]; }
                    elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i]; }
                    elsif ($CIGAR_op[$i] eq "S") {
                        if (($start > 0) or ($length > 0)) {last;}
                        $start=$CIGAR_va[$i];
                    }
                }
            }
            else {
                for(my $i=$Nr-1; $i>=0; $i--){
                    if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i];}
                    elsif ($CIGAR_op[$i] eq "I") {$length+=$CIGAR_va[$i];}
                    elsif ($CIGAR_op[$i] eq "S") {
                        if (($start > 0) or ($length > 0)) {last;}
                        $start=$CIGAR_va[$i];
                    }
                }
            }
            if (exists $anno{$tmp[0]}{$start}) {}
            else {
                $Chr{$tmp[0]}++;
                $anno{$tmp[0]}{$start}=join("\t",$length,$tmp[1],($tmp[1]+$length),$strandy);
            }
        }
        # pick the most-hit chromosome
        my $chromo="";
        foreach my $id(sort{$Chr{$b} <=> $Chr{$a}} keys %Chr) {
            if ($Chr{$id} > 1) {$chromo=$id; last;}
        }
        # each part hit once
        if ($chromo eq "") {
            my $result=$Len;
            my $cov1=0; 
            my $cov2=0; 
            my @Cov;
            for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
            foreach my $id(sort keys %anno) {
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$id}}) {
                    $result=$result."\t".$id."\t".$pos."\t".$anno{$id}{$pos};
                    my @t=split("\t",$anno{$id}{$pos});
                    for(my $i=0; $i<$t[0]; $i++) {$Cov[$i+$pos]=1;}
                }
            }
            my $first=-1;
            my $last=-1;
            for(my $i=0; $i<$Len; $i++) {
                if ($Cov[$i] eq 1) {
                    $cov1++;
                    if ($first eq -1) {$first = $i;}
                    $last=$i;
                }
            }
            $cov2=$last - $first + 1;
            if ($cov1 > 0) { print OUT2 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\n"; }
        }
        # exactly two-hitter, interesting !!!
        elsif ($Chr{$chromo} eq 2) {
            my $flag=1;
            my $gap=0;
            foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                my @tmp=split("\t",$anno{$chromo}{$pos});
                if ($tmp[3] eq "-") {$flag+=2;}
            }
            if ($flag eq 1) {
                # both segments on the plus strand
                my $result=$Len."\t".$chromo;
                my $cov1=0; 
                my $cov2=0; 
                my @Cov;
                for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                    if ($gap eq 0) {$gap = $tmp[2];}
                    else {$gap=$tmp[1] - $gap;}
                }
                my $first=-1;
                my $last=-1;
                for(my $i=0; $i<$Len; $i++) {
                    if ($Cov[$i] eq 1) {
                        $cov1++;
                        if ($first eq -1) {$first = $i;}
                        $last=$i;
                    }
                }
                $cov2=$last - $first + 1;
                print OUT21 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$gap,"\n";
                print OUTUID $a[0],"\n";
                if (($cov2 eq  $cov1) and ($cov2 > $cutoff*$Len)) {
                    print OUT2S $a[0],"\t",$result,"\t",$gap,"\n";
                }
            }
            elsif ($flag eq 5) {
                # both segments on the minus strand
                my $result=$Len."\t".$chromo;
                my $cov1=0; 
                my $cov2=0; 
                my @Cov;
                for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                foreach my $pos (sort{$b <=> $a} keys %{$anno{$chromo}}) {
                    $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                    if ($gap eq 0) {$gap = $tmp[2];}
                    else {$gap=$tmp[1] - $gap;}
                }
                my $first=-1;
                my $last=-1;
                for(my $i=0; $i<$Len; $i++) {
                    if ($Cov[$i] eq 1) {
                        $cov1++;
                        if ($first eq -1) {$first = $i;}
                        $last=$i;
                    }
                }
                $cov2=$last - $first + 1;
                print OUT21 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$gap,"\n";
                print OUTUID $a[0],"\n";
                if (($cov2 eq  $cov1) and ($cov2 > $cutoff*$Len)) {
                    print OUT2S $a[0],"\t",$result,"\t",$gap,"\n";
                }
            }
            else {
                # two segments on different strands
                my $result=$Len."\t".$chromo;
                my $cov1=0; 
                my $cov2=0; 
                my @Cov;
                for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
                foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                    $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                    my @tmp=split("\t",$anno{$chromo}{$pos});
                    for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
                    if ($gap eq 0) {$gap = $tmp[1];}
                    else {$gap=$tmp[1] - $gap;}
                }
                my $first=-1;
                my $last=-1;
                for(my $i=0; $i<$Len; $i++) {
                    if ($Cov[$i] eq 1) {
                        $cov1++;
                        if ($first eq -1) {$first = $i;}
                        $last=$i;
                    }
                }
                $cov2=$last - $first + 1;
                print OUT22 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\t",$gap,"\n";
            }
        }
        # more than two-hitter, whassap ??
        elsif ($Chr{$chromo} > 2) {
            my $result=$Len."\t".$chromo;
            my $cov1=0;
            my $cov2=0; 
            my @Cov;
            for(my $i=0; $i<$Len; $i++) { $Cov[$i]=0; }
            foreach my $pos (sort{$a <=> $b} keys %{$anno{$chromo}}) {
                $result=$result."\t".$pos."\t".$anno{$chromo}{$pos};
                my @tmp=split("\t",$anno{$chromo}{$pos});
                for(my $i=0; $i<$tmp[0]; $i++) {$Cov[$i+$pos]=1;}
            }
            my $first=-1;
            my $last=-1;
            for(my $i=0; $i<$Len; $i++) {
                if ($Cov[$i] eq 1) {
                    $cov1++;
                    if ($first eq -1) {$first = $i;}
                    $last=$i;
                }
            }
            $cov2=$last - $first + 1;
            print OUT4 $a[0],"\t",$cov2,"\t",$cov1,"\t",$result,"\n";
            print OUTUID $a[0],"\n";
        }
    }
}
close IN;


open(OUTFLAG,">Step1_finished");
print OUTFLAG "Step1_finished\n";
close OUTFLAG;
