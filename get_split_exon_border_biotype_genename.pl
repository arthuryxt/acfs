#!/usr/bin/perl -w
use strict;
die "Usage: $0 \"input\"  \"output\"  \"\(optional\)debug\" " if (@ARGV < 2);
my $filein=$ARGV[0];    
my $fileout=$ARGV[1];   
my $debug=0;
if (scalar(@ARGV) > 2) {$debug=$ARGV[2];}
open IN,$filein;
open OUT,">".$fileout;
open OUT2,">".$fileout."_gene";
my %uniq;
my %biotype;
my %Gname;
my %Terminals;
my %Internals;
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    # MT      protein_coding  exon    3307    4262    .       +       .        gene_id "ENSG00000198888"; transcript_id "ENST00000361390"; exon_number "1"; gene_name "MT-ND1"; gene_biotype "protein_coding"; transcript_name "MT-ND1-201";
    if (($a[2] eq "exon") and ($a[0]!~m/NT/)){
        my @b=split(/\"/,$a[8]);
        my $Nr=scalar(@b);
        my $gene_id="";
        my $transcript_id="";
        #my $exon_number="";
        my $gene_name="";
        my $gene_biotype=$a[5];
        for(my $i=0; $i<$Nr; $i=$i+2){
            $b[$i]=~s/\s//g; $b[$i]=~s/\;//g;
            if ($b[$i] eq "gene_id") {$gene_id=$b[$i+1];}
            elsif ($b[$i] eq "transcript_id") {$transcript_id=$b[$i+1];}
            elsif ($b[$i] eq "gene_name") {$gene_name=$b[$i+1];}
            elsif ($b[$i] eq "gene_biotype") {$gene_biotype=$b[$i+1];}
        }
        if ($gene_name eq "") {$gene_name=$gene_id;}
        my $id=$gene_id."\t".$a[0]."\t".$a[6];
        if (!exists $uniq{$id}) {
            $biotype{$id}=$gene_biotype;
            $Gname{$id}=$gene_name;
        }
        $uniq{$id}{$a[3]}{$a[4]}=join("\t",@a);
        if ((exists $Terminals{$id}) and (exists $Terminals{$id}{$transcript_id})) {
            my @T=split("\t",$Terminals{$id}{$transcript_id});
            if ($a[3] < $T[0]) { $T[0]=$a[3]; }
            if ($a[4] > $T[1]) { $T[1]=$a[4]; }
            $Terminals{$id}{$transcript_id}=join("\t",@T);
            $Internals{$id}{$transcript_id}=$Internals{$id}{$transcript_id}."\t".join("\t",$a[3],$a[4]);
        }
        else {
            $Terminals{$id}{$transcript_id}=join("\t",$a[3],$a[4]);
            $Internals{$id}{$transcript_id}=join("\t",$a[3],$a[4]);
        }
    }
}

foreach my $id (keys %uniq) {
    my %INTEx;
    foreach my $tid (keys %{$Internals{$id}}) {
        my @Ttmp=split("\t",$Internals{$id}{$tid});
        my @Stmp=sort {$a <=> $b} @Ttmp;
        my $Nr=scalar(@Stmp);
        for(my $i=1; $i<($Nr-1); $i++) { $INTEx{$Stmp[$i]}=1; }
    }
    my %T;
    my $minT=999999999;
    my $maxT=-1;
    foreach my $tid (keys %{$Terminals{$id}}) {
        my @Ttmp=split("\t",$Terminals{$id}{$tid});
        if (!exists $INTEx{$Ttmp[0]}) { $T{$Ttmp[0]}=$Ttmp[1]; }
        if (!exists $INTEx{$Ttmp[1]}) { $T{$Ttmp[1]}=$Ttmp[0]; }
        if ($minT > $Ttmp[0]) { $minT=$Ttmp[0]; }
        if ($maxT < $Ttmp[1]) { $maxT=$Ttmp[1]; }
    }
    my %tmp;
    my %End;
    foreach my $left (sort {$a <=> $b} keys %{$uniq{$id}}) {
        foreach my $right (sort {$a <=> $b} keys %{$uniq{$id}{$left}}) {
            $tmp{$right}=1;
            $End{$right}=1;
        }
        $tmp{$left}=2;
    }
    
    foreach my $left (sort {$a <=> $b} keys %{$uniq{$id}}) {
        foreach my $right (sort {$a <=> $b} keys %{$uniq{$id}{$left}}) {
            foreach my $test_end (sort {$a <=> $b} keys %End) {
                if (($left < $test_end) and ($test_end < $right)) {$End{$test_end}++;}
            }
        }
    }
    if ($debug eq 1) {
        my $tmp="";
        foreach my $i (sort{$a <=> $b} keys %tmp) {
            if ($tmp eq "") {$tmp=$i."\t".$tmp{$i};}
            else {$tmp=$tmp."\n".$i."\t".$tmp{$i}}
        }
        print "tmp:\n",$tmp,"\n";
        $tmp="";
        foreach my $i (sort{$a <=> $b} keys %End) {
            if ($tmp eq "") {$tmp=$i."\t".$End{$i};}
            else {$tmp=$tmp."\n".$i."\t".$End{$i}}
        }
        print "End:\n",$tmp,"\n";
        $tmp="";
        foreach my $i (sort{$a <=> $b} keys %T) {
            if ($tmp eq "") {$tmp=$i;}
            else {$tmp=$tmp."\n".$i;}
        }
        print "Terminal:\n",$tmp,"\n";
    }
    
    my @border=();
    my $Nr=0;
    my $last=-1;
    my $usedrightborder=-1;
    my %used;
    foreach my $pos (sort{$a <=> $b} keys %tmp) {
        if ($debug eq 1) { print join("\t","debug_loop_start",$Nr,$last,$pos),"\n"; }
        if ($last eq -1) { $last=$pos; }
        #elsif ($usedrightborder eq $last) { $last=$pos; }
        else {
            if (exists $End{$last}) {   # $last could be a right_border, therefore $pos could be a left_border
                if ($debug eq 1) {print $tmp{$last},"\t",$End{$last},"\t",$tmp{$pos},"\n";}
                if ($End{$last} > 1) {  # $last is middle_border
                    if ((exists $T{$pos}) and ($pos ne $maxT) and ($pos ne $minT)) {
                        #$last=$usedrightborder;
                    }
                    elsif ((exists $T{$last}) and ($last ne $maxT) and ($last ne $minT)) {
                        $last=$pos;
                    }
                    elsif ($tmp{$pos} eq 2) {  # $pos is a left_border
                        #if ($usedrightborder eq $last) {
                        #}
                        #else{
                            my $p5=($last+1);
                            my $p3=($pos-1);
                            if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; $usedrightborder=$p3; }
                        #}
                        $last=$pos;
                    }
                    else { 
                        my $p5=($last+1);
                        my $p3=($pos);
                        if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; $usedrightborder=$p3; }
                        $last=$pos;
                    }
                }
                else { $last=$pos; }
            }
            else {  # $last is a left_border; and $pos could be a right_border
                if ($debug eq 1) {print $tmp{$last},"\tNA\t",$tmp{$pos},"\n";}
                if ((exists $T{$pos}) and ($pos ne $maxT) and ($pos ne $minT)) {  # $pos is a transcript terminal, and should not serve as a splicing site
                    # don't update $last, find a new $pos,
                }
                elsif ((exists $T{$last}) and ($last ne $maxT) and ($last ne $minT)) {
                    $last=$pos;
                }
                elsif ($tmp{$pos} eq 2) {  # $pos is a left_border, minus 1 to prepare for the adjacent exon
                    my $p5=($last);
                    my $p3=($pos-1);
                    if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; $usedrightborder=$p3; }
                    $last=$pos;
                } 
                else {  # $pos is a right_border
                    my $p5=($last);
                    my $p3=($pos);
                    if ($p5 <= $p3) { $border[$Nr]=$p5."\t".$p3; $Nr++; $usedrightborder=$p3; }
                    $last=$pos;
                }
            }
        }
        if (($debug eq 1) and ($Nr > 0)){
            my $tmp=$border[0];
            for(my $i=1; $i<$Nr; $i++) { $tmp=$tmp."\t".$border[$i];}
            print $tmp,"\n debug_loop_end \n";
        }
    }
    if ($debug eq 1) {print $tmp{$last},"\t",$End{$last},"\n";}
    if (($debug eq 2) and (!exists $End{$last})) {
        foreach my $left (sort {$a <=> $b} keys %{$uniq{$id}}) {
            foreach my $right (sort {$a <=> $b} keys %{$uniq{$id}{$left}}) {
                print join("\t",$id,$left,$right,$uniq{$id}{$left}{$right}),"\n";
            }
        }
    }
    if (($tmp{$last} eq 2) and ($End{$last} eq 1)) { $border[$Nr]=$last."\t".$last; $Nr++; }
    
    my @a=split("\t",$id);
    my $start_g=-1;
    my $end_g=-1;
    my $biotype_g="";
    my $gene_id="";
    my $chr="";
    my $strand="";
    if ($a[2] eq "+") {
        for (my $i=0; $i<$Nr; $i++) {
            print OUT join("\t",$a[1],"split","exon",$border[$i],$biotype{$id},$a[2],$Gname{$id},$a[0]."___".($i+1)."___".$Nr),"\n";
            my @b=split("\t",$border[$i]);
            if ($start_g eq -1) {$start_g=$b[0]; $chr=$a[1]; $biotype_g=$biotype{$id}; $gene_id=$a[0]; $strand=$a[2];}
            if ($end_g < $b[1]) {$end_g=$b[1];}
        }
    }
    else {
        for (my $i=0; $i<$Nr; $i++) {
            print OUT join("\t",$a[1],"split","exon",$border[$i],$biotype{$id},$a[2],$Gname{$id},$a[0]."___".($Nr - $i)."___".$Nr),"\n";
            my @b=split("\t",$border[$i]);
            if ($start_g eq -1) {$start_g=$b[0]; $chr=$a[1]; $biotype_g=$biotype{$id}; $gene_id=$a[0]; $strand=$a[2];}
            if ($end_g < $b[1]) {$end_g=$b[1];}
        }
    }
    print OUT2 join("\t",$chr,".","gene",$start_g,$end_g,$biotype_g,$strand,$Gname{$id},$gene_id),"\n";
}
