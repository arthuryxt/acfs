#!/usr/bin/perl -w
use strict;
# check splicing sites for 2-segment-diff_chr_strand
die "Usage: $0   \"CB_splice_folder\"  \"unmap.parsed.tmp\"  \"genome_location\"  \"output basename\"  \"\(optional\)cutoff\"  \"\(optional\) extend N bases\"  \"\(optional\)min_AS\"  \"\(optional\)maxJump\"   \"\(DEBUG set to 1\)\"" if (@ARGV < 4);
my $DIR=$ARGV[0];		# /home/CB_splice/
my $anno=$ARGV[1];      # unmap.parsed.tmp
my $genome=$ARGV[2];    # /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes/
my $fileout=$ARGV[3];
my $cutoff=0.9;
if (scalar(@ARGV) > 4) {$cutoff=$ARGV[4];}
my $Extend=15;           # 15nt by default. as 3' splice strength need 23nt, 20nt from intron and 3 from exon.
if (scalar(@ARGV) > 5) {$Extend=$ARGV[5];}
my $MAS=0;
if (scalar(@ARGV) > 6) {$MAS=$ARGV[6];}
my $maxJump=1000000;
if (scalar(@ARGV) > 7) {$maxJump=$ARGV[7];}
my $debug=0;
if (scalar(@ARGV) > 8) {$debug=$ARGV[8];}
die "Usage: $0   \"CB_splice_folder\"  \"unmap.parsed.tmp\"  \"genome_location\"  \"output basename\"  \"\(optional\)cutoff\"  \"\(optional\) extend N bases\"  \"\(optional\)min_AS\"  \"\(optional\)maxJump\"   \"\(DEBUG set to 1\)\"" if (scalar(@ARGV) > 9);

my %uniq;
my %READ;
my %genome;
open IN,$anno;
open OUT,">".$fileout;
open OUT1,">".$fileout.".0";	# splice-motif found				trust hit  CB prediction
open OUT2,">".$fileout.".1";	# closure == 0,						trust best CB prediction
open OUT3,">".$fileout.".01";	# combined files of the above 2
open OUTerr,">".$fileout."_withN";
open OUTp,">".$fileout.".S1";
while(<IN>) {
    chomp;
    # modified to read from "unmap.parsed.tmp" directly, since the trans-splicing overhand could contain more low complexity sequences and hense fail to pass the MAS critera
    # the MAS can be re-enact here
    my @a0=split("\t",$_);
    # newid-1__1	75	11,-118352428,36S39M,48,0	9,+20365667,37S38M,48,0
    if (scalar(@a0) ne 4) { next; }
    my $Len=$a0[1];
    my @a2=split(/,/,$a0[2]);
    my @a3=split(/,/,$a0[3]);
    my $strand1=substr($a2[1],0,1);
    my $strand2=substr($a3[1],0,1);
    my $pos1=substr($a2[1],1);
    my $pos2=substr($a3[1],1);
    my $multi="";
    if (($a2[0] ne $a3[0]) or ($strand1 ne $strand2) or (abs($pos1 - $pos2) > $maxJump)) {
        if (($a2[3] + $a3[3]) >= (2 * $MAS)) {
            my $info=$a0[2]."\t".$a0[3];
            # process $info
            my %anno;
            my %Chr;
            my @p=split("\t",$info);
            for(@p) {
                my @tmp=split(/\,/,$_);
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
                        elsif ($CIGAR_op[$i] eq "D") {$length+=$CIGAR_va[$i]; }
                        elsif ($CIGAR_op[$i] eq "I") {$length-=$CIGAR_va[$i]; }
                        elsif ($CIGAR_op[$i] eq "S") {
                            if (($start > 0) or ($length > 0)) {last;}
                            $start=$CIGAR_va[$i];
                        }
                    }
                }
                else {
                    for(my $i=$Nr-1; $i>=0; $i--){
                        if ($CIGAR_op[$i] eq "M") {$length+=$CIGAR_va[$i];}
                        elsif ($CIGAR_op[$i] eq "D") {$length+=$CIGAR_va[$i]; }
                        elsif ($CIGAR_op[$i] eq "I") {$length-=$CIGAR_va[$i]; }
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
                if ($cov1 > 0) { $multi=$a0[0]."\t".$cov2."\t".$cov1."\t".$result; }
            } else {
                # hit on the same chromosome
                # to-be added
            }
            
        }
    }
    if ($multi eq "") { next; }
    my @a=split("\t",$multi);
    
    # remove the above if start from "unmap.parsed.multi"
    # and enable the following line:
    # my @a=split("\t",$_);
    # 0                       1       2       3       4       5			6       7    		8    		9     	10      11		12		13			14			15
    #						  in_len  out_len len     chr     seq_s		len     pos_s	  	pos_e	 	str     chr     seq_s   len     pos_s		pos_e		strand
    # newid-231525/1__69      151     151     151     3       0			115     144260022	144260137	-       X       113     38      108011492	108011530	-
    # newid-1__1	          75	  75	  75	  11	  0	        39	    118352428	118352467	-	    9	    37	    38	    20365667	20365705	+
    if(($a[1] eq $a[2]) and ($a[1] > $cutoff*$a[3]) and (scalar(@a) eq 16)){
		#print OUTp join("\t",@a),"\n";
		$a[8]=$a[8]-1;
		$a[14]=$a[14]-1;
		if ($a[4]=~m/^chromosome/i) {$a[4]=~s/chromosome//i;}
		    if ($a[4]=~m/^chr/i) {$a[4]=~s/chr//i;}
		if ($a[10]=~m/^chromosome/i) {$a[10]=~s/chromosome//i;}
		    if ($a[10]=~m/^chr/i) {$a[10]=~s/chr//i;}
		if ($a[5] < $a[11]) {
			$uniq{$a[4]}++;
			$uniq{$a[10]}++;
			$READ{$a[0]}=join("\t",@a);
			print OUTp $READ{$a[0]},"\n";
		    #$uniq{$a[4]}{$a[7]}{$a[8]}{$a[0]}=join("\t",1,$a[5],$a[6],$a[9]);
		    #$uniq{$a[10]}{$a[13]}{$a[14]}{$a[0]}=join("\t",2,$a[11],$a[12],$a[15]);
		}
		else {
			$uniq{$a[4]}++;
			$uniq{$a[10]}++;
			#$READ{$a[0]}=join("\t",@a);
			$READ{$a[0]}=join("\t",$a[0],$a[1],$a[2],$a[3],$a[10],$a[11],$a[12],$a[13],$a[14],$a[15],$a[4],$a[5],$a[6],$a[7],$a[8],$a[9]);
			print OUTp $READ{$a[0]},"\n";
		    #$uniq{$a[4]}{$a[7]}{$a[8]}{$a[0]}=join("\t",2,$a[5],$a[6],$a[9]);
		    #$uniq{$a[10]}{$a[13]}{$a[14]}{$a[0]}=join("\t",1,$a[11],$a[12],$a[15]);
		}
		#$READ{$a[0]}=join("\t",@a);
    }
}
close IN;

##### calculate 5p splice strength #####

my $tmpfile1=$DIR."/me2x5";
#my %me2x5 = &makescorematrix("~/CB_splice/me2x5");
my %me2x5 = &makescorematrix($tmpfile1);
my $tmpfile2=$DIR."/splicemodels/splice5sequences";
#my %seq = &makesequencematrix("/home/CB_splice/splicemodels/splice5sequences");
my %seq = &makesequencematrix($tmpfile2);

my %bgd;
$bgd{'A'} = 0.27;
$bgd{'C'} = 0.23;
$bgd{'G'} = 0.23;
$bgd{'T'} = 0.27; 
  
sub makesequencematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$_} = $n;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub makescorematrix{
    my $file = shift;
    my %matrix;my $n=0;
    open(SCOREF, $file) || die "Can't open $file!\n";
    while(<SCOREF>) { 
	chomp;
	$_=~ s/\s//;
	$matrix{$n} = $_;
	$n++;
    }
    close(SCOREF);
    return %matrix;
}

sub getrest5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    return $seqa[0].$seqa[1].$seqa[2].$seqa[5].$seqa[6].$seqa[7].$seqa[8];
}

sub scoreconsensus5{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.004; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.9896; $cons1{'T'} = 0.0032;
    my %cons2;
    $cons2{'A'} = 0.0034; $cons2{'C'} = 0.0039; $cons2{'G'} = 0.0042; $cons2{'T'} = 0.9884;
    my $addscore = $cons1{$seqa[3]}*$cons2{$seqa[4]}/($bgd{$seqa[3]}*$bgd{$seqa[4]}); 
    return $addscore;
}

sub log2{
    my ($val) = @_;
    return log($val)/log(2);
}
##### calculate 5p splice strength #####

##### calculate 3p splice strength #####
sub hashseq{
    #returns hash of sequence in base 4
    # &hashseq('CAGAAGT') returns 4619
    my $seq = shift;
    $seq = uc($seq);
    $seq =~ tr/ACGT/0123/;
    my @seqa = split(//,$seq);
    my $sum = 0;
    my $len = length($seq);
    my @four = (1,4,16,64,256,1024,4096,16384);
    my $i=0;
    while ($i<$len) {
        $sum+= $seqa[$i] * $four[$len - $i -1] ;
	$i++;
    }
    return $sum;
}

my @metables = &makemaxentscores;
sub makemaxentscores{
    my $dir = $DIR."/splicemodels/";
    my @list = ('me2x3acc1','me2x3acc2','me2x3acc3','me2x3acc4','me2x3acc5','me2x3acc6','me2x3acc7','me2x3acc8','me2x3acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makewmmscores{
    my $dir = $DIR."/splicemodels/";
    my @list = ('me1s0acc1','me1s0acc2','me1s0acc3','me1s0acc4','me1s0acc5','me1s0acc6','me1s0acc7','me1s0acc8','me1s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub makemmscores{
    my $dir = $DIR."/splicemodels/";
    my @list = ('me2s0acc1','me2s0acc2','me2s0acc3','me2s0acc4','me2s0acc5','me2s0acc6','me2s0acc7','me2s0acc8','me2s0acc9');
    my @metables;
    my $num = 0 ;
    foreach my $file (@list) {
	my $n = 0;
	open (SCOREF,"<".$dir.$file) || die "Can't open $file!\n";
	while(<SCOREF>) {
	    chomp;
	    $_=~ s/\s//;
	    $metables[$num]{$n} = $_;
	    $n++;
	}
	close(SCOREF);
	#print STDERR $file."\t".$num."\t".$n."\n";
	$num++;
    }
    return @metables;
}
sub maxentscore{
    my $seq = shift;
    my $table_ref = shift;
    my @metables = @$table_ref;
    my @sc;
    $sc[0] = $metables[0]{&hashseq(substr($seq,0,7))};
    $sc[1] = $metables[1]{&hashseq(substr($seq,7,7))};
    $sc[2] = $metables[2]{&hashseq(substr($seq,14,7))};
    $sc[3] = $metables[3]{&hashseq(substr($seq,4,7))};
    $sc[4] = $metables[4]{&hashseq(substr($seq,11,7))};
    $sc[5] = $metables[5]{&hashseq(substr($seq,4,3))};
    $sc[6] = $metables[6]{&hashseq(substr($seq,7,4))};
    $sc[7] = $metables[7]{&hashseq(substr($seq,11,3))};
    $sc[8] = $metables[8]{&hashseq(substr($seq,14,4))};
    my $finalscore = $sc[0] * $sc[1] * $sc[2] * $sc[3] * $sc[4] / ($sc[5] * $sc[6] * $sc[7] * $sc[8]);
    return $finalscore;
}    
    
sub getrest3{
    my $seq = shift;
    my $seq_noconsensus = substr($seq,0,18).substr($seq,20,3);
    return $seq_noconsensus;
}

sub scoreconsensus3{
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my %bgd; 
    $bgd{'A'} = 0.27; $bgd{'C'} = 0.23; $bgd{'G'} = 0.23; $bgd{'T'} = 0.27;  
    my %cons1;
    $cons1{'A'} = 0.9903; $cons1{'C'} = 0.0032; $cons1{'G'} = 0.0034; $cons1{'T'} = 0.0030;
    my %cons2;
    $cons2{'A'} = 0.0027; $cons2{'C'} = 0.0037; $cons2{'G'} = 0.9905; $cons2{'T'} = 0.0030;
    my $addscore = $cons1{$seqa[18]} * $cons2{$seqa[19]}/ ($bgd{$seqa[18]} * $bgd{$seqa[19]}); 
    return $addscore;
}
##### calculate 3p splice strength #####

sub sanity {
    my $seq = shift;
    my @seqa = split(//,uc($seq));
    my $Nr=scalar(@seqa);
    my $flag=0;
    for(my $i=0; $i<$Nr; $i++) {
        if ($seqa[$i] !~m/[ATCG]/) {$flag++;}
    }
    return $flag;
}

my %SMotif;
$SMotif{"GTAG"}=3;
$SMotif{"GCAG"}=2;
$SMotif{"ATAC"}=1;
$SMotif{"CTAC"}=-3;
$SMotif{"CTGC"}=-2;
$SMotif{"GTAT"}=-1;

my %reported;
my $f=0; 
foreach my $chr (sort keys %uniq) {
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

foreach my $id (keys %READ) {
	my @a=split("\t",$READ{$id});
    if ($debug eq -1) { print $READ{$id},"\n"; }
	# 0                       1       2       3       4       5			6       7    		8    		9     	10      11		12		13			14			15
    #						  in_len  out_len len     chr     seq_s		len     pos_s	  	pos_e	 	str     chr     seq_s   len     pos_s		pos_e		strand
    # newid-1__1	          75	  75	  75	  11	  0	        39	    118352428	118352466	-	    9	    37	    38	    20365667	20365704	+
    # 11      split   exon    118352430       118352807       protein_coding  +       MLL     ENSG00000118058___19___67
    # 9       split   exon    20365667        20365742        protein_coding  -       MLLT3   ENSG00000171843___16___30
    # 11      split   exon    118353137       118353210       protein_coding  +       MLL     ENSG00000118058___22___67
    # newid-2__1	          75	  75	  75	  9	      0	        41	    20365706	20365746	+	    11	    37	    38	    118353173	118353210	-
    # $left_seq  => 5' splice site
    # $right_seq => 3' splice site
    if ($a[9] eq "-") {
		my $overlap=$a[5]+$a[6]-$a[11];
		my $start=$a[7];
		my $end=0;
		my $chr1=$a[4];
		my $chr2=$a[10];
		my $left_seq="";
		my $left_max=-99999;
		my $left_pos=0;
		my $left_flag=0;
		my $right_seq="";
		my $right_max=-99999;
		my $right_pos=0;
		my $right_flag=0;
        my $SMscoren=-99999;
		# R-m+	[9876543AG|0123456789]		R-m+	[9876543210|GU3456789]		pos = 0 on exon
		if ($a[15] eq "-") {
			$end=$a[14];
			$left_seq=substr($genome{$chr2},($end-2*$Extend-1),(4*$Extend + 1));
			if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$right_seq,"\t",$left_seq,"\t",$overlap,"\n"; next; }
			if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr2,$end,"right_404"),"\t",join("\t",@a),"\n"; next;}
            $right_seq=substr($genome{$chr1},($start-2*$Extend-1),(4*$Extend + 1));
			if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr1,$start,"left_404"),"\t",join("\t",@a),"\n"; next;}
		}
		# R-m+	[9876543AG|0123456789]		R+m-	[9876543AC|0123456789]		pos = 0 on exon
		else {
			$end=$a[13];
			$left_seq=substr($genome{$chr2},($end-2*$Extend-1),(4*$Extend + 1));
			if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$right_seq,"\t",$left_seq,"\t",$overlap,"\n"; next; }
			if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr2,$end,"right_404"),"\t",join("\t",@a),"\n"; next;}
			my $tmp1=$left_seq;
			$tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
			my $tmp2=scalar reverse $tmp1;
			$left_seq=$tmp2;
            $right_seq=substr($genome{$chr1},($start-2*$Extend-1),(4*$Extend + 1));
			if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr1,$start,"left_404"),"\t",join("\t",@a),"\n"; next;}
		}
		{
			# sanity check for non-[ATCG] nt
			if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
			if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
			# check canonical splice-motif
			my $SMS=0;	# score of motif
			my $PMS=0;	# position of best motif
			if ($debug eq 1){print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n";}
			for(my $k=0; $k<=$overlap; $k++) {
				my $ml="";
				my $mr="";
				$ml=substr($left_seq,2*$Extend+1-$k,2);
				$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
				if (($ml ne "") and ($mr ne "")) {
					my $motif=$ml.$mr;
					if ($debug eq 1){print $k,"\t",$motif,"\n";}
                    if ((exists $SMotif{$motif}) and ($SMotif{$motif} > 0)) {
                        my $str1="";
                        $str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
                        my $str2="";
                        $str2=uc substr($right_seq,(2*$Extend-2+$overlap-$k-20+2),23);
                        if ((length($str1) eq 9) and (length($str2) eq 23)) {
                            my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                            my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                            my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                            if ($SMscoren < $sumt) { $SMscoren=$sumt; $SMS=$SMotif{$motif}; $PMS=$k; }
                        }
                    }
					#if ((exists $SMotif{$motif}) and ($SMotif{$motif} > 0)) { if (abs($SMotif{$motif}) > abs($SMS) ){ $SMS=$SMotif{$motif}; $PMS=$k; } }
				}
			}
			if ($debug eq 1){print "\n",$SMS,"\t",$PMS,"\n";}
			if ($SMS > 0) {
				# found junction on the plus strand
				my $str1="";
				$str1=uc substr($left_seq,(2*$Extend+1-$PMS-3),9);
				my $str2="";
				$str2=uc substr($right_seq,(2*$Extend-2+$overlap-$PMS-20+2),23);
				if ($debug eq 1) { print "str1 = $str1\nstr2 = $str2\n"; }
				if ((length($str1) eq 9) and (length($str2) eq 23)) {
					$left_max=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
					$right_max=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
					my $sum=sprintf("%.2f",($left_max + $right_max));
					$a[9]=$a[9]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."+5";
					$a[15]=$a[15]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."+3";
					$a[15]=$a[15]."\t".$sum."\t".$SMS."\t".$PMS;
					print OUT join("\t",@a),"\n";
					my @a=split("\t",$READ{$id});
					$a[11]=$a[11]+$PMS;
					$a[12]=$a[12]-$PMS;
					$a[6]=$a[6]-($overlap-$PMS);
					$a[7]=$a[7]+($overlap-$PMS);
                    if ($a[15] eq "-") { $a[14]=$a[14]-$PMS;}
                    else { $a[13]=$a[13]+$PMS; }
					$a[15]=$a[15]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."+"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
					print OUT1 join("\t",@a),"\n";
                    print OUT3 join("\t",@a),"\n";
					next;
				}
			}
			# 0                       1       2       3       4       5			6       7    		8    		9     	10      11		12		13			14			15
			#						  in_len  out_len len     chr     seq_s		len     pos_s	  	pos_e	 	str     chr     seq_s   len     pos_s		pos_e		strand
			# newid-32171/1__3        151     151     151     13      0         70      94111709    94111778    -       13      67      84      94125952    94126035    -       
			else {
                $SMscoren=-99999;
                $left_max=0;
                $left_pos=0;
                $right_max=0;
                $right_pos=0;
                for(my $k=0; $k<=$overlap; $k++) {
                    my $str1="";
                    $str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
                    my $str2="";
                    $str2=uc substr($right_seq,(2*$Extend-2+$overlap-$k-20+2),23);
                    if ((length($str1) eq 9) and (length($str2) eq 23)) {
                        my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                        my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                        my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                        if ($SMscoren < $sumt) {
                            $SMscoren=$sumt;
                            $PMS=$k;
                            $left_pos=$k;
                            $left_max=$left_maxt;
                            $right_pos=$overlap-$k;
                            $right_max=$right_maxt;
                            $left_flag=+5;
                            $right_flag=+3;
                        }
                    }    
                }
                if ($debug eq 1){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
                $a[9]=$a[9]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."+5";
                $a[15]=$a[15]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."+3";
                $a[15]=$a[15]."\t".$SMscoren."\t".$SMS."\t".$PMS;
                print OUT join("\t",@a),"\n";
                my @a=split("\t",$READ{$id});
                $a[11]=$a[11]+$PMS;
                $a[12]=$a[12]-$PMS;
                $a[6]=$a[6]-($overlap-$PMS);
                $a[7]=$a[7]+($overlap-$PMS);
                if ($a[15] eq "-") { $a[14]=$a[14]-$PMS;}
                else { $a[13]=$a[13]+$PMS; }
                $a[15]=$a[15]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."+"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
                print OUT2 join("\t",@a),"\n";
                print OUT3 join("\t",@a),"\n";
			}
		}
	}
	else {
		my $overlap=$a[5]+$a[6]-$a[11];
		my $start=$a[8];
		my $end=0;
		my $chr1=$a[4];
		my $chr2=$a[10];
		my $left_seq="";
		my $left_max=-99999;
		my $left_pos=0;
		my $left_flag=0;
		my $right_seq="";
		my $right_max=-99999;
		my $right_pos=0;
		my $right_flag=0;
		# R+m-	[9876543210|CT2345678]		R+m-	[9876543AC|0123456789]		pos = 0 on exon
		if ($a[15] eq "+") {
			$end=$a[13];
            $left_seq=substr($genome{$chr2},($end-2*$Extend-1),(4*$Extend + 1));
            if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr1,$start,"left_404"),"\t",join("\t",@a),"\n"; next;}
            $right_seq=substr($genome{$chr1},($start-2*$Extend-1),(4*$Extend + 1));
            if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr2,$end,"right_404"),"\t",join("\t",@a),"\n"; next;}
			if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
            my $tmp1=$left_seq;
			$tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
			my $tmp2=scalar reverse $tmp1;
			$left_seq=$tmp2;
			$tmp1=$right_seq;
			$tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
			$tmp2=scalar reverse $tmp1;
			$right_seq=$tmp2;
		}
		# R+m-	[9876543210|CT2345678]		R-m+	[9876543210|GU3456789]		pos = 0 on exon
		else {
			$end=$a[14];
            $left_seq=substr($genome{$chr2},($end-2*$Extend-1),(4*$Extend + 1));
            if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr1,$start,"left_404"),"\t",join("\t",@a),"\n"; next;}
            $right_seq=substr($genome{$chr1},($start-2*$Extend-1),(4*$Extend + 1));
            if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr2,$end,"right_404"),"\t",join("\t",@a),"\n"; next;}
			if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
			my $tmp1=$right_seq;
			$tmp1=~tr/[atcgATCG]/[TAGCTAGC]/;
			my $tmp2=scalar reverse $tmp1;
			$right_seq=$tmp2;
		}
		{
			# sanity check for non-[ATCG] nt
			if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
			if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
			# check canonical splice-motif
			my $SMS=0;	# score of motif
			my $PMS=0;	# position of best motif
            my $SMscoren=-99999;
			if ($debug eq 1){print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n";}
			for(my $k=0; $k<=$overlap; $k++) {
				my $ml="";
				my $mr="";
				$ml=substr($left_seq,2*$Extend+1-$k,2);
				$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
				if (($ml ne "") and ($mr ne "")) {
					my $motif=$ml.$mr;
					if ($debug eq 1){print $k,"\t",$motif,"\n";}
                    if ((exists $SMotif{$motif}) and ($SMotif{$motif} > 0)) {
                        my $str1="";
                        $str1=uc substr($left_seq,(2*$Extend+1-$k-3),9);
                        my $str2="";
                        $str2=uc substr($right_seq,(2*$Extend-2+$overlap-$k-20+2),23);
                        if ((length($str1) eq 9) and (length($str2) eq 23)) {
                            my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                            my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                            my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                            if ($SMscoren < $sumt) { $SMscoren=$sumt; $SMS=$SMotif{$motif}; $PMS=$k; }
                        }
                    }
                    #if ((exists $SMotif{$motif}) and ($SMotif{$motif} < 0)) { if (abs($SMotif{$motif}) > abs($SMS) ){ $SMS=$SMotif{$motif}; $PMS=$k; } }
				}
			}
			if ($debug eq 1){print "\n",$SMS,"\t",$PMS,"\n";}
			if ($SMS > 0) {
                my $str1="";
				$str1=uc substr($left_seq,(2*$Extend+1-$PMS-3),9);
				my $str2="";
				$str2=uc substr($right_seq,(2*$Extend-2+$overlap-$PMS-20+2),23);
				if ($debug eq 1) { print "str1 = $str1\nstr2 = $str2\n"; }
				if ((length($str1) eq 9) and (length($str2) eq 23)) {
                    $left_max=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                    $right_max=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
					my $sum=sprintf("%.2f",($left_max + $right_max));
					$a[9]=$a[9]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."-5";
					$a[15]=$a[15]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."-3";
					$a[15]=$a[15]."\t".$sum."\t".$SMS."\t".$PMS;
					print OUT join("\t",@a),"\n";
					my @a=split("\t",$READ{$id});
					$a[6]=$a[6]-($overlap-$PMS);
					$a[8]=$a[8]-($overlap-$PMS);
					$a[11]=$a[11]+$PMS;
					$a[12]=$a[12]-$PMS;
					if ($a[15] eq "-") { $a[14]=$a[14]-$PMS; }
                    else { $a[13]=$a[13]+$PMS; }
					$a[15]=$a[15]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."-"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
					print OUT1 join("\t",@a),"\n";
                    print OUT3 join("\t",@a),"\n";
					next;
				}
			}
			# 0                       1       2       3       4       5			6       7    		8    		9     	10      11		12		13			14			15
			#						  in_len  out_len len     chr     seq_s		len     pos_s	  	pos_e	 	str     chr     seq_s   len     pos_s		pos_e		strand
			# newid-16829/1__5        151     151     151     1       0         75      157809618   157809692   +       1       71      80      157800282   157800361   +
			# newid-16829/1__5        151     151     151     1       0         73      157809618   157809690   +       1       73      78      157800284   157800361   +       -9411   14.50 6.04     8.46    -       4       2       2       -3      2
			else {
                $SMscoren=-99999;
                $left_max=0;
                $left_pos=0;
                $right_max=0;
                $right_pos=0;
                for(my $k=0; $k<=$overlap; $k++) {
                    my $str1="";
                    $str1=uc substr($left_seq,(2*$Extend+1-$PMS-3),9);
                    my $str2="";
                    $str2=uc substr($right_seq,(2*$Extend-2+$overlap-$PMS-20+2),23);
                    if ((length($str1) eq 9) and (length($str2) eq 23)) {
                        my $left_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
                        my $right_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
                        my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                        if ($SMscoren < $sumt) {
                            $SMscoren=$sumt;
                            $PMS=$k;
                            $left_pos=$k;
                            $left_max=$left_maxt;
                            $right_pos=$overlap-$k;
                            $right_max=$right_maxt;
                            $left_flag=-5;
                            $right_flag=-3;
                        }
                    }    
                }
                if ($debug eq 1){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
                $a[9]=$a[9]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."-5";
                $a[15]=$a[15]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."-3";
                $a[15]=$a[15]."\t".$SMscoren."\t".$SMS."\t".$PMS;
                print OUT join("\t",@a),"\n";
                my @a=split("\t",$READ{$id});
                $a[6]=$a[6]-($overlap-$PMS);
                $a[8]=$a[8]-($overlap-$PMS);
                $a[11]=$a[11]+$PMS;
                $a[12]=$a[12]-$PMS;
                if ($a[15] eq "-") { $a[14]=$a[14]-$PMS; }
                else { $a[13]=$a[13]+$PMS; }  
                $a[15]=$a[15]."\t".$SMscoren."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."-"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
                print OUT2 join("\t",@a),"\n";
                print OUT3 join("\t",@a),"\n";
			}
		}
	}
}
close IN;
close OUT;
close OUT1;
close OUT2;
close OUT3;
close OUTerr;
close OUTp;







