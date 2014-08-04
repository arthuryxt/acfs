#!/usr/bin/perl -w
use strict;
# check splicing sites for every single circRNA candidates. Require stranded sequencing and the sequencing reads are reverse-complementary to mRNA
die "Usage: $0  \"CB_splice_folder\"  \"selected_circ_reads\"  \"genome_location\"  \"output basename\"  \"\(optional\) extend N bases\"  " if (@ARGV < 4);
my $DIR=$ARGV[0];		# /home/arthur/CB_splice/
my $anno=$ARGV[1];      # unmap.parsed.2pp
my $genome=$ARGV[2];    # /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes/
my $fileout=$ARGV[3];
my $Extend=15;           # 15nt by default. as 3' splice strength need 23nt, 20nt from intron and 3 from exon.
if (scalar(@ARGV) > 4) {$Extend=$ARGV[3];}
my %uniq;

my $command="rm -f Step2_finished";
system($command);

open IN,$anno;
open OUT,">".$fileout;
open OUT1,">".$fileout.".0";	# splice-motif found				trust hit  CB prediction
open OUT2,">".$fileout.".1";	# closure == 0,						trust best CB prediction
open OUT3,">".$fileout.".2";	# closure != 0 but overlap == 0,	report original_positioned CB prediction 
open OUT4,">".$fileout.".3";	# closure != 0 and overlap != 0,	report original_positioned CB prediction
open OUTerr,">".$fileout."_withN";
open OUTp,">".$fileout.".pgap";
while(<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[13] < 0) {
        if ($a[2]=~m/^chromosome/i) {$a[2]=~s/chromosome//i;}
        if ($a[2]=~m/^chr/i) {$a[2]=~s/chr//i;}
		$a[6]=$a[6]-1;
		$a[11]=$a[11]-1;
        $uniq{$a[2]}{$a[6]}{$a[10]}{$a[0]}=join("\t",@a);
    }
    else {print OUTp join("\t",@a),"\n";}
}
close IN;

my $tmpfile1=$DIR."/me2x5";
#my %me2x5 = &makescorematrix("~/CB_splice/me2x5");
my %me2x5 = &makescorematrix($tmpfile1);
my $tmpfile2=$DIR."/splicemodels/splice5sequences";
#my %seq = &makesequencematrix("/home/arthur/CB_splice/splicemodels/splice5sequences");
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
sub hashseq{
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
    #my $dir = "~/CB_splice/splicemodels/";
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
	$num++;
    }
    return @metables;
}
sub makewmmscores{
	my $dir = $DIR."/splicemodels/";
    #my $dir = "/home/arthur/CB_splice/splicemodels/";
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
	$num++;
    }
    return @metables;
}
sub makemmscores{
	my $dir = $DIR."/splicemodels/";
    #my $dir = "/home/arthur/CB_splice/splicemodels/";
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
    close IN1;
    foreach my $start (sort keys %{$uniq{$chr}}) {
        foreach my $end (sort keys %{$uniq{$chr}{$start}}) {
            foreach my $seqid (keys %{$uniq{$chr}{$start}{$end}}) {
                my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
				if ($a[7] eq "+") {
					# R+m-
					my $overlap=$a[3]+$a[4]-$a[8];
					my $left_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
					if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,"left_404"),"\n"; next;}
					my $left_max=-99999;
					my $left_pos=0;
					my $left_flag=0;
					my $right_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
					if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
					if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,$end,$seqid,"right_404"),"\n";next;}
					my $right_max=-99999;
					my $right_pos=0;
					my $right_flag=0;
					if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
					if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
					my $SMS=0;	
					my $PMS=0;	
					for(my $k=0; $k<=$overlap; $k++) {
						my $ml="";
						my $mr="";
						$ml=substr($left_seq,2*$Extend+1-$k,2);
						$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
						if (($ml ne "") and ($mr ne "")) {
							my $motif=$ml.$mr;
							if ((exists $SMotif{$motif}) and ($SMotif{$motif} < 0)) { if (abs($SMotif{$motif}) > abs($SMS) ){ $SMS=$SMotif{$motif}; $PMS=$k; } }
						}
					}
					if ($SMS < 0) {
						{
						my $t=$left_seq;					
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_left_seq=scalar reverse $t;		
						$t=$right_seq;					
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_right_seq=scalar reverse $t;		
						my $str1="";
						$str1=uc substr($rc_left_seq,(2*$Extend+$PMS-20),23);
						my $str2="";
						$str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$PMS),9);
						if ((length($str1) eq 23) and (length($str2) eq 9)) {
							$left_max=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
							$right_max=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
							my $sum=sprintf("%.2f",($left_max + $right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."-5";
							$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."-3";
							$a[13]=$a[13]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							$a[4]=$a[4]-$PMS;
							$a[6]=$a[6]-$PMS;
							$a[8]=$a[8]+($overlap-$PMS);
							$a[9]=$a[9]-($overlap-$PMS);
							$a[10]=$a[10]+($overlap-$PMS);
							$a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."-"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
							print OUT1 join("\t",@a),"\n";
							next;
						}
						}
					}
					else {
						#
						my $t=$left_seq;
					    $t=~tr/[atcgATCG]/[TAGCTAGC]/;
					    my $rc_left_seq=scalar reverse $t;		
					    $t=$right_seq;
					    $t=~tr/[atcgATCG]/[TAGCTAGC]/;
					    my $rc_right_seq=scalar reverse $t;		
					    my $rc_left_max=-99999;
					    my $rc_left_pos=0;
					    my $rc_left_flag=0;
					    my $rc_right_max=-99999;
					    my $rc_right_pos=0;
					    my $rc_right_flag=0;
					    for (my $i=0; $i<((4*$Extend +1) - 23); $i++) {
					        my $str=substr($rc_left_seq,$i,23);
					        if (length($str) ne 23) {print OUTerr join("\t",@a),"\n";last;}
					        my $tmp_score=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
					        if ($tmp_score > $rc_left_max) {$rc_left_max = $tmp_score; $rc_left_pos=$i; $rc_left_flag=-5;}   
					    }
					    for (my $i=0; $i<((2*$Extend +1) - 9); $i++) {
					        my $str=uc substr($rc_right_seq,($Extend + $i),9);
					        if (length($str) ne 9) {print OUTerr join("\t",@a),"\n";last;}
					        my $tmp_score=sprintf("%.2f",&log2(&scoreconsensus5($str)*$me2x5{$seq{&getrest5($str)}}));
					        if ($tmp_score > $rc_right_max) {$rc_right_max = $tmp_score; $rc_right_pos=$i; $rc_right_flag=-3;}   
					    }
						{
							my $sum=sprintf("%.2f",($rc_left_max + $rc_right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$rc_left_max."\t".$rc_left_pos."\t".$rc_left_flag;
							$a[12]=$a[12]."\t".$right_seq."\t".$rc_right_max."\t".$rc_right_pos."\t".$rc_right_flag;
							$a[13]=$a[13]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							my $strand="+";
							if ($rc_left_flag < 0) {$strand="-";}
							my $left_move=2*$Extend - ($rc_left_pos + 20);
							my $right_move=($Extend + $rc_right_pos + 3) - (2*$Extend + 1);
							my $closure=0;
							if ($a[7] eq "+") { $closure = $left_move + $right_move - $overlap; }
							else { $closure = $left_move + $right_move - $overlap;}
							if ($closure eq 0) {
							    $a[4]=$a[4]-$left_move;
							    $a[6]=$a[6]-$left_move;
								$a[8]=$a[8]+$right_move;
							    $a[9]=$a[9]-$right_move;
							    $a[10]=$a[10]+$right_move;
							    $a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$rc_left_max)."\t".sprintf("%.2f",$rc_right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							    print OUT2 join("\t",@a),"\n";
							}
							else {
								if ($overlap eq 0) {
							        my $str=uc substr($rc_left_seq,(2*$Extend - 20),23);
							        $rc_right_max=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
							        my $str2=uc substr($rc_right_seq,(2*$Extend + 1 - 3),9);
							        $rc_left_max=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
							        $sum=sprintf("%.2f",($rc_left_max + $rc_right_max));
							        $a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$rc_left_max)."\t".sprintf("%.2f",$rc_right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							        print OUT3 join("\t",@a),"\n";
							    }
							    else {
							        my $str=uc substr($rc_left_seq,(2*$Extend - 20),23);
							        $rc_right_max=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
							        my $str2=uc substr($rc_right_seq,(2*$Extend + 1 - 3),9);
							        $rc_left_max=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
							        $sum=sprintf("%.2f",($rc_left_max + $rc_right_max));
									$a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$rc_left_max)."\t".sprintf("%.2f",$rc_right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							        print OUT4 join("\t",@a),"\n";
							    }
							}
						}
					}
				}
				else {
					# R-m+
					my $overlap=$a[8]+$a[9]-$a[3];
					my $left_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
					if (length($left_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,"left_404"),"\n"; next;}
					my $left_max=-99999;
					my $left_pos=0;
					my $left_flag=0;
					my $right_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
					if ($overlap > $Extend) { print OUTerr join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
					if (length($right_seq) ne (4*$Extend + 1)) {print OUTerr join("\t",$chr,$start,$end,$seqid,"right_404"),"\n";next;}
					my $right_max=-99999;
					my $right_pos=0;
					my $right_flag=0;
					if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
					if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
					my $SMS=0;	
					my $PMS=0;	
					for(my $k=0; $k<=$overlap; $k++) {
						my $ml="";
						my $mr="";
						$ml=substr($left_seq,2*$Extend+1-$k,2);
						$mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
						if (($ml ne "") and ($mr ne "")) {
							my $motif=$ml.$mr;
							if ((exists $SMotif{$motif}) and ($SMotif{$motif} > 0)) { if (abs($SMotif{$motif}) > abs($SMS) ){ $SMS=$SMotif{$motif}; $PMS=$k; } }
						}
					}
					if ($SMS > 0) {
						my $str1="";
						$str1=uc substr($left_seq,(2*$Extend+1-$PMS-3),9);
						my $str2="";
						$str2=uc substr($right_seq,(2*$Extend-2+$overlap-$PMS-20+2),23);
						if ((length($str1) eq 9) and (length($str2) eq 23)) {
							$left_max=sprintf("%.2f",&log2(&scoreconsensus5($str1)*$me2x5{$seq{&getrest5($str1)}}));
							$right_max=sprintf("%.2f", &log2(&scoreconsensus3($str2)*&maxentscore(&getrest3($str2),\@metables)));
							my $sum=sprintf("%.2f",($left_max + $right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$PMS."\t"."+5";
							$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".($overlap-$PMS)."\t"."+3";
							$a[13]=$a[13]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							$a[3]=$a[3]+$PMS;
							$a[4]=$a[4]-$PMS;
							$a[6]=$a[6]-$PMS;
							$a[9]=$a[9]-($overlap-$PMS);
							$a[10]=$a[10]+($overlap-$PMS);
							$a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t"."+"."\t".$overlap."\t".$PMS."\t".($overlap-$PMS)."\t".$SMS."\t".$PMS;
							print OUT1 join("\t",@a),"\n";
							next;
						}
					}
					else {
						#
						for (my $i=0; $i<((2*$Extend +1) - 9); $i++) {
							my $str=uc substr($left_seq,($Extend + $i),9);
							if (length($str) ne 9) {print OUTerr join("\t",@a),"\n"; last;}
							my $tmp_score=sprintf("%.2f",&log2(&scoreconsensus5($str)*$me2x5{$seq{&getrest5($str)}}));
							if ($tmp_score > $left_max) {$left_max = $tmp_score; $left_pos=$i; $left_flag=+5;}   
						}
						for (my $i=0; $i<((4*$Extend +1) - 23); $i++) {
						    my $str=uc substr($right_seq,$i,23);
						    if (length($str) ne 23) {print OUTerr join("\t",@a),"\n";last;}
						    my $tmp_score=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
						    if ($tmp_score > $right_max) {$right_max = $tmp_score; $right_pos=$i; $right_flag=+3;}   
						}
						{
							my $sum=sprintf("%.2f",($left_max + $right_max));
							$a[7]=$a[7]."\t".$left_seq."\t".$left_max."\t".$left_pos."\t".$left_flag;
							$a[12]=$a[12]."\t".$right_seq."\t".$right_max."\t".$right_pos."\t".$right_flag;
							$a[13]=$a[13]."\t".$sum."\t".$SMS."\t".$PMS;
							print OUT join("\t",@a),"\n";
							my @a=split("\t",$uniq{$chr}{$start}{$end}{$seqid});
							my $strand="+";
							if ($left_flag < 0) {$strand="-";}
							my $left_move=(2*$Extend + 1) - ($Extend + $left_pos + 3);
							my $right_move=($right_pos + 20 + 1)  - (2*$Extend - 1);
							my $closure=0;
							if ($a[7] eq "+") { $closure = $left_move + $right_move - $overlap; }
							else { $closure = $left_move + $right_move - $overlap; }
							if ($closure eq 0) {
								$a[3]=$a[3]+$left_move;
							    $a[4]=$a[4]-$left_move;
							    $a[6]=$a[6]-$left_move;
							    $a[9]=$a[9]+$right_move;
							    $a[10]=$a[10]+$right_move;
							    $a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							    print OUT2 join("\t",@a),"\n";
							}
							else {
								if ($overlap eq 0) {
							        my $str=uc substr($left_seq,(2*$Extend + 1 - 3),9);
							        $left_max=sprintf("%.2f",&log2(&scoreconsensus5($str)*$me2x5{$seq{&getrest5($str)}}));
							        $str=uc substr($right_seq,(2*$Extend -20),23);
							        $right_max=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
							        $sum=sprintf("%.2f",($left_max + $right_max));
							        $a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							        print OUT3 join("\t",@a),"\n";
							    }
							    else {
									my $str=uc substr($left_seq,(2*$Extend + 1 - 3),9);
							        $left_max=sprintf("%.2f",&log2(&scoreconsensus5($str)*$me2x5{$seq{&getrest5($str)}}));
							        $str=uc substr($right_seq,(2*$Extend -20),23);
							        $right_max=sprintf("%.2f", &log2(&scoreconsensus3($str)*&maxentscore(&getrest3($str),\@metables)));
							        $sum=sprintf("%.2f",($left_max + $right_max));
							        $a[13]=$a[13]."\t".$sum."\t".sprintf("%.2f",$left_max)."\t".sprintf("%.2f",$right_max)."\t".$strand."\t".$overlap."\t".$left_move."\t".$right_move."\t".$SMS."\t".$PMS;
							        print OUT4 join("\t",@a),"\n";
							    }
							}
						}
					}
				}

			}
	    }
    }
}
if ($f eq 1) {print "\t=======\tyes\n";}
elsif ($f eq -1) {print "\t=======\tno\n";}

my $file0=$fileout.".0";
my $file1=$fileout.".1";
my $file2=$fileout.".sum";
$command="cat $file0 $file1 | sort -k3,3 -k7,7n -k11,11n > $file2";
system($command);
#cat human_unmap.parsed.2pp.S2.0 human_unmap.parsed.2pp.S2.1 | sort -k3,3 -k7,7n -k11,11n > human_unmap.parsed.2pp.S2.sum01

open(OUTFLAG,">Step2_finished");
print OUTFLAG "Step2_finished\n";
close OUTFLAG;



