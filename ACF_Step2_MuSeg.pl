#!/usr/bin/perl -w
use strict;
# rescue circ spanning very short exons
die "Usage: $0  \"CB_splice_folder\"   \"multiple_segments\"  \"split_exon_gtf\"   \"genome_folder\"   \"output\"  \"\(optional\) extend N bases\"  \"\(DEBUG\)\"" if (@ARGV < 5);
my $DIR=$ARGV[0];		# /home/arthur/CB_splice/
my $filein1=$ARGV[1];
my $filein2=$ARGV[2];
my $genome=$ARGV[3];
my $fileout=$ARGV[4];
my $Ext=4;
if (scalar(@ARGV) > 5) {$Ext=$ARGV[5];}
die "Usage: $0  \"multiple_segments\"  \"split_exon_gtf\"   \"output\"  \"\(optional\) extend positive N bases\"  \"\(DEBUG set to 1\)\"" if ($Ext < 0);
print "Extension = $Ext\n";
my $debug=-1;
if (scalar(@ARGV) > 6) {$debug=$ARGV[6];}
my $cutoff=0.9;
my $Extend=15;
my $closegap=50;
my $command="rm -f Step2_MuSeg_finished";
system($command);

my %uniq;
if ($filein2 ne "no") {
    open IN1,$filein2;
    while(<IN1>) {
        chomp;
        my @a=split("\t",$_);
        if ($a[0]=~m/^chromosome/i) {$a[0]=~s/chromosome//i;}
        if ($a[0]=~m/^chr/i) {$a[0]=~s/chr//i;}
        my $id=$a[7]."___".$a[8];
        if ($a[6] eq "+") {
        # 5' of exon
        my $base=int($a[3]/1000);
        if (exists $uniq{$a[0]}{$base}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base}{$a[3]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base}{$a[3]}=join("___",10,$id); }
        
        if (exists $uniq{$a[0]}{$base+1}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base+1}{$a[3]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base+1}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base+1}{$a[3]}=join("___",10,$id); }
        
        if (exists $uniq{$a[0]}{$base-1}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base-1}{$a[3]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base-1}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base-1}{$a[3]}=join("___",10,$id); }
        # 3' of exon
        $base=int($a[4]/1000);
        if (exists $uniq{$a[0]}{$base}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base}{$a[4]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base}{$a[4]}=join("___",1,$id); }
        
        if (exists $uniq{$a[0]}{$base+1}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base+1}{$a[4]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base+1}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base+1}{$a[4]}=join("___",1,$id); }
        
        if (exists $uniq{$a[0]}{$base-1}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base-1}{$a[4]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base-1}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base-1}{$a[4]}=join("___",1,$id); }
        }
        else {
        # 5' of exon
        my $base=int($a[4] / 1000);
        if (exists $uniq{$a[0]}{$base}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base}{$a[4]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base}{$a[4]}=join("___",10,$id); }
        
        if (exists $uniq{$a[0]}{$base+1}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base+1}{$a[4]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base+1}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base+1}{$a[4]}=join("___",10,$id); }
        
        if (exists $uniq{$a[0]}{$base-1}{$a[4]}) {
            my @b=split("___",$uniq{$a[0]}{$base-1}{$a[4]});
            if ($b[0] < 10) { $b[0]+=10; $uniq{$a[0]}{$base-1}{$a[4]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base-1}{$a[4]}=join("___",10,$id); }
        # 3' of exon
        $base=int($a[3] / 1000);
        if (exists $uniq{$a[0]}{$base}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base}{$a[3]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base}{$a[3]}=join("___",1,$id); }
        
        if (exists $uniq{$a[0]}{$base+1}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base+1}{$a[3]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base+1}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base+1}{$a[3]}=join("___",1,$id); }
        
        if (exists $uniq{$a[0]}{$base-1}{$a[3]}) {
            my @b=split("___",$uniq{$a[0]}{$base-1}{$a[3]});
            if (($b[0] < 1) or ($b[0] eq 10)) { $b[0]+=1; $uniq{$a[0]}{$base-1}{$a[3]}=join("___",@b)}
        }
        else { $uniq{$a[0]}{$base-1}{$a[3]}=join("___",1,$id); }
        }
    }
    close IN1;
}

open IN2,$filein1;
open OUT1,">".$fileout.".S1.pp";
open OUT11,">".$fileout.".S1.ppparsed";
open OUT2,">".$fileout.".S1.pm";
open OUT22,">".$fileout.".S1.pmparsed";
my %READ;
my %Overlap;
while(<IN2>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[4]=~m/^chromosome/i) {$a[4]=~s/chromosome//i;}
    if ($a[4]=~m/^chr/i) {$a[4]=~s/chr//i;}
	$READ{$a[0]}=join("\t",$_);
    my $NR=scalar(@a)/5;
    my %strand;
	my $gap=-1;
    for(my $i=1; $i<$NR; $i++) {
        $a[5*$i+3]--;
		$strand{$a[5*$i + 4]}++;
		if ($debug eq 1) { print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4]),"\n"; }
		my $base=int(($a[5*$i + 2])/1000);
		if ($debug eq 1) { print $base,"\n"; }
		my $info="NA";
		if (exists $uniq{$a[4]}) {
		    if (exists $uniq{$a[4]}{$base}) {
				foreach my $pos (keys %{$uniq{$a[4]}{$base}}) {
				    if ($debug eq 1) {print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4],$pos,$uniq{$a[4]}{$base}{$pos}),"\n";}
				    if (abs($pos - $a[5*$i + 2]) <= $Ext) { $info=$uniq{$a[4]}{$base}{$pos}; last; }
				}
			}
		}
		$a[5*$i+2]=$a[5*$i + 2]."\t".$info;
	
		$base=int(($a[5*$i + 3])/1000);
		if ($debug eq 1) { print $base,"\n"; }
		$info="NA";
		if (exists $uniq{$a[4]}) {
		    if (exists $uniq{$a[4]}{$base}) {
				foreach my $pos (keys %{$uniq{$a[4]}{$base}}) {
				    if ($debug eq 1) {print join("\t",$a[0],$a[4],$a[5*$i],$a[5*$i+1],$a[5*$i+2],$a[5*$i+3],$a[5*$i+4],$pos,$uniq{$a[4]}{$base}{$pos}),"\n";}
				    if (abs($pos - $a[5*$i + 3]) <= $Ext) { $info=$uniq{$a[4]}{$base}{$pos}; last; }
				}
		    }
		}
		$a[5*$i+3]=$a[5*$i + 3]."\t".$info;
		#if ($i > 1) {
		#	my $tmp_gap=$a[5*$i-5] + $a[5*$i-5+1] - $a[5*$i];
		#	if ($tmp_gap > $gap) { $gap=$tmp_gap;}
		#}
		#$Overlap{$a[0]}=$gap;
    }
    my $ss=0;
    foreach my $i (keys %strand) { $ss++; }
    if ($ss eq 1) {
		print OUT1 join("\t",@a),"\n";
		
		my @Outer_left=qw(999999999 999999999 NA);
		my @Outer_right=qw(-1 -1 NA);	# where the Back Splice took place
		my @Inner_left=qw(-1 -1 NA);
		my @Inner_right=qw(-1 -1 NA);	# should be the first and the last segments, but NO (left < right) relation assumed!!!
		my %POS;
        my $leftmostrank=-1;
        my $rightmostrank=-1;
		for(my $i=1; $i<$NR; $i++) {
		    my @tmpa=split("\t",$a[5*$i+2]);
		    my @tmpb=split("\t",$a[5*$i+3]);
		    $POS{$tmpa[0]}{$i}=join("____",$tmpa[1],$tmpb[0],$tmpb[1],$a[5*$i+4],$a[5*$i],$a[5*$i+1]);
		    if ($Inner_left[0] eq -1) { $Inner_left[0]=$tmpa[0]; $Inner_left[1]=$tmpb[0]; $Inner_left[2]=$tmpa[1];}
		    $Inner_right[0]=$tmpa[0]; $Inner_right[1]=$tmpb[0]; $Inner_right[2]=$tmpb[1];
		    if ($Outer_left[0] > $tmpa[0]) {$Outer_left[0] = $tmpa[0]; $Outer_left[1]=$tmpb[0]; $Outer_left[2]=$tmpa[1]; $leftmostrank=$i;}
		    if ($Outer_right[1] < $tmpb[0]) {$Outer_right[0] = $tmpa[0]; $Outer_right[1]=$tmpb[0]; $Outer_right[2]=$tmpb[1]; $rightmostrank=$i;}
            if ($debug eq 20) {
                print join("\t",@Outer_left),"\t\t",join("\t",@Outer_right),"\t\t",join("\t",@Inner_left),"\t\t",join("\t",@Inner_right),"\n";
            }
            
		}
	
		# check if Inner border is truely inner
		my $flag=0;	# 0:circle;	1:polymerase backslide or templateswitch;	10-or-100: backslide	110:linear
		#if ( (($Inner_right[0] <= $Inner_left[1]) and ($Inner_left[1] <= $Inner_right[1]))  or  (($Inner_left[0] <= $Inner_right[1]) and ($Inner_right[1] <= $Inner_left[1])) ) {$flag++;}
		#for(my $i=2; ($i<($NR-1))and($flag eq 0); $i++) {
		#    my @tmpa=split("\t",$a[5*$i+2]);
		#    my @tmpb=split("\t",$a[5*$i+3]);
		#    if ( $Inner_left[1]  <  $Inner_right[0]) {
		#		if (($Inner_left[1] <= $tmpa[0]) and ($tmpa[0] <= $Inner_right[0])) {$flag+=10;}
		#		if (($Inner_left[1] <= $tmpb[0]) and ($tmpb[0] <= $Inner_right[0])) {$flag+=100;}
		#    }
		#    else {
		#		if (($Inner_right[1] <= $tmpa[0]) and ($tmpa[0] <= $Inner_left[0])) {$flag+=10;}
		#		if (($Inner_right[1] <= $tmpb[0]) and ($tmpb[0] <= $Inner_left[0])) {$flag+=100;}
		#    }
		#}
		# check for circularity
		#if ($flag eq 0) {
		#	my $last=-1;
        #    my $breakpoint=0;
		#	foreach my $left (sort {$a <=> $b} keys %POS) {
		#		foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
		#			if ($last eq -1) {}
		#			elsif (($last - $rank) > 1) {
		#				if ($a[-1] eq "+") {$breakpoint++;}
		#				else {$flag++;}
		#			}
		#			elsif (($rank - $last) > 1) {
		#				if ($a[-1] eq "-") {$breakpoint++;}
		#				else {$flag++;}
		#			}
		#			$last=$rank;
		#		}
		#	}
            #if($breakpoint ne 1){$flag=1000+$breakpoint;}
		#}
		# report the segmental order with sticky-end info
		my $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$a[-1],$flag,$Outer_left[0],$Outer_right[1]);
		if ($debug eq 3) { $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$a[-1],$flag,$Outer_left[0],$Outer_right[1],$Inner_left[0],$Inner_left[1],$Inner_right[0],$Inner_right[1]) }
		foreach my $left (sort {$a <=> $b} keys %POS) {
		    foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
                if (($rank eq $rightmostrank) or ($rank eq $leftmostrank)) {
                    $info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
                }
				#$info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
		    }
		}
		# report
		print OUT11 $info,"\n";
    }
    else {
		print OUT2 join("\t",@a),"\n";
		
		my $Plus_left=999999999;
		my $Plus_right=-1;
		my $Minus_left=999999999;
		my $Minus_right=-1;
		my %POS;
		for(my $i=1; $i<$NR; $i++) {
		    my @tmpa=split("\t",$a[5*$i+2]);
		    my @tmpb=split("\t",$a[5*$i+3]);
		    $POS{$tmpa[0]}{$i}=join("____",$tmpa[1],$tmpb[0],$tmpb[1],$a[5*$i+4]);
		    if ($a[5*$i+4] eq "+") {
				if ($Plus_left >= $tmpa[0]) {$Plus_left = $tmpa[0];}
				if ($Plus_right <= $tmpb[0]) {$Plus_right = $tmpb[0];}
		    }
		    else {
				if ($Minus_left >= $tmpa[0]) {$Minus_left = $tmpa[0];}
				if ($Minus_right <= $tmpb[0]) {$Minus_right = $tmpb[0];}
		    }
		}
    	# Not necessary to check if Inner border is truely inner
		my $flag=1;	# 0:circle;	1:polymerase backslide or templateswitch;	10-or-100: backslide	110:linear
		# report the segmental order with sticky-end info
		my $info=join("\t",$a[0],$a[1],$a[2],$a[3],$a[4],$Plus_left,$Plus_right,$Minus_left,$Minus_right);
		foreach my $left (sort {$a <=> $b} keys %POS) {
		    foreach my $rank (sort{$a <=> $b} keys %{$POS{$left}}) {
				$info=$info."\t".$rank."____".$left."____".$POS{$left}{$rank};
		    }
		}
		# report
		print OUT22 $info,"\n";
    }
}
close IN2;
close OUT1;		# ".S1.pp"
close OUT11;	# ".S1.ppparsed"
close OUT2;		# ".S1.pm"
close OUT22;	# ".S1.pmparsed"


my %SMotif;
$SMotif{"GTAG"}=3;
$SMotif{"GCAG"}=2;
$SMotif{"ATAC"}=1;
$SMotif{"CTAC"}=-3;
$SMotif{"CTGC"}=-2;
$SMotif{"GTAT"}=-1;

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

open INf, $fileout.".S1.ppparsed";
open OUTf, ">".$fileout.".S1.ppparsed.id";
open OUTfm, ">".$fileout.".S1.ppparsed.loop";
open OUTf1, ">".$fileout.".S1";
open OUTfa1, ">".$fileout.".S1.fa";
open OUTf2, ">".$fileout.".S2";
open OUTfa2, ">".$fileout.".S2.fa";
open OUTf12, ">".$fileout.".S1.bed12";
open OUTf22, ">".$fileout.".S2.bed12";
open OUTerr,">".$fileout."_withN";

my %uniqf;
my %uniqF;
while(<INf>) {
	chomp;
	next if (m/^#/);
	my @a=split("\t",$_);
	if (($a[6] eq 0) and ($a[1] eq $a[2]) and ($a[1] > $cutoff*$a[3])){
        my $strand="+";
        if ($a[5] eq "+"){ $strand="-";}
		my $id=join("_",$a[4],$a[8],$a[7],$strand.abs($a[7]-$a[8]));
		# report seperately looped circles
		# newid-83501973/1__1     151     151     151     10      +       0       111901776       111901810       2____111901776____NA____111901863____NA____+    3____111901776____NA____111901810____NA____+    1____111901827____NA____111901863____NA____+
		my $last=-1;
		my $f=0;
		my $Nr=scalar(@a);
		for(my $i=9; $i<$Nr; $i++) {
			my @c=split(/\_\_\_\_/,$a[$i]);
			if ($last eq -1) { $last=$c[3]; }
			else {
				if ($last > $c[1]) { $f++; }
				$last=$c[3];
			}
		}
		if ($f > 0) {
			print OUTfm join("\t",@a),"\n";
			#next;
		}
		
		if (exists $uniqf{$id}) {
			$uniqf{$id}=$uniqf{$id}.",".$a[0];
		}
		else {
			$uniqf{$id}=$a[0];
		}
		$uniqF{$a[0]}=join("\t",@a);
		
        if ($a[5] eq "-") {
            # R-m+
            my @tmpLeft=split(/\_\_\_\_/,$a[9]);
            my @tmpRight=split(/\_\_\_\_/,$a[-1]);
            $Overlap{$a[0]}=$tmpLeft[6]+$tmpLeft[7]-$tmpRight[6];
        }
        else{
            # R+m-
            my @tmpLeft=split(/\_\_\_\_/,$a[9]);
            my @tmpRight=split(/\_\_\_\_/,$a[-1]);
            $Overlap{$a[0]}=$tmpRight[6]+$tmpRight[7]-$tmpLeft[6];
        }
	}
}

my %candy;
foreach my $id (sort keys %uniqf) {
	my @a=split(",",$uniqf{$id});
	print OUTf $id,"\t",scalar(@a),"\t",$uniqf{$id},"\n";
	my %POS;
	my %break_left;
	my %break_right;
	my $f=0;
	my $strand=0;
	my $min_gap=999;
	my $gstart=999999999;
	my $gend=-1;
	for(@a) {
		my $tmpid=$_;
		my @b=split("\t",$uniqF{$tmpid});
		if ($b[5] eq "+") { $strand++; }
		else {$strand--;}
		if ($min_gap > $Overlap{$tmpid}) { $min_gap=$Overlap{$tmpid}; }
		
		my $Nr=scalar(@b);
		my $last_rank=-1;
		my $last_pos=-1;
		if ($f eq 0) {  # deal with the first seq
			for(my $i=9; $i<$Nr; $i++) {
				my @c=split(/\_\_\_\_/,$b[$i]);
				if ($c[1] < $gstart) { $gstart=$c[1];}
				if ($c[1] > $gend) { $gend=$c[1];}
				$POS{$c[1]}=$c[3];
				$f++;
				if ($last_rank eq -1) {	}
				else {
					if(abs($last_rank - $c[0]) > 1) {
						$break_left{$last_pos}=1;
						$break_right{$c[1]}=1;
					}
					else {
						$break_left{$last_pos}=2;
						$break_right{$c[1]}=2;
					}
				}
				$last_rank=$c[0];
				$last_pos=$c[1];
			}
		}
		else {  # update with following seqs
			for(my $i=9; $i<$Nr; $i++) {
				my @c=split(/\_\_\_\_/,$b[$i]);
				if ($c[1] < $gstart) { $gstart=$c[1];}
				if ($c[1] > $gend) { $gend=$c[1];}
				my $overlap=0;
				my $this_left;
				my $this_right;
				foreach my $left (sort{$a <=> $b} keys %POS) {
					my $right=$POS{$left};
					if (($c[3] < $left) or ($right < $c[1])) {
						# no-overlap
					}
					else {
						$overlap=$left;
						last;
					}
				}
				my $known=0;
				if ($debug eq 10) { print $b[$i],"\n",$overlap,"\t",$POS{$overlap},"\n"; }
				if ($overlap eq 0) {
					$POS{$c[1]}=$c[3];
					$f++;
					$this_left=$c[1];
					$this_right=$c[3];
				}
				else {
					$known=$overlap;
					$this_left=$c[1] <= $overlap ? $c[1] : $overlap;
					$this_right=$c[3] <= $POS{$overlap} ? $POS{$overlap} : $c[3];
					if (($this_left ne $overlap) or ($this_right ne $POS{$overlap})) {
						if (exists $break_left{$overlap}) {
							my $tmp=$break_left{$overlap};
							delete $break_left{$overlap};
							$break_left{$this_left}=$tmp;
						}
						if (exists $break_right{$overlap}) {
							my $tmp=$break_right{$overlap};
							delete $break_right{$overlap};
							$break_right{$this_left}=$tmp
						}
						delete $POS{$overlap};
						$POS{$this_left}=$this_right;
					}
					$overlap=0;
					# check if this expansion would overlap with the next exon
					foreach my $left (sort{$a <=> $b} keys %POS) {
						if ($left ne $this_left) {
							my $right=$POS{$left};
							if (($c[3] < $left) or ($right < $c[1])) {
								# no-overlap
							}
							else {
								$overlap=$left;
								last;
							}
						}
					}
					if ($overlap ne 0) {
						if ($debug eq 10) { print "Another merge","\n",$this_left,"\t",$this_right,"\n",$overlap,"\t",$POS{$overlap},"\n"; }
						# merge exon {$this_left, $this_right} and {$overlap, $POS{$overlap}}
						my $new_left=$this_left <= $overlap ? $this_left : $overlap;
						my $new_right=$this_right <= $POS{$overlap} ? $POS{$overlap} : $this_right;
						$break_left{$this_left}=2;
						$break_right{$overlap}=2;
						delete $POS{$this_left};
						delete $POS{$overlap};
						$POS{$new_left}=$new_right;
						if (exists $break_right{$this_left}) {
							my $tmp=$break_right{$this_left};
							delete $break_right{$this_left};
							$break_right{$new_left}=$tmp;
						}
						if (exists $break_left{$overlap}) {
							my $tmp=$break_left{$overlap};
							delete $break_left{$overlap};
							$break_left{$new_left}=$tmp;
						}
						$this_left=$new_left;
						$this_right=$new_right;
					}
				}
				if ($last_rank eq -1) {
					$last_rank=$c[0];
					$last_pos=$this_left;
				}
				else {
					if(abs($last_rank - $c[0]) > 1) {
						if (($known eq 0) and ((!exists $break_left{$last_pos}) or ($break_left{$last_pos} ne 2)) and ($last_pos ne $gend)) { $break_left{$last_pos}=1; }
						if (($known eq 0) and ((!exists $break_right{$this_left}) or ($break_right{$this_left} ne 2)) and ($this_left ne $gstart)) { $break_right{$this_left}=1; }
					}
					else {
						$break_left{$last_pos}=2;
						$break_right{$this_left}=2;
						if ($known eq 0) {
							if (((!exists $break_left{$this_left}) or ($break_left{$this_left} ne 2)) and ($this_left ne $gend)){
								$break_left{$this_left}=1
							}
						}
						
					}
					$last_rank=$c[0];
					$last_pos=$this_left;
				}
				if ($debug eq 10) {
					foreach my $ttt (sort keys %break_left) { print join("\t","break_left",$ttt,$break_left{$ttt}),"\n"; }
					foreach my $ttt (sort keys %break_right) { print join("\t","break_right",$ttt,$break_right{$ttt}),"\n"; }
					foreach my $ttt (sort keys %POS) { print join("\t","POS",$ttt,$POS{$ttt}),"\n"; }
				}
			}
		}
		if ($debug eq 10) {
			print $uniqF{$tmpid},"\n";
			my $debug=$id."\t".$tmpid."\t".$f;
			foreach my $left (sort{$a <=> $b} keys %POS) { $debug=$debug."\t".$left."\t".$POS{$left}; }
			print $debug,"\n";
			foreach my $ttt (sort keys %break_left) { print join("\t","break_left",$ttt,$break_left{$ttt}),"\n"; }
			foreach my $ttt (sort keys %break_right) { print join("\t","break_right",$ttt,$break_right{$ttt}),"\n"; }
			print join("\t",$gstart,$gend),"\n";
			print "\n";
		}
	}
	# picked the most used borders instead of simply the most-outer ones
	my %Freq;
	my %POSs;
	foreach my $left (sort{$a <=> $b} keys %POS) {
		my $right=$POS{$left};
		for(@a) {
			my $tmpid=$_;
			my @b=split("\t",$uniqF{$tmpid});
			my $Nr=scalar(@b);
			for(my $i=9; $i<$Nr; $i++) {
				my @c=split(/\_\_\_\_/,$b[$i]);
				if (($left <= $c[1]) and ($c[3] <= $right)) {
					$Freq{$left}{$c[1]}++;
					$Freq{$right}{$c[3]}++;
				}
			}
		}
		my $new_left=-1;
		foreach my $i (sort{$Freq{$left}{$b} <=> $Freq{$left}{$a}} keys %{$Freq{$left}}) {
			if ($new_left eq -1) { $new_left=$i; }
			elsif (($Freq{$left}{$i} >= $Freq{$left}{$new_left}) and ($i < $new_left)) { $new_left=$i; }
		}
		my $new_right=-1;
		foreach my $i (sort{$Freq{$right}{$b} <=> $Freq{$right}{$a}} keys %{$Freq{$right}}) {
			if ($new_right eq -1) { $new_right=$i; }
			elsif (($Freq{$right}{$i} >= $Freq{$right}{$new_right}) and ($i > $new_right)) { $new_right=$i; }
		}
        if (($debug eq 100) and (($left ne $new_left) or ($right ne $new_right))){
            print $id,"\t",$left,"\t",$right,"\t",$new_left,"\t",$new_right,"\t",join("\t",@a),"\n";
        }
        
		$POSs{$new_left}=$new_right;
		if (exists $break_left{$left}) { 
			my $tmp=$break_left{$left};
			delete $break_left{$left};
			$break_left{$new_left}=$tmp;
		}
		if (exists $break_right{$left}) { 
			my $tmp=$break_right{$left};
			delete $break_right{$left};
			$break_right{$new_left}=$tmp;
		}
	}
	if ($debug eq 10) {
		foreach my $ttt (sort keys %break_left) { print join("\t","break_left",$ttt,$break_left{$ttt}),"\n"; }
		foreach my $ttt (sort keys %break_right) { print join("\t","break_right",$ttt,$break_right{$ttt}),"\n"; }
	}
	
	# close gap when gap is shorter than 100, then NOT extend the break-points by 50 on both side
	my $last=-1;
	foreach my $left (sort{$a <=> $b} keys %POSs) {
		my $right=$POSs{$left};
		if ($last eq -1) { $last=$left; }
		else {
			if ((exists $break_left{$last}) and (exists $break_right{$left}) and ($break_left{$last} eq 1) and ($break_right{$left} eq 1)) {
				if (($left - $POSs{$last}) < 2*$closegap) {
					$break_left{$last}=2;
					if (exists $break_left{$left}) { $break_left{$last}=$break_left{$left};	delete $break_left{$left};}
					$break_right{$left}=2;
					delete $POSs{$left};
					delete $POSs{$last};
					$POSs{$last}=$right;
				}
				else {
					#$break_right{$left-$closegap}=$break_right{$left};
					#$break_right{$left}=2;
					#$POSs{$last}+=$closegap;
					#$POSs{$left-$closegap}=$POSs{$left};
					#delete $POSs{$left};
				}
			}
			$last=$left; 
		}
	}
	
	my @ID=split(/\_/,$id);
	my $break="|";
	my $blockSizes="";
	my $blockStarts="";
	my $blockcounts=0;
	$gstart=999999999;
    $gend=-1;
	foreach my $left (sort{$a <=> $b} keys %POSs) {
		if ((exists $break_left{$left}) and ($break_left{$left} eq 1)) { $break=$break.$POSs{$left}."|"; }
		if ((exists $break_right{$left}) and ($break_right{$left} eq 1)) { $break=$break.$left."|"; }
		if ($blockStarts eq "") {
			$blockStarts=0;
			$blockSizes=$POSs{$left} - $left;
		}
		else {
			$blockStarts=$blockStarts.",".($left - $ID[2]);
			$blockSizes=$blockSizes.",".($POSs{$left} - $left);
		}
		$blockcounts++;
		if ($left < $gstart) { $gstart=$left;}
        if ($POSs{$left} > $gend) { $gend=$POSs{$left};}
	}
	
	my $report=join("\t","chr".$ID[0],$gstart-1,$gend-1,$id."_".$break,scalar(@a));
	if ($strand > 0) { $report=$report."\t-\t".($gstart-1)."\t".($gstart-1)."\t0\t".$blockcounts; }
	else { $report=$report."\t+\t".($gstart-1)."\t".($gstart-1)."\t0\t".$blockcounts; }
	$report=$report."\t".$blockSizes."\t".$blockStarts;
	
	if (($gend - $gstart) < 1000000) { print OUTf12 $report,"\n"; }
	
	$report=join("\t",$id,$id."_".$break,$ID[0]);
	if ($strand > 0) { $report=$report."\t-";}
	else { $report=$report."\t+"; }
	$report=join("\t",$report,$gstart,$gend,$gstart,$gstart,$blockcounts);
	my $ExonL="";
	my $ExonR="";
	foreach my $left (sort{$a <=> $b} keys %POSs) {
		if ($ExonL eq "") {
			$ExonL=$left;
			$ExonR=$POSs{$left};
		}
		else {
			$ExonL=$ExonL.",".$left;
			$ExonR=$ExonR.",".$POSs{$left};
		}
	}
	$report=join("\t",$report,$ExonL,$ExonR,$min_gap);
	#print OUTf1 $report,"\n";
	
	$candy{$ID[0]}{$gstart}{$gend}{$strand}=join("\t",$id,$strand,$blockcounts,$ExonL,$ExonR,$min_gap,$break,scalar(@a),join(",",@a));
	
}
close INf;
close OUTf;

my %reported;
my $f=0; 
foreach my $chr (sort keys %candy) {
    if ($f eq 1) {print "\t=======\tyes\n";}
    elsif ($f eq -1) {print "\t=======\tno\n";}
    $f=-1;
    my $ChrFile=$genome."/".$chr.".fa";
    if (-e $ChrFile) { $f=$ChrFile; }
    $ChrFile=$genome."/chr".$chr.".fa";
    if (-e $ChrFile) { $f=$ChrFile; }
    $ChrFile=$genome."/chromosome".$chr.".fa";
    if (-e $ChrFile) { $f=$ChrFile; }
    if ($f ne -1) {
        $ChrFile=$f;
        $f=1;
    } else {next;}
    print "searching chromosome : $chr : ".$ChrFile;
    if (exists $reported{$chr}) {next;}
    open (IN1, $ChrFile) or next;
    <IN1>;
    my $SEQ;
    $f=1;
    $reported{$chr}=1;
    while(<IN1>) {chomp; $SEQ=$SEQ.$_;}
    close IN1;
	foreach my $start (sort keys %{$candy{$chr}}) {
        foreach my $end (sort keys %{$candy{$chr}{$start}}) {
            foreach my $mystrand (sort keys %{$candy{$chr}{$start}{$end}}) {
            my @a=split("\t",$candy{$chr}{$start}{$end}{$mystrand});
			# R+m-
			if ($a[1] > 0) {
                my %SScore5;
				my %SScore3;
                my $overlap=$a[5];
                my $left_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
                if ((length($left_seq) ne (4*$Extend + 1)) and ($debug eq 20)) {print join("\t",$a[0],$chr,$start,"left_404"),"\n"; next;}
				if ((length($left_seq) ne (4*$Extend + 1))) { next;}
                my $right_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
                if (($overlap > $Extend) and ($debug eq 20)){ print join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
                if ($overlap > $Extend) { next; }
                if ((length($right_seq) ne (4*$Extend + 1)) and ($debug eq 20)) {print join("\t",$a[0],$chr,$start,$end,"right_404"),"\n";next;}
				if ((length($right_seq) ne (4*$Extend + 1))) { next;}
                my $SMS=0;	
                my $PMS=0;
                my $SMscoren=-99999;
                if ($debug eq 20){ print "\n"; }
                if ($debug eq 20){ print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n";}
                for(my $k=0; $k<=$overlap; $k++) {
                    my $ml="";
                    my $mr="";
                    $ml=substr($left_seq,2*$Extend+1-$k,2);
                    $mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
                    if (($ml ne "") and ($mr ne "")) {
                        my $motif=$ml.$mr;
                        if ($debug eq 20){print $k,"\t",$motif,"\n";}
                        if ((exists $SMotif{$motif}) and ($SMotif{$motif} < 0)) {
                            my $t=$left_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_left_seq=scalar reverse $t;		
                            $t=$right_seq;					
                            $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                            my $rc_right_seq=scalar reverse $t;		
                            my $str1="";
                            $str1=uc substr($rc_left_seq,(2*$Extend+$k-20),23);
                            my $str2="";
                            $str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$k),9);
                            if ((length($str1) eq 23) and (length($str2) eq 9)) {
                                my $left_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
                                my $right_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
                                my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                                if ($SMscoren < $sumt) { $SMscoren=$sumt; $SMS=$SMotif{$motif}; $PMS=$k; }
                                $SScore5{$k}=$right_maxt;
                                $SScore3{$k}=$left_maxt;
                            }
                            if ($debug eq 20){print $str1,"\t",$str2,"\n";}
                        }
                    }
                }
                if ($debug eq 20){print $SMS,"\t",$PMS,"\n";}
                if ($SMS < 0) {}
                else {
                    $SMscoren=-99999;
                    for(my $k=0; $k<=$overlap; $k++) {
                        my $t=$left_seq;					
                        $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                        my $rc_left_seq=scalar reverse $t;		
                        $t=$right_seq;					
                        $t=~tr/[atcgATCG]/[TAGCTAGC]/;
                        my $rc_right_seq=scalar reverse $t;		
                        my $str1="";
                        $str1=uc substr($rc_left_seq,(2*$Extend+$k-20),23);
                        my $str2="";
                        $str2=uc substr($rc_right_seq,(2*$Extend+1-3-$overlap+$k),9);
                        if ((length($str1) eq 23) and (length($str2) eq 9)) {
                            my $left_maxt=sprintf("%.2f", &log2(&scoreconsensus3($str1)*&maxentscore(&getrest3($str1),\@metables)));
                            my $right_maxt=sprintf("%.2f",&log2(&scoreconsensus5($str2)*$me2x5{$seq{&getrest5($str2)}}));
                            my $sumt=sprintf("%.2f",($left_maxt + $right_maxt));
                            if ($SMscoren < $sumt) {
                                $SMscoren=$sumt;
                                $PMS=$k;
                            }
                            $SScore5{$k}=$right_maxt;
                            $SScore3{$k}=$left_maxt;
                        }    
                    }
                    if ($debug eq 20){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
                }
                
				my $Start=$start;
				my $End=$end;
				my @ExonL=split(/\,/,$a[3]);
				my @ExonR=split(/\,/,$a[4]);
				my @ID=split(/\_/,$a[0]);
				my $report=join("\t",$a[0]."_".$a[6],$a[0],"chr".$ID[0]);
				if ($a[1] > 0) { $report=$report."\t-";}
				else { $report=$report."\t+"; }
				$report=join("\t",$report,$Start,$End,$Start,$Start,$a[2]);
                my $tmpPMS=0;
                if (!exists $SScore5{$tmpPMS}) { $SScore5{$tmpPMS}=-99; $SScore3{$tmpPMS}=-99; }
				$report=join("\t",$report,join(",",@ExonL),join(",",@ExonR),$a[5],0,$SScore5{$tmpPMS},0,$SScore3{$tmpPMS},($SScore5{$tmpPMS}+$SScore3{$tmpPMS}));
				print OUTf1 $report,"\n";
				my $fasta="";
				for(my $k=0; $k<$a[2]; $k++) {
					if ($k eq 0) {
						my $t=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k]));
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_t=scalar reverse $t;
						$fasta=$rc_t;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k]-1,$ExonR[$k]-1,$rc_t),"\n";}
					}
					else {
						my $t=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k]));
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_t=scalar reverse $t;
						$fasta=$rc_t.$fasta;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k]-1,$ExonR[$k]-1,$rc_t),"\n";}
					}
				}
				print OUTfa1 ">".$a[0],"\n",$fasta,"\n";

				$Start=$start+($overlap-$PMS);
				$End=$end-$PMS;
				my $fixedid=join("_",$ID[0],$End,$Start,"-".abs($Start-$End));	# replace $a[0]
				#$report=join("\t",$a[0]."_".$a[6],$a[0],"chr".$ID[0]);
				$report=join("\t",$fixedid."_".$a[6],$fixedid,"chr".$ID[0]);
				if ($a[1] > 0) { $report=$report."\t-";}
				else { $report=$report."\t+"; }
				$report=join("\t",$report,$Start,$End,$Start,$Start,$a[2]);
				$ExonL[0]=$Start;
				$ExonR[-1]=$End;
				$report=join("\t",$report,join(",",@ExonL),join(",",@ExonR),($SScore5{$PMS}+$SScore3{$PMS}),$SScore5{$PMS},$SScore3{$PMS},$a[5],($a[5]-$PMS),$PMS,$a[-1]);
				print OUTf2 $report,"\n";
				my $blockSizes="";
				my $blockStarts="";
				$fasta="";
				for(my $k=0; $k<$a[2]; $k++) {
					if ($blockStarts eq "") {
						$blockStarts=0;
						$blockSizes=$ExonR[$k] - $ExonL[$k];
						my $t=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k]));
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_t=scalar reverse $t;
						$fasta=$rc_t;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k]-1,$ExonR[$k]-1,$rc_t),"\n";}
					}
					else {
						$blockStarts=$blockStarts.",".($ExonL[$k] - $Start);
						$blockSizes=$blockSizes.",".($ExonR[$k] - $ExonL[$k]);
						my $t=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k]));
						$t=~tr/[atcgATCG]/[TAGCTAGC]/;
						my $rc_t=scalar reverse $t;
						$fasta=$rc_t.$fasta;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k]-1,$ExonR[$k]-1,$rc_t),"\n";}
					}
				}
				$report=join("\t","chr".$ID[0],$Start-1,$End-1,$fixedid."_".$a[6]."_".($SScore5{$PMS}+$SScore3{$PMS}),$a[7],"-",$Start-1,$Start-1,0,$a[2],$blockSizes,$blockStarts);
				if (($End - $Start) < 1000000) { print OUTf22 $report,"\n"; }
				print OUTfa2 ">".$fixedid,"\n",$fasta,"\n";
			}
			# R-m+
			else {
				my %SScore5;
				my %SScore3;
                my $overlap=$a[5];
                my $left_seq=substr($SEQ,($end-2*$Extend-1),(4*$Extend + 1));
                if ((length($left_seq) ne (4*$Extend + 1)) and ($debug eq 20)) {print join("\t",$a[0],$chr,$start,"left_404"),"\n"; next;}
				if ((length($left_seq) ne (4*$Extend + 1))) { next;}
                my $right_seq=substr($SEQ,($start-2*$Extend-1),(4*$Extend + 1));
                if (($overlap > $Extend) and ($debug eq 20)){ print join("\t",@a),"\tlarge_overlap\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next; }
                if ($overlap > $Extend) { next; }
                if ((length($right_seq) ne (4*$Extend + 1)) and ($debug eq 20)) {print join("\t",$a[0],$chr,$start,$end,"right_404"),"\n";next;}
				if ((length($right_seq) ne (4*$Extend + 1))) { next;}
                if ((sanity($left_seq) > 0) or (sanity($right_seq) > 0)) { print OUTerr join("\t",@a),"\tsanity_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
                if ((length($left_seq) ne (4*$Extend+1)) or (length($right_seq) ne (4*$Extend+1))) { print OUTerr join("\t",@a),"\tlength_fail\t",$left_seq,"\t",$right_seq,"\t",$overlap,"\n"; next;}
                my $SMS=0;	
                my $PMS=0;
                my $SMscoren=-99999;
                if ($debug eq 20){ print "\n"; }
                if ($debug eq 20){print join("\t",@a),"\n",$left_seq,"\t",$right_seq,"\n";}
                
                for(my $k=0; $k<=$overlap; $k++) {
                    my $ml="";
                    my $mr="";
                    $ml=substr($left_seq,2*$Extend+1-$k,2);
                    $mr=substr($right_seq,2*$Extend-2+$overlap-$k,2);
                    if (($ml ne "") and ($mr ne "")) {
                        my $motif=$ml.$mr;
                        if ($debug eq 20){print $k,"\t",$motif,"\n";}
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
                                $SScore5{$k}=$right_maxt;
                                $SScore3{$k}=$left_maxt;
                            }
                            if ($debug eq 20){print $str1,"\t",$str2,"\n";}
                        }
                    }
                }
                if ($debug eq 20){print $SMS,"\t",$PMS,"\n";}
                if ($SMS > 0) {}
                else {
                    $SMscoren=-99999;
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
                            }
                            $SScore5{$k}=$right_maxt;
                            $SScore3{$k}=$left_maxt;
                        }    
                    }
                    if ($debug eq 20){print $SMscoren,"\t",$PMS,"\t",($overlap-$PMS),"\t",$overlap,"\n";}
                }
                
				# refFlat format
				my @ExonL=split(/\,/,$a[3]);
				my @ExonR=split(/\,/,$a[4]);
				#for(my $k=0; $k<$a[2]; $k++) { $ExonR[$k]++; }
				my $Start=$start;
				my $End=$end;
				my @ID=split(/\_/,$a[0]);
				my $report=join("\t",$a[0]."_".$a[6],$a[0],"chr".$ID[0]);
				if ($a[1] > 0) { $report=$report."\t-";}
				else { $report=$report."\t+"; }
				$report=join("\t",$report,$Start,$End,$Start,$Start,$a[2]);
				$ExonL[0]=$Start;
				$ExonR[-1]=$End;
                my $tmpPMS=0;
                if (!exists $SScore5{$tmpPMS}) { $SScore5{$tmpPMS}=-99; $SScore3{$tmpPMS}=-99; }
				$report=join("\t",$report,join(",",@ExonL),join(",",@ExonR),$a[5],0,$SScore5{$tmpPMS},0,$SScore3{$tmpPMS},($SScore3{$tmpPMS}+$SScore5{$tmpPMS}));
				print OUTf1 $report,"\n";
				my $fasta="";
				for(my $k=0; $k<$a[2]; $k++) {
					if ($k eq 0) {
						my $tmpseq=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k] - 0));
						$fasta=$tmpseq;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k],$ExonR[$k],$tmpseq),"\n";}
					}
					else {
						my $tmpseq=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k] - 0));
						$fasta=$fasta.$tmpseq;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k],$ExonR[$k],$tmpseq),"\n";}
					}
				}
				print OUTfa1 ">".$a[0],"\n",$fasta,"\n";
                
                $Start=$start+($overlap-$PMS);
				$End=$end-$PMS;
				my $fixedid=join("_",$ID[0],$End,$Start,"+".abs($Start-$End));	# replace $a[0]
				#$report=join("\t",$a[0]."_".$a[6],$a[0],"chr".$ID[0]);
				$report=join("\t",$fixedid."_".$a[6],$fixedid,"chr".$ID[0]);
				if ($a[1] > 0) { $report=$report."\t-";}
				else { $report=$report."\t+"; }
				$report=join("\t",$report,$Start,$End,$Start,$Start,$a[2]);
				$ExonL[0]=$Start;
				$ExonR[-1]=$End;
                $report=join("\t",$report,join(",",@ExonL),join(",",@ExonR),($SScore5{$PMS}+$SScore3{$PMS}),$SScore5{$PMS},$SScore3{$PMS},$a[5],($a[5]-$PMS),$PMS,$a[-1]);
				print OUTf2 $report,"\n";
				
				my $blockSizes="";
				my $blockStarts="";
				$fasta="";
				for(my $k=0; $k<$a[2]; $k++) {
					if ($blockStarts eq "") {
						$blockStarts=0;
						$blockSizes=$ExonR[$k] - $ExonL[$k];
						my $tmpseq=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k] - 0));
						$fasta=$tmpseq;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k],$ExonR[$k],$tmpseq),"\n";}
					}
					else {
						$blockStarts=$blockStarts.",".($ExonL[$k] - $Start);
						$blockSizes=$blockSizes.",".($ExonR[$k] - $ExonL[$k]);
						my $tmpseq=substr($SEQ,$ExonL[$k]-1,($ExonR[$k] - $ExonL[$k] - 0));
						$fasta=$fasta.$tmpseq;
						if ($debug eq 30) {print join("\t",$k,$ExonL[$k],$ExonR[$k],$tmpseq),"\n";}
					}
				}
				$report=join("\t","chr".$ID[0],$Start-1,$End-1,$fixedid."_".$a[6]."_".($SScore5{$PMS}+$SScore3{$PMS}),$a[7],"+",$Start-1,$Start-1,0,$a[2],$blockSizes,$blockStarts);
				if (($End - $Start) < 1000000) { print OUTf22 $report,"\n"; }
				print OUTfa2 ">".$fixedid,"\n",$fasta,"\n";
			}
            }
		}
	}
}
if ($f eq 1) {print "\t=======\tyes\n";}
elsif ($f eq -1) {print "\t=======\tno\n";}



open(OUTFLAG,">Step2_MuSeg_finished");
print OUTFLAG "Step2_MuSeg_finished\n";
close OUTFLAG;


# flag = 110, could be reads spanning >=2 junctions that tophat2 failed to report
# e.g. newid-338078/1__46      151     151     151     5       110     43557213        43582556        1____43557213____NA____43557268____NA____+      2____43578770____NA____43578847____NA____+      3____43582533____NA____43582556____NA____+
#      TGGCAGGTAGGCTCAGGTGGTGAGCAATCCAGGGTACAGGCAACAGAGCAGGAACGTCATTTGGGGTTTGAGCTTCAGGTTACAGTCCATCATTAAAGGAAGTCTTGGCAGGAAGCTGGAAGCAGGAACTGTCACTTAGAAGAATGTCAGG

# flag = 1, could be back-slide
# e.g. newid-342383/1__45      76      76      76      7       1       25775554        25779181        2____25775554____NA____25775578____10___Atp1a3___ENSMUSG00000040907___19___31____+      3____25779151____10___Atp1a3___ENSMUSG00000040907___18___31____25779181____NA____+      1____25779156____10___Atp1a3___ENSMUSG00000040907___18___31____25779186____NA____+
#      ATTTGTTGGTCGAGTTGAACGGAATCTCAGTCTCATGGATGGATAGCTGGTATTTGTTGGTCGAGTTGAACGGAAT
#      newid-342383/1__45      76      7,+25779151,46S30M,33,0 7,+25779156,30M46S,33,0 7,+25775554,25S24M27S,14,0
