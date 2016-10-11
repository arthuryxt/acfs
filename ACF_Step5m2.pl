#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"parsed.tmp\"   \"MEA\"    \"CBR\"   \"expr\"    \"\(optional\) min_AS\"" if (@ARGV < 4);
# to estimate the expresison of circs after parsed remap
my $filetmp=$ARGV[0];	# parsed.tmp
my $filein1=$ARGV[1];	# circle_candidates_MEA
my $filein3=$ARGV[2];	# circle_candidates_CBR
my $filein4=$ARGV[3];	# UNMAP_expr    # if this argument is set to "no", then assume all read have read_count == 1 for the only one sample; so that merging reads is not needed.
my $MAS=0;
if (scalar(@ARGV) > 4) {$MAS=$ARGV[4];}
my %OK;
my %OKinfo;
open IN0, $filetmp;
while(<IN0>) {
    chomp;
    my @a=split("\t",$_);
    my $Nr=scalar(@a);
    my $info="";
    my $ASsum=0;
    my $ASc=0;
    for(my $i=2; $i<$Nr; $i++){
		my @b=split(/\,/,$a[$i]);
		#if ($b[3] > $MAS){		# the AS is useless here
		    if($info ne ""){$info=$info."\t".join(",",@b);}
		    else{$info=join(",",@b)}
		    $ASsum+=$b[3];
		    $ASc++;
		#}
    }
    if ($info ne "") {
		$OK{$a[0]}=0;
		$OKinfo{$a[0]}=$a[0]."\t".int($ASsum/$ASc)."\t".$info;
    }
    else { $OK{$a[0]}=0; } 
}
close IN0;

my %Gname;
# circle_candidates_MEA
my %uniq1;
open IN1,$filein1.".p1.2";
while(<IN1>) {
    chomp;
    my @a=split("\t",$_);
	my @b=split(/\_\_\_/,$a[2]);
    #if (exists $OK{$a[0]}) {
		if (exists $uniq1{$b[0]}) { $uniq1{$b[0]}=$uniq1{$b[0]}."\t".$a[0]; }
		else { $uniq1{$b[0]}=$a[0]; }
		$OK{$a[0]}++;
    #}
}
close IN1;
open IN11,$filein1.".refFlat";
while(<IN11>) {
    chomp;
    my @a=split("\t",$_);
    $Gname{$a[1]}=$a[0];
}
close IN11;

# circle_candidates_CBR
my %uniq3;
open IN3,$filein3.".p1.2";
while(<IN3>) {
    chomp;
    my @a=split("\t",$_);
	my @b=split(/\_\_\_/,$a[2]);
    if ($OK{$a[0]} eq 0) {
		if (exists $uniq3{$b[0]}) { $uniq3{$b[0]}=$uniq3{$b[0]}."\t".$a[0]; }
		else { $uniq3{$b[0]}=$a[0]; }
		$OK{$a[0]}++;
    }
}
close IN3;
open IN31,$filein3.".refFlat";
while(<IN31>) {
    chomp;
    my @a=split("\t",$_);
    $Gname{$a[1]}=$a[0];
}
close IN31;

my %Anno;
my $header="";
my @Header;
if (($filein4 ne "no") and (-e $filein4)) {
    open IN4, $filein4;
    $header=<IN4>;
    chomp $header;
    @Header=split("\t",$header);
    my $tmpid=$Header[0];
    $Header[0]=$Header[0]."\tGname";
    $header=join("\t",@Header);
    while(<IN4>) {
        chomp;
        my @a=split("\t",$_);
        $Anno{$a[0]}=join("\t",@a);
    }
    close IN4;
}
else {
    $header=join("\t","newid","Sample");
    @Header=split("\t",$header);
    $Header[0]=$Header[0]."\tGname";
    $header=join("\t",@Header);
    foreach my $id(keys %OK) { $Anno{$id}=join("\t",1,1); }
}

my $Nr=scalar(@Header);
my $template=0;
for(my $i=2; $i<=$Nr; $i++) {$template=$template."\t0";}

# report the minimal number of read from prediction, in case that the second-step alignment missed due to heuristics
my %total_MEA;
my %info_MEA;
open(INtmp1, $filein1.".refFlat");
while (<INtmp1>) {
    chomp;
    if (m/^#/) { next;}
    my @a=split("\t",$_);
    $total_MEA{$a[1]}=1;
    $info_MEA{$a[1]}=join("\t",@a);
}
close INtmp1;
open(INtmp11, $filein1);
while (<INtmp11>) {
    chomp;
    my @a=split("\t",$_);
    if (exists $total_MEA{$a[0]}) {
        my @c=split(/\,/,$a[-1]);
        my @info=split("\t",$a[0]."\t".$template);
        $info[1]=$a[6];
        for(@c){
            my $read=$_;
            if (exists $Anno{$read}) {
                my @b=split("\t",$Anno{$read});
                for(my $i=2; $i<=$Nr; $i++) { $info[$i]+=$b[$i-1]; }
            }
        }
        $total_MEA{$a[0]}=join("\t",@info);
    }
}
close INtmp11;

my %total_CBR;
my %info_CBR;
open(INtmp3, $filein3.".refFlat");
while (<INtmp3>) {
    chomp;
    if (m/^#/) { next;}
    my @a=split("\t",$_);
    $total_CBR{$a[1]}=1;
    $info_CBR{$a[1]}=join("\t",@a);
}
close INtmp3;
open(INtmp31, $filein3);
while (<INtmp31>) {
    chomp;
    my @a=split("\t",$_);
    if (exists $total_CBR{$a[0]}) {
        my @c=split(/\,/,$a[-1]);
        my @info=split("\t",$a[0]."\t".$template);
        $info[1]=$a[6];
        if (($a[6] eq "na") and ($a[12] ne "na")) { $info[1]=$a[12]; }
        for(@c){
            my $read=$_;
            if (exists $Anno{$read}) {
                my @b=split("\t",$Anno{$read});
                for(my $i=2; $i<=$Nr; $i++) { $info[$i]+=$b[$i-1]; }
            }
        }
        $total_CBR{$a[0]}=join("\t",@info);
    }
}
close INtmp31;


#my %expP;
#my $cntP="unmap.parsed.2pp.S4";
#if (-e $cntP) {
#    my $tmpcnt=0;
#    open(INP, $cntP);
#    while (<INP>) {
#        chomp;
#        my @a=split("\t",$_);
#        if ((!exists $total_MEA{$a[0]}) and (!exists $total_CBR{$a[0]})) { next; }
#        my @c=split(/\,/,$a[14]);
#        my @info=split("\t",$a[0]."\t".$template);
#        for(@c){
#            my $read=$_;
#            if (exists $Anno{$read}) {
#                my @b=split("\t",$Anno{$read});
#                for(my $i=2; $i<=$Nr; $i++) { $info[$i]+=$b[$i-1]; }
#            }
#        }
#        $expP{$a[0]}=join("\t",@info);
#        $tmpcnt++;
#    }
#    close INP;
#    print "reading $tmpcnt lines in unmap.parsed.2pp.S4 finished.\n";
#}

# report the minimal number of read from prediction, in case that the second-step alignment missed due to heuristics

open OUT1,">".$filein1.".expr";
open OUT11,">".$filein1.".newid";
print OUT1 $header,"\n";
foreach my $id (sort keys %total_MEA) {
    if (exists $uniq1{$id}) {
        my @a=split("\t",$uniq1{$id});
        my @info=split("\t",$id."\t".$template);
        for(@a) {
            my $read=$_;
            if (exists $OKinfo{$read}) { print OUT11 $OKinfo{$read},"\n"; }
            else  {print OUT11 $read,"\n"; }
            my @b=split("\t",$Anno{$read});
            $info[1]=$Gname{$id};
            for(my $i=2; $i<=$Nr; $i++) {
                $info[$i]+=$b[$i-1];
            }
        }
        if (exists $total_MEA{$id}) {
            my @org=split("\t",$total_MEA{$id});
            for(my $i=2; $i<=$Nr; $i++) {
                if ($info[$i] < $org[$i]) { $info[$i]=$org[$i];}
            }
        }
        print OUT1 join("\t",@info),"\n";
    }
    else{
        print OUT1 $total_MEA{$id},"\n";
    }
}
close OUT1;
close OUT11;


open OUT3,">".$filein3.".expr";
open OUT31,">".$filein3.".newid";
print OUT3 $header,"\n";
foreach my $id (sort keys %total_CBR) {
    if (exists $uniq3{$id}) {
        my @a=split("\t",$uniq3{$id});
        my @info=split("\t",$id."\t".$template);
        for(@a) {
            my $read=$_;
            if (exists $OKinfo{$read}) { print OUT31 $OKinfo{$read},"\n"; }
            else  {print OUT31 $read,"\n"; }
            my @b=split("\t",$Anno{$read});
            $info[1]=$Gname{$id};
            for(my $i=2; $i<=$Nr; $i++) {
                $info[$i]+=$b[$i-1];
            }
        }
        if (exists $total_CBR{$id}) {
            my @org=split("\t",$total_CBR{$id});
            for(my $i=2; $i<=$Nr; $i++) {
                if ($info[$i] < $org[$i]) { $info[$i]=$org[$i];}
            }
        }
        print OUT3 join("\t",@info),"\n";
    }
    else{
        print OUT3 $total_CBR{$id},"\n";
    }
}

close OUT3;
close OUT31;

