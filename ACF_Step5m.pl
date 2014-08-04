#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"parsed.tmp\"   \"MEA\"    \"MuS\"   \"CBR\"   \"expr\"    \"\(optional\) min_AS\"" if (@ARGV < 5);
# to estimate the expresison of circs after parsed remap
my $filetmp=$ARGV[0];	# parsed.tmp
my $filein1=$ARGV[1];	# circle_candidates_MEA
my $filein2=$ARGV[2];	# circle_candidates_MuS
my $filein3=$ARGV[3];	# circle_candidates_CBR
my $filein4=$ARGV[4];	# UNMAP_expr
my $MAS=30;
if (scalar(@ARGV) > 5) {$MAS=$ARGV[5];}
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
		if ($b[3] > $MAS){
		    if($info ne ""){$info=$info."\t".join(",",@b);}
		    else{$info=join(",",@b)}
		    $ASsum+=$b[3];
		    $ASc++;
		}
    }
    if ($info ne "") {
		$OK{$a[0]}=1000;
		$OKinfo{$a[0]}=$a[0]."\t".int($ASsum/$ASc)."\t".$info;
    }
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
    if (exists $OK{$a[0]}) {
		if (exists $uniq1{$b[0]}) { $uniq1{$b[0]}=$uniq1{$b[0]}."\t".$a[0]; }
		else { $uniq1{$b[0]}=$a[0]; }
		$OK{$a[0]}++;
    }
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
    if ((exists $OK{$a[0]}) and ($OK{$a[0]} eq 1000)) {
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

my %uniq2;
open IN2,$filein2.".p1.2";
while(<IN2>) {
    chomp;
    my @a=split("\t",$_);
	my @b=split(/\_\_\_/,$a[2]);
    if ((exists $OK{$a[0]}) and ($OK{$a[0]} eq 1000)) {
		if (exists $uniq2{$b[0]}) { $uniq2{$b[0]}=$uniq2{$b[0]}."\t".$a[0]; }
		else { $uniq2{$b[0]}=$a[0]; }
		$OK{$a[0]}++;
		if (!exists $Gname{$b[0]}) { $Gname{$b[0]}="na"; }
    }
}
close IN2;

my %Anno;
open IN4, $filein4;
while(<IN4>) {
    chomp;
    my @a=split("\t",$_);
    if ($a[0] eq "newid") { my $tmpid=$a[0]; $a[0]=$a[0]."\tGname"; $Anno{$tmpid}=join("\t",@a);}
    elsif ((exists $OK{$a[0]}) and ($OK{$a[0]} > 1000)) {
		$Anno{$a[0]}=join("\t",@a);
    }
}

my $header=$Anno{"newid"};
my @Header=split("\t",$header);
my $Nr=scalar(@Header);
my $template=0;
for(my $i=2; $i<$Nr; $i++) {$template=$template."\t0";}

open OUT1,">".$filein1.".expr";
open OUT11,">".$filein1.".newid";
print OUT1 $header,"\n";
foreach my $id (sort keys %uniq1) {
    my @a=split("\t",$uniq1{$id});
    my @info=split("\t",$id."\t".$template);
    for(@a) {
		my $read=$_;
		print OUT11 $OKinfo{$read},"\n";
		my @b=split("\t",$Anno{$read});
		$info[1]=$Gname{$id};
		for(my $i=2; $i<$Nr; $i++) {
		    $info[$i]+=$b[$i-1];
		}
    }
    print OUT1 join("\t",@info),"\n";
}
close OUT1;
close OUT11;

open OUT2,">".$filein2.".expr";
open OUT21,">".$filein2.".newid";
print OUT2 $header,"\n";
foreach my $id (sort keys %uniq2) {
    my @a=split("\t",$uniq2{$id});
    my @info=split("\t",$id."\t".$template);
    for(@a) {
		my $read=$_;
		print OUT21 $OKinfo{$read},"\n";
		my @b=split("\t",$Anno{$read});
		$info[1]=$Gname{$id};
		for(my $i=2; $i<$Nr; $i++) {
		    $info[$i]+=$b[$i-1];
		}
    }
    print OUT2 join("\t",@info),"\n";
}
close OUT2;
close OUT21;

open OUT3,">".$filein3.".expr";
open OUT31,">".$filein3.".newid";
print OUT3 $header,"\n";
foreach my $id (sort keys %uniq3) {
    my @a=split("\t",$uniq3{$id});
    my @info=split("\t",$id."\t".$template);
    for(@a) {
		my $read=$_;
		print OUT31 $OKinfo{$read},"\n";
		my @b=split("\t",$Anno{$read});
		$info[1]=$Gname{$id};
		for(my $i=2; $i<$Nr; $i++) {
		    $info[$i]+=$b[$i-1];
		}
    }
    print OUT3 join("\t",@info),"\n";
}
close OUT3;
close OUT31;

