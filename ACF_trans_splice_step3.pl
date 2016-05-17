#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"parsed.tmp\"   \"unmap.trans.splicing\"   \"expr\"    \"\(optional\) min_AS\"" if (@ARGV < 3);
# use Scalar::Util qw(looks_like_number);
# to estimate the expresison of circs after parsed remap
my $filetmp=$ARGV[0];	# unmap.parsed.tmp
my $filein1=$ARGV[1];	# unmap.trans.splicing
my $filein4=$ARGV[2];	# UNMAP_expr
my $MAS=0;
if (scalar(@ARGV) > 3) {$MAS=$ARGV[3];}
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
    #if (exists $OK{$a[0]}) {
		if (exists $uniq1{$b[0]}) { $uniq1{$b[0]}=$uniq1{$b[0]}."\t".$a[0]; }
		else { $uniq1{$b[0]}=$a[0]; }
		$OK{$a[0]}++;
    #}
}
close IN1;
open IN11,$filein1.".tsloci.anno";
while(<IN11>) {
    chomp;
    my @a=split("\t",$_);
    $Gname{$a[0]}=$a[3]."|".$a[2];
}
close IN11;


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

open OUT1,">".$filein1.".expr";
open OUT11,">".$filein1.".newid";
print OUT1 join("\t",@Header),"\n";
foreach my $id (sort keys %uniq1) {
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
            #if (!looks_like_number($b[$i-1])) {print $id,"\t",$uniq1{$id},"\n",$Anno{$read},"\n";}
		}
    }
    print OUT1 join("\t",@info),"\n";
}
close OUT1;
close OUT11;


