#!/usr/bin/perl -w
use strict;
die "Usage: $0  \"config_file\"   \"output_sh\"  " if (@ARGV < 2);
my $filein=$ARGV[0];
my $fileout=$ARGV[1];
my %SPEC;
open(IN,$filein) or die "Cannot open config_file $filein";
while (<IN>) {
    chomp;
    my @a=split("\t",$_);
    if ((scalar(@a) < 2) or ($a[0] eq "")) { next; }
    $SPEC{$a[0]}=$a[1];
}
close IN;
my $bwa_seed_len=16;
my $bwa_min_score=20;
my $thread=1;
my $minJump=100;
my $maxJump=2500000;
my $minSSSum=10;
my $minSamplecnt=1;
my $minReadcnt=1;
my $MAS=20;
my $coverage=0.9;
my $Junc=6;
my $ER=0.05;
my $stranded="no";
my $pre_defined_circRNA="no";
my $search_trans_splicing="no";
my $blat_path="";
my $do_blat_search="no";
my $ts_coverage=$coverage;
my $ts_MAS=0;
my $ts_minSSSum=$minSSSum;
my $ts_maxSpan=$maxJump;
# check if all parameters are set
if (!exists $SPEC{"BWA_folder"}) { die "BWA_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"BWA_genome_Index"}) { die "BWA_genome_Index must by specified in the config_file $filein";}
if (!exists $SPEC{"BWA_genome_folder"}) { die "BWA_genome_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"ACF_folder"}) { die "ACF_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"CBR_folder"}) { die "CBR_folder must by specified in the config_file $filein";}
if (!exists $SPEC{"Agtf"}) { $SPEC{"Agtf"}="no"; print "No gene annotation file is provided. Providing annotation would enhance the performance.\n";}
if (!exists $SPEC{"UNMAP"}) { die "UNMAP must by specified in the config_file $filein";}
if (!exists $SPEC{"UNMAP_expr"}) { die "UNMAP_expr must by specified in the config_file $filein";}
if (!exists $SPEC{"Seq_len"}) { die "Seq_len must by specified in the config_file $filein";}
if (exists $SPEC{"Thread"}) { $thread=$SPEC{"Thread"}; }
if (exists $SPEC{"minJump"}) { $minJump=$SPEC{"minJump"}; }
if (exists $SPEC{"maxJump"}) { $maxJump=$SPEC{"maxJump"}; }
if (exists $SPEC{"minSplicingScore"}) { $minSSSum=$SPEC{"minSplicingScore"}; }
if (exists $SPEC{"minSampleCnt"}) { $minSamplecnt=$SPEC{"minSampleCnt"}; }
if (exists $SPEC{"minReadCnt"}) { $minReadcnt=$SPEC{"minReadCnt"}; }
if (exists $SPEC{"minMappingQuality"}) { $MAS=$SPEC{"minMappingQuality"}; }
if (exists $SPEC{"Coverage"}) { $coverage=$SPEC{"Coverage"}; }
if (exists $SPEC{"minSpanJunc"}) { $Junc=$SPEC{"minSpanJunc"}; }
if (exists $SPEC{"ErrorRate"}) { $ER=$SPEC{"ErrorRate"}; }
if (exists $SPEC{"Strandness"}) { $stranded=$SPEC{"Strandness"}; }
if (exists $SPEC{"pre_defined_circle_bed"}) { $pre_defined_circRNA=$SPEC{"pre_defined_circle_bed"}; }
if (exists $SPEC{"BWA_seed_length"}) { $bwa_seed_len=$SPEC{"BWA_seed_length"}; }
if (exists $SPEC{"BWA_min_score"}) { $bwa_min_score=$SPEC{"BWA_min_score"}; }
if (exists $SPEC{"Search_trans_splicing"}) { $search_trans_splicing=$SPEC{"Search_trans_splicing"}; }
if (exists $SPEC{"blat_path"}) { $blat_path=$SPEC{"blat_path"}; }
if (exists $SPEC{"blat_search"}) { $do_blat_search=$SPEC{"blat_search"}; }
if ($blat_path eq "") { $do_blat_search="no"; }
if (exists $SPEC{"trans_splicing_coverage"}) { $ts_coverage=$SPEC{"trans_splicing_coverage"}; }
if (exists $SPEC{"trans_splicing_minMappingQuality"}) { $ts_MAS=$SPEC{"trans_splicing_minMappingQuality"}; }
if (exists $SPEC{"trans_splicing_minSplicingScore"}) { $ts_minSSSum=$SPEC{"trans_splicing_minSplicingScore"}; }
if (exists $SPEC{"trans_splicing_maxSpan"}) { $ts_maxSpan=$SPEC{"trans_splicing_maxSpan"}; }


my $command="";
open(OUT, ">".$fileout) or die "Cannot open output_sh file $fileout";
print OUT "#!/bin/bash\n\n";
print OUT "date\n";
print OUT "#Step1\n";
print OUT "echo \"Step1 maping_the_unmapped_reads_to_genome Started\" \n";
$command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." ".$SPEC{"BWA_genome_Index"}." ".$SPEC{"UNMAP"}." \> unmap.sam";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step1.pl unmap.sam unmap.parsed $MAS $coverage";
print OUT $command,"\n";
#$command="perl ".$SPEC{"ACF_folder"}."/get_selected_fa_from_pool.pl unmap.parsed.UID ".$SPEC{"UNMAP"}." unmap.parsed.UID.fa";
$command="ln -s ".$SPEC{"UNMAP"}." unmap.parsed.UID.fa";
print OUT $command,"\n";
print OUT "echo \"Step1 maping_the_unmapped_reads_to_genome Finished\" \n\n\n";


print OUT "#Step2\n";
print OUT "date\n";
print OUT "echo \"Step2 find_circle_supporting_sequences Started\" \n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step2.pl ".$SPEC{"CBR_folder"}."/ unmap.parsed.2pp.S1 ".$SPEC{"BWA_genome_folder"}."/ unmap.parsed.2pp.S2";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step2_MuSeg.pl ".$SPEC{"CBR_folder"}."/ unmap.parsed.segs ".$SPEC{"Agtf"}." ".$SPEC{"BWA_genome_folder"}."/ unmap.parsed.segs 10";
print OUT $command,"\n";
print OUT "echo \"Step2 find_circle_supporting_sequences Finished\" \n\n\n";


print OUT "#Step3\n";
print OUT "date\n";
print OUT "echo \"Step3 define_circle Started\" \n";
if ($pre_defined_circRNA eq "no"){
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step3.pl unmap.parsed.2pp.S3 unmap.parsed.2pp.S2.sum";
    print OUT $command,"\n";
}
else {
    $command="perl ".$SPEC{"ACF_folder"}."/get_circRNA_from_bed.pl ".$pre_defined_circRNA;
    print OUT $command,"\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step3.pl unmap.parsed.2pp.S3 unmap.parsed.2pp.S2.sum pre_defined_circRNA.sum";
    print OUT $command,"\n";
}

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step3_MuSeg.pl unmap.parsed.2pp.S3 unmap.parsed.segs.S2";
print OUT $command,"\n";
$command="cat unmap.parsed.segs.S2.matched unmap.parsed.segs.S2.novel2 > unmap.parsed.2pp.S4";
print OUT $command,"\n";
print OUT "echo \"Step3 define_circle Finished\" \n\n\n";


print OUT "#Step4\n";
print OUT "date\n";
print OUT "echo \"Step4 annotate_select_and_make_pseudo_sequences_for_circles Started\" \n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4.pl unmap.parsed.2pp.S4 ".$SPEC{"Agtf"}." circle_candidates 0 $minJump $maxJump $minSSSum";
print OUT $command,"\n";

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4_MEA.pl circle_candidates_MEA ".$SPEC{"Agtf"}." circle_candidates_MEA";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_split_exon_border_biotype_genename.pl circle_candidates_MEA.gtf circle_candidates_MEA.agtf";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_seq_from_agtf.pl circle_candidates_MEA.agtf ".$SPEC{"BWA_genome_folder"}."/ circle_candidates_MEA.pseudo";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_pseudo_circle.pl circle_candidates_MEA.pseudo.gene.fa circle_candidates_MEA.CL ".$SPEC{"Seq_len"};
print OUT $command,"\n";

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step4_CBR.pl circle_candidates_CBR ".$SPEC{"Agtf"}." circle_candidates_CBR";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_seq_from_agtf.pl circle_candidates_CBR.agtf ".$SPEC{"BWA_genome_folder"}."/ circle_candidates_CBR.pseudo";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_pseudo_circle.pl circle_candidates_CBR.pseudo.gene.fa circle_candidates_CBR.CL ".$SPEC{"Seq_len"};
print OUT $command,"\n";

print OUT "#Step5\n";
print OUT "date\n";
print OUT "echo \"Step5 caliberate_the_expression_of_circles Started\" \n";
$command=$SPEC{"BWA_folder"}."/bwa index circle_candidates_MEA.CL";
print OUT $command,"\n";
$command=$SPEC{"BWA_folder"}."/bwa index circle_candidates_CBR.CL";
print OUT $command,"\n";


$command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." circle_candidates_MEA.CL unmap.parsed.UID.fa \> circle_candidates_MEA.sam";
print OUT $command,"\n";
$command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." circle_candidates_CBR.CL unmap.parsed.UID.fa \> circle_candidates_CBR.sam";
print OUT $command,"\n";

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl circle_candidates_MEA.sam circle_candidates_MEA.CL circle_candidates_MEA.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl circle_candidates_CBR.sam circle_candidates_CBR.CL circle_candidates_CBR.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
print OUT $command,"\n";

$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5m2.pl unmap.parsed.tmp circle_candidates_MEA circle_candidates_CBR ".$SPEC{"UNMAP_expr"};
print OUT $command,"\n";
print OUT "echo \"Step5 caliberate_the_expression_of_circles Finished\" \n\n\n";


print OUT "#Step6\n";
print OUT "echo \"Step6 generating_UCSC_Genome_Browser_files Started\" \n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step6.pl circle_candidates_MEA.refFlat circle_candidates_MEA.expr circle_candidates_MEA.bed12 circle_candidates_MEA $minSSSum $minSamplecnt $minReadcnt 0 0";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/ACF_Step6.pl circle_candidates_CBR.refFlat circle_candidates_CBR.expr circle_candidates_CBR.bed12 circle_candidates_CBR $minSSSum $minSamplecnt $minReadcnt 0 0";
print OUT $command,"\n";
$command="cat circ*bed12 | grep -v 'track name' | grep -v '^-' | cut -f4 | cut -d\\| -f1 > circle_candidates_id";
print OUT $command,"\n";
$command="head -n 1 circle_candidates_MEA.expr > tmpheader";
print OUT $command,"\n";
$command="cat circle_candidates_MEA.expr circle_candidates_CBR.expr | grep -v 'newid' | grep -v '^-' > circle_candidates_expr_1";
print OUT $command,"\n";
$command="perl ".$SPEC{"ACF_folder"}."/get_selected_from_pool_singleline.pl circle_candidates_id circle_candidates_expr_1 circle_candidates_expr_2";
print OUT $command,"\n";
$command="cat tmpheader circle_candidates_expr_2 > circle_candidates_expr";
print OUT $command,"\n";
$command="rm -rf tmpheader";
print OUT $command,"\n";
$command="rm -rf circle_candidates_expr_?";
print OUT $command,"\n";
print OUT "echo \"Step6 generating_UCSC_Genome_Browser_files Finished\" \n\n\n";

print OUT "date\n";


if ($search_trans_splicing eq "yes") {
    print OUT "\n#Extra Step\: finding trans_splicing events\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step1.pl ".$SPEC{"CBR_folder"}."/ unmap.parsed.tmp ".$SPEC{"BWA_genome_folder"}."/ unmap.trans.splicing $ts_coverage 15 $ts_MAS $maxJump";
    print OUT $command,"\n";
    $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step2.pl unmap.trans.splicing ".$SPEC{"BWA_genome_folder"}."/ unmap.trans.splicing ".$SPEC{"Seq_len"}." $ts_minSSSum ".$SPEC{"Agtf"};
    print OUT $command,"\n";
    if (($do_blat_search eq "yes") and (-e $SPEC{"BWA_genome_Index"})) {
        #increasing the -minScore parameter value beyond one-half of the query size has no further effect
        #my $newlen=int(1.3*($SPEC{"Seq_len"}));
        my $newlen=$SPEC{"Seq_len"};
        $command=$blat_path." -minScore=".$newlen." -noHead ".$SPEC{"BWA_genome_Index"}." unmap.trans.splicing.tsloci.fa unmap.trans.splicing.tsloci.psl";
        print OUT $command,"\n";
        #$command="perl -ne \'chomp; my \@a=split(\"\\t\",\$_); if(abs(\$a[11]-\$a[12]) > $newlen){print \$a[9],\"\\n\"\;}\' unmap.trans.splicing.tsloci.psl \> unmap.trans.splicing.tsloci.psl.badid";
        #print OUT $command,"\n";
        #$command="sort \-u unmap.trans.splicing.tsloci.psl.badid \> unmap.trans.splicing.tsloci.psl.badids";
        #print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_linear_id_from_psl.pl unmap.trans.splicing.tsloci.psl unmap.trans.splicing.tsloci.psl.badids";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_from_pool_singleline.pl unmap.trans.splicing.tsloci.psl.badids unmap.trans.splicing.tsloci unmap.trans.splicing.tsloci.good 0 -1";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_from_pool_singleline.pl unmap.trans.splicing.tsloci.psl.badids unmap.trans.splicing.tsloci.anno unmap.trans.splicing.tsloci.anno.good 0 -1";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/get_selected_fa_from_pool.pl unmap.trans.splicing.tsloci.psl.badids unmap.trans.splicing.tsloci.fa unmap.trans.splicing.tsloci.fa.good -1";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa index unmap.trans.splicing.tsloci.fa.good";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." unmap.trans.splicing.tsloci.fa.good unmap.parsed.UID.fa \> unmap.trans.splicing.sam";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl unmap.trans.splicing.sam unmap.trans.splicing.tsloci.fa unmap.trans.splicing.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step3.pl unmap.parsed.tmp unmap.trans.splicing ".$SPEC{"UNMAP_expr"};
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_fusion_circRNAs.pl fusion_circRNAs unmap.trans.splicing ".$ts_maxSpan." unmap.trans.splicing.expr";
        print OUT $command,"\n";
    }
    else {
        $command=$SPEC{"BWA_folder"}."/bwa index unmap.trans.splicing.tsloci.fa";
        print OUT $command,"\n";
        $command=$SPEC{"BWA_folder"}."/bwa mem -t ".$thread." -k ".$bwa_seed_len." -T ".$bwa_min_score." unmap.trans.splicing.tsloci.fa unmap.parsed.UID.fa \> unmap.trans.splicing.sam";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_Step5.pl unmap.trans.splicing.sam unmap.trans.splicing.tsloci.fa unmap.trans.splicing.p1 ".$SPEC{"Seq_len"}." $Junc $ER $stranded";
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_trans_splice_step3.pl unmap.parsed.tmp unmap.trans.splicing ".$SPEC{"UNMAP_expr"};
        print OUT $command,"\n";
        $command="perl ".$SPEC{"ACF_folder"}."/ACF_fusion_circRNAs.pl fusion_circRNAs unmap.trans.splicing ".$ts_maxSpan." unmap.trans.splicing.expr";
        print OUT $command,"\n";
    }

}

print OUT "date\n";


close OUT;

