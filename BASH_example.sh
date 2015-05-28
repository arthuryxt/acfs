#!/bin/bash

#Step1
echo "Step1 maping_the_unmapped_reads_to_genome Started" 
/home/username/bin/bwa037a//bwa mem -t 30 -k 16 -T 20 /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/BWAIndex/genome.fa UNMAP > unmap.sam
perl /home/username/ACF//ACF_Step1.pl unmap.sam unmap.parsed 30 0.9
perl /home/username/ACF//get_selected_fa_from_pool.pl unmap.parsed.UID UNMAP unmap.parsed.UID.fa
echo "Step1 maping_the_unmapped_reads_to_genome Finished" 


#Step2
echo "Step2 find_circle_supporting_sequences Started" 
perl /home/username/ACF//ACF_Step2.pl /home/username/CB_splice// unmap.parsed.2pp.S1 /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes// unmap.parsed.2pp.S2
perl /home/username/ACF//ACF_Step2_MuSeg.pl /home/username/CB_splice// unmap.parsed.segs /data/iGenome/mouse/Ensembl/NCBIM37/Annotation/Genes/Mus_musculus.NCBIM37.67_split_exon.gtf /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes// unmap.parsed.segs 10
echo "Step2 find_circle_supporting_sequences Finished" 


#Step3
echo "Step3 define_circle Started" 
perl /home/username/ACF//ACF_Step3.pl unmap.parsed.2pp.S3 unmap.parsed.2pp.S2.sum
perl /home/username/ACF//ACF_Step3_MuSeg.pl unmap.parsed.2pp.S3 unmap.parsed.segs.S2
echo "Step3 define_circle Finished" 


#Step4
echo "Step4 annotate_select_and_make_pseudo_sequences_for_circles Started" 
perl /home/username/ACF//ACF_Step4.pl unmap.parsed.2pp.S3 /data/iGenome/mouse/Ensembl/NCBIM37/Annotation/Genes/Mus_musculus.NCBIM37.67_split_exon.gtf circle_candidates 10 100 1000000 10
perl /home/username/ACF//ACF_Step4_MEA.pl circle_candidates_MEA /data/iGenome/mouse/Ensembl/NCBIM37/Annotation/Genes/Mus_musculus.NCBIM37.67_split_exon.gtf circle_candidates_MEA
perl /home/username/ACF//get_split_exon_border_biotype_genename.pl circle_candidates_MEA.gtf circle_candidates_MEA.agtf
perl /home/username/ACF//get_seq_from_agtf.pl circle_candidates_MEA.agtf /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes// circle_candidates_MEA.pseudo
perl /home/username/ACF//get_pseudo_circle.pl circle_candidates_MEA.pseudo.gene.fa circle_candidates_MEA.CL 76
perl /home/username/ACF//ACF_Step4_CBR.pl circle_candidates_CBR /data/iGenome/mouse/Ensembl/NCBIM37/Annotation/Genes/Mus_musculus.NCBIM37.67_split_exon.gtf circle_candidates_CBR
perl /home/username/ACF//get_seq_from_agtf.pl circle_candidates_CBR.agtf /data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes// circle_candidates_CBR.pseudo
perl /home/username/ACF//get_pseudo_circle.pl circle_candidates_CBR.pseudo.gene.fa circle_candidates_CBR.CL 76
perl /home/username/ACF//get_pseudo_circle.pl unmap.parsed.segs.S2.novel.fa circle_candidates_MuS.CL 76
echo "Step4 annotate_select_and_make_pseudo_sequences_for_circles Finished" 


#Step5
echo "Step5 caliberate_the_expression_of_circles Started" 
/home/username/bin/bwa037a//bwa index circle_candidates_MEA.CL
/home/username/bin/bwa037a//bwa index circle_candidates_CBR.CL
/home/username/bin/bwa037a//bwa index circle_candidates_MuS.CL
/home/username/bin/bwa037a//bwa mem -t 30 -k 16 -T 20 circle_candidates_MEA.CL unmap.parsed.UID.fa > circle_candidates_MEA.sam
/home/username/bin/bwa037a//bwa mem -t 30 -k 16 -T 20 circle_candidates_CBR.CL unmap.parsed.UID.fa > circle_candidates_CBR.sam
/home/username/bin/bwa037a//bwa mem -t 30 -k 16 -T 20 circle_candidates_MuS.CL unmap.parsed.UID.fa > circle_candidates_MuS.sam
perl /home/username/ACF//ACF_Step5.pl circle_candidates_MEA.sam circle_candidates_MEA.CL circle_candidates_MEA.p1 76 6 0.05 - 
perl /home/username/ACF//ACF_Step5.pl circle_candidates_CBR.sam circle_candidates_CBR.CL circle_candidates_CBR.p1 76 6 0.05 - 
perl /home/username/ACF//ACF_Step5.pl circle_candidates_MuS.sam circle_candidates_MuS.CL circle_candidates_MuS.p1 76 6 0.05 - 
perl /home/username/ACF//ACF_Step5m.pl unmap.parsed.tmp circle_candidates_MEA circle_candidates_MuS circle_candidates_CBR UNMAP_expr 30
echo "Step5 caliberate_the_expression_of_circles Finished" 


#Step6
echo "Step6 generating_UCSC_Genome_Browser_files Started" 
perl /home/username/ACF//ACF_Step6.pl circle_candidates_MEA.refFlat circle_candidates_MEA.p1.2.expr circle_candidates_MEA.bed12 circle_candidates_MEA 10 1 2 0 0
perl /home/username/ACF//ACF_Step6.pl circle_candidates_CBR.refFlat circle_candidates_CBR.p1.2.expr circle_candidates_CBR.bed12 circle_candidates_CBR 10 1 2 0 0
perl /home/username/ACF//ACF_Step6.pl unmap.parsed.segs.S2 circle_candidates_MuS.p1.2.expr circle_candidates_MuS.bed12 circle_candidates_MuS 10 1 2 0 0
echo "Step6 generating_UCSC_Genome_Browser_files Finished" 


