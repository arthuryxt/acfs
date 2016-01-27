========== 0. NOTE ==========
0. Contact arthur.yxt@gmail.com
1. ACF = Arthurian Circle Finder
2. ACF takes fasta and fastq, though the quality information is NOT used.
3. Fasta/q header (such as ">Default_header") MUST be converted to ACF style so that multiple samples can be processes in one run.
   Change ">Default_header" into ">Truseq_sample_Default_header", where "sample" are the names of your samples.
   Make sure there is No underline within the sample name.
   e.g. ">Truseq_ctrl.1_Default_header" and ">Truseq_AGO2KO.1_Default_header" are OK; ">Truseq_ctrl_1_Default_header" is BAD because ACF will register "ctrl" instead of "ctrl_1".
4. Assume reads are obtained by Directional Sequencing, and the reads are assumed to be reverse-complementary to mRNAs.
   So  if the reads are actually sense (the same 5'->3' direction as mRNA), please reverse-complement all reads.
   And if the reads are actually stransless or pair-ended, please run in parallel original reads and reverse-complemented reads.
   No pair-end information is used, as the exact junction-site must be supported by a single read. Seeing is believing.
5. Pipelining through ACF_MAKE.pl is highly recommended. To see the arguments needed in each perl script, simply run the script without any arguments.
6. Splicing Strength module is adpoted from MaxEnt module from here : http://genes.mit.edu/burgelab/maxent/download/. Many thanks to Christoph Burge lab.
7. For non-commercial purposes the code is released under the General Public License (version 3). See the file LICENSE for more detail. If you are interested in commercial use of this software contact the author.
8. Please kindly cite the following paper : Nat Neurosci 2015 Apr;18(4):603-10 [PMID: 25714049]
by Arthur 2013-11-01
========== 0. NOTE ==========


Scroll down to section "6. Example" for a real-world example.


========= 1. pipeline description ========
1. Map all Tophat2-unmapped-reads to genome using BWA, seperate :
    1-part
    2-part-same-chromosome-same-strand  => true positive
    2-part-same-chromosome-diff-strand
    >2-part-same-chromosome             => some might be true positive, if they are originated from short exons
    >=2-part-diff-chromosome            => true negative
2. Estimate splicing strength using Christoph Berge's method for "2-part-same-chromosome-same-strand". report:
    forward-splice
    back-splice and canonical splice-motif {[GT-AT],[GC-AG],[AT-AC]}, calculate strength using CB
    back-splice and non-canonical splice-motif, using Christoph Berge's method, slide to find the maximal score and border
   Try to rescue circles using ">2-part-same-chromosome"
3. Check if both splice sites of the cirRNAs are on known exon-borders. report:
    both known exon-border (termed : MEA. [Match with Existing Annotation. somehow more reliable])
    at least one splice site sits on unknown exon-border (termed : CBR. [Chris Berge Rescued])
4. Build gtf and pseudo-transcript for results from step3
5. Map all Tophat2-unmapped-reads to pseudo-transcipts and estimate the abundance of circRNAs
6. Generate bed track for visulization
========= 1. pipeline description ========


========== 2. pre-process ==========
#1# copy folders "CB_splice" and "ACF_code" to local server

#2# Change fasta/fastq header format
#MUST-DO#
perl change_fastq_header.pl SRR650317_1.fasta SRR650317_1.fa Truseq_SRR650317left
perl change_fastq_header.pl SRR650317_2.fasta SRR650317_2.fa Truseq_SRR650317right
#MUST-DO#

#3# merge sequences from multiple fasta/fastq files into one fasta file
perl Truseq_merge_unique_fa.pl UNMAP newid SRR650317_1.fa SRR650317_2.fa
#alternatively, if there are many files to merge, generate a file (named filelist for example) contains the full-path of each file
perl Truseq_merge_unique_filelist.pl UNMAP newid filelist
# the aggregated fasta file "UNMAP" will be generated together with the readcount file "UNMAP_expr"

#4# build BWA index, using verion 0.73a (support for higher versions will be added soon)
#   please use the included package OR download at http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.3a.tar.bz2/download
/bin/bwa073a/bwa index /data//iGenome/human/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa

#5# prepare for annotation
# use the gtf file from iGenome package    or    download ensembl gtf here : ftp://ftp.ensembl.org/pub/current_gtf/
perl get_split_exon_border_biotype_genename.pl /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/genes.gtf /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf
========== 2. pre-process ==========



========== 3. How_to_Find_Your_Circles ==========
# make tab-delimited SPEC file
# These 9 parameters are mandatory :

#1# where you put BWA aligner
perl -e 'print join("\t","BWA_folder","/home/bin/bwa037a/"),"\n";' > SPEC_example.txt
#2# where you put BWA genome index
perl -e 'print join("\t","BWA_genome_Index","/data/iGenome/mouse/Ensembl/NCBIM37/Sequence/BWAIndex/genome.fa"),"\n";' >> SPEC_example.txt
#3# where you put genome sequences. Note that there MUST be one fasta sequence per chromosome. A concatenated fasta file does NOT work. 
perl -e 'print join("\t","BWA_genome_folder","/data/iGenome/mouse/Ensembl/NCBIM37/Sequence/Chromosomes/"),"\n";' >> SPEC_example.txt
#4# where you put ACF scripts
perl -e 'print join("\t","ACF_folder","/home/ACF/"),"\n";' >> SPEC_example.txt
#5# where you put Splicing Strength calculation module
perl -e 'print join("\t","CBR_folder","/home/CB_splice/"),"\n";' >> SPEC_example.txt
#6# where you put processed gtf. This is different from the gtf file from Ensembl or UCSC website. See Point-5 in Session "0. pre-process"
perl -e 'print join("\t","Agtf","/data/iGenome/mouse/Ensembl/NCBIM37/Annotation/Genes/Mus_musculus.NCBIM37.67_split_exon.gtf"),"\n";' >> SPEC_example.txt
#7# where you put reads that cannot be mapped to either genome or transcriptome directly
perl -e 'print join("\t","UNMAP","UNMAP"),"\n";' >> SPEC_example.txt
#8# where you put expression file of the collasped reads
perl -e 'print join("\t","UNMAP_expr","UNMAP_expr"),"\n";' >> SPEC_example.txt
#9# the length of your sequencing reads
perl -e 'print join("\t","Seq_len","150"),"\n";' >> SPEC_example.txt

# These 13 parameters are optional:
#10# the minimum distance of a back-splice (default = 100). The smaller this value, the more likely you can find circles from short exons
perl -e 'print join("\t","minJump","100"),"\n";' >> SPEC_example.txt
#11# the maximum distance of a back-splice (default = 1000000). The larger this value, the more likely you can find circles from long genes
perl -e 'print join("\t","maxJump","1000000"),"\n";' >> SPEC_example.txt
#12# the minimum score for the sum of splicing strength at both splice site (default = 10). A value small than 6 might suggest the circle is not generated from splicing, providing it is true
perl -e 'print join("\t","minSplicingScore","10"),"\n";' >> SPEC_example.txt
#13# the minimum number of samples that detect any given circle (default = 1).
perl -e 'print join("\t","minSampleCnt","1"),"\n";' >> SPEC_example.txt
#14# the minimum number of reads (from all samples) that detect any given circle (default = 2).
perl -e 'print join("\t","minReadCnt","2"),"\n";' >> SPEC_example.txt
#15# the minimum mapping quality of any given sequence (default = 30). A value larger than 20 is recommended.
perl -e 'print join("\t","minMappingQuality","30"),"\n";' >> SPEC_example.txt
#16# the minimum percentage of any given read is aligned (default = 0.9). The larger the better.
perl -e 'print join("\t","Coverage","0.9"),"\n";' >> SPEC_example.txt
#17# the minimum number of bases reach beyond the back-splice-site (default = 6). The larger the less likely of false-positive.
perl -e 'print join("\t","minSpanJunc","6"),"\n";' >> SPEC_example.txt
#18# the maximum error rate for re-alignment (default = 0.05). The smaller the more better.
perl -e 'print join("\t","ErrorRate","0.05"),"\n";' >> SPEC_example.txt
#19# the strand information of sequencing (default = "no"). "+" if the reads are sense to transcripts; "-" if the reads are anti-sense to transcripts; "no" if the sequencing is un-stranded.
perl -e 'print join("\t","Strandness","no"),"\n";' >> SPEC_example.txt
#20# the number of threads used for BWA alignment (default = 30).
perl -e 'print join("\t","Thread","30"),"\n";' >> SPEC_example.txt
#21# pre-defined circle annotation in bed6 or bed12 format (default = "no", to increase sensitivity for lowly expressed circRNAs please include bed12 files of annotated one, e.g. merge the bed files from GSE61991)
perl -e 'print join("\t","pre_defined_circle_bed","no"),"\n";' >> SPEC_example.txt

# make Pipeline Bash file
perl ~/ACF/ACF_MAKE.pl SPEC_example.txt BASH_example.sh

# find circles
bash BASH_example.sh
========== 3. How_to_Find_Your_Circles ==========


========== 4. How_to_Visulize_Your_Circles ==========
# circles are stored in two files:
circle_candidates_MEA.bed12
circle_candidates_CBR.bed12
# the higher the value in 5th column, the more "likely" that circle is true
# the name in 4th column can be seperated by "|" into four segments: circle-ID, sum-of-Splicing-Strength, number-of-supporting-samples, number-of-supporting-reads 
# the sequence for the "middle exons" in almost every circle is hypothetically filled in using the annotations provided, true sequences could be determined by integrating mRNA-Seq data and validated by inward and outward PCRs.

# their expression per sample is stored here:
circle_candidates_MEA.expr
circle_candidates_CBR.expr
# the second column denote the name of the gene from which this circRNA is derived 
========== 4. How_to_Visulize_Your_Circles ==========


========== 5. Tutorial ==========
#_1_# pre-processing for sequencing reads

#1A# for single-end RNA-Seq only/ (Feasible but not a good idea in practice, since you don't want to map all reads again. Use unmapped reads instead.)
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX852%2FSRX852583/SRR1772422/SRR1772422.sra
fastq-dump.2 --fasta 0 --split-files SRR1772422.sra
perl change_fastq_header.pl SRR1772422.fasta SRR1772422.fa Truseq_HippoSyn
perl Truseq_merge_unique_fa.pl UNMAP newid SRR1772422.fa
# Note this sample is from mouse, therefore mouse annotation should be used.

#1B# for paired-end RNA-Seq only. (Feasible but not a good idea in practice, since you don't want to map all reads again. Use unmapped reads instead.)
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX218%2FSRX218203/SRR650317/SRR650317.sra
fastq-dump.2 --fasta 0 --split-files SRR650317.sra
perl change_fastq_header.pl SRR650317_1.fasta SRR650317_1.fa Truseq_SRR650317left
perl change_fastq_header.pl SRR650317_2.fasta SRR650317_2.fa Truseq_SRR650317right
perl Truseq_merge_unique_fa.pl UNMAP newid SRR650317_1.fa SRR650317_2.fa
# As acfs require sequencing reads to be stranded, reverse complement all reads to fake strandness
perl -ne 'chomp; if(m/^>/){s/^>//; print ">rc".$_,"\n";}else{my $t=$_; $t=~tr/[ATCG]/[TAGC]/; my $r=scalar reverse $t; print $r,"\n";}' UNMAP1 > UNMAP1.rc
cat UNMAP1 UNMAP1.rc > UNMAP
perl -e 'open IN,"UNMAP_expr1"; my $line=<IN>; print $line; while(<IN>){chomp; my @a=split("\t",$_); print join("\t",@a),"\n"; $a[0]="rc".$a[0]; print join("\t",@a),"\n";}' > UNMAP_expr
# Note this sample is from human, therefore human annotation should be used.

#1C# for single-end and paired-end RNA-Seq unmapped reads.
# First align raw RNA-Seq reads to some reference with tools (such as Bowtie/BWA/STAR), obtain a SAM file containing all unmapped reads.
# For example if you choose Tophat, please convert <unmapped.bam> to <unmapped.sam>
# Then run the following command:
perl convert_unmapped_SAM_to_fa_for_acfs.pl unmapped.sam newName newNAME.fa
# where <newName> is the new sample ID that is used to change the fasta/fastq header to ACF style. It could be "Ctrl12" or "TreatA.rep2", your choice. Just remember, there should be NO underscore as "_" in it.
# and <newNAME.fa> is the name of the converted file. Name it whatever you like.
perl Truseq_merge_unique_fa.pl UNMAP newid newNAME.fa

#_2_# prepare for annotation
# use the gtf provided in iGenome package    or   get Homo_sapiens.GRCh37.71.gtf from  ftp://ftp.ensembl.org/pub/current_gtf/
perl get_split_exon_border_biotype_genename.pl /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/genes.gtf /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf

#_3_# generate parameter file
perl -e 'print join("\t","BWA_folder","/home/bin/bwa073a/"),"\n";' > SPEC_example.txt
perl -e 'print join("\t","BWA_genome_Index","/data/iGenome/human/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","BWA_genome_folder","/data/iGenome/human/Ensembl/GRCh37/Sequence/Chromosomes/"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","ACF_folder","/home/ACF/"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","CBR_folder","/home/CB_splice/"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","Agtf","/data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","UNMAP","UNMAP"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","UNMAP_expr","UNMAP_expr"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","Seq_len","76"),"\n";' >> SPEC_example.txt

# if you don't want to change the default parameters, you don't need to run the belowing command-lines.
perl -e 'print join("\t","Thread","30"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minJump","100"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","maxJump","1000000"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minSplicingScore","0"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minSampleCnt","1"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minReadCnt","1"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minMappingQuality","20"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","Coverage","0.9"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","minSpanJunc","6"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","ErrorRate","0.05"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","Strandness","-"),"\n";' >> SPEC_example.txt
perl -e 'print join("\t","pre_defined_circle_bed","no"),"\n";' >> SPEC_example.txt

#_4_# generate pipeline
perl ACF_MAKE.pl SPEC_example.txt BASH_example.sh

#_5_# get your circRNAs
bash BASH_example.sh

#_6_# take a look at these :
# bed12 files for visualization of circRNAs, such as using G-Browser
circle_candidates_MEA.bed12      # high-confident circles that cross annotated boundarys of exon(s)
circle_candidates_CBR.bed12      # low-confident circles (could still be true)
# expression(readcounts) table for circRNAs, with circRNAs in rows and samples in columns
circle_candidates_MEA.expr       
circle_candidates_CBR.expr
========== 5. Tutorial ==========


========== 6. Example ==========
1) DIY from scratch
Let's assume you have downloaded SRR1772422.sra (It's a single-end RNA-Seq from the mouse hippocampus synaptosome sample, read-length=101)
All reads are aligned to mouse genome (and transcriptome) using Tophat2.
All unmapped reads are extracted from "unmapped.bam" into "unmapped.fasta".
2) Alternatively, you can unzip "test.tar.gz" and use the "unmapped.fasta" directly. To save space, only circRNA supporting reads are included here.
3) change fasta header
perl change_fastq_header.pl unmapped.fasta unmapped.fa Truseq_F4
4) make aggregated file
perl Truseq_merge_unique_fa.pl UNMAP newid unmapped.fa
5) locate the acfs-version of annotation. (Otherwise, please refer to Point 5 in section "2. pre-process")
6) generate parameter file. (Please refer to Point-3 in section "5. Tutorial")
7) generate pipeline
perl ACF_MAKE.pl SPEC_example.txt BASH_example.sh
8) identify circRNAs
bash BASH_example.sh
9) the newly generated "circle_candidates_???.expr" and "circle_candidates_???.bed12" should be the same as the ones with "old_" prefix
========== 6. Example ==========


========== Change Log ==========
Update on 2015-03-09 :
   Now ACF can include pre-defined circRNA annotations from a bed6 or bed12 file (and their authenticity will be checked, so please adjust (minJump, maxJump, minSplicingScore) accordingly ).
   This way, you can both predict novel circRNAs in your data and estimate the abundance of annotated circRNAs.
Update on 2015-08-11 :
   Corrected the Tutorial section in README, thanks to Zol.
Update on 2015-08-20 :
   Added support for paired-end reads.
Update on 2015-09-17 :
   Added a small real-world example.
Update on 2016-01-27
   Performance improvement and add scripts for simulation
Expecting a much faster standalone version soon...
========== Change Log ==========
