# ACFS
Accurate CircRNA Finder Suite. Discovering circRNAs from RNA-Seq data.


# Overview
CircRNAs are generated through splicing, or to be precise, back-splicing where the downstream splice donor attact an upstream splice acceptor. Identifying the exact site of back-splice lies in the heart of circRNA discovery.

ACFS first examines and pinpoints the back-splice site from RNA-Seq alignment using a maximal entropy model. The expression of circRNAs is estimated from a second round of alignment to the inferred pseudo circular sequences.

No prior knowledge of gene annotation is needed for ACFS, but reading the coordinates is far less interesting than reading gene names, so circRNAs are annotated using the gene annotation if provided.


# Change Log
- Update on 2016-05-06 :
   Added extension for fining trans-splicing evidence, which can be used to identify fusion-circRNAs
- Update on 2016-01-27 :
   Performance improvement and add scripts for simulation
- Update on 2015-09-17 :
   Added a small real-world example
- Update on 2015-08-20 :
   Added support for paired-end reads
- Update on 2015-08-11 :
   Corrected the Tutorial section in README, thanks to Zol
- Update on 2015-03-09 :
   Now ACF can include pre-defined circRNA annotations from a bed6 or bed12 file (and their authenticity will be checked, so please adjust (minJump, maxJump, minSplicingScore) accordingly ).
   This way, you can both predict novel circRNAs in your data and estimate the abundance of annotated circRNAs.
- First release on 2013-11-01


# Installation
Simply unpack the ACFS package.


# Requirement
- bwa-0.7.3a (included in the package, but you **_need_** to do "make")
- perl
- blat (not necessary, sometime it helps to rule out false positive fusion-circRNAs when there is a gene-duplication or gene-pseudogene dilemma)


# Pipeline scheme
1. Map all Tophat2-unmapped-reads to genome using BWA, seperate :
    - 1-part
    - 2-part-same-chromosome-same-strand  (contain circRNAs)
    - 2-part-same-chromosome-diff-strand  (possibly PolII backwalk)
    - \>2-part-same-chromosome             (contain circRNAs, if they are originated from short exons)
    - \>=2-part-diff-chromosome            (contain trans-splicing, or even fusion-circRNAs)
2. Estimate splicing strength using Christoph Berge's method for "2-part-same-chromosome-same-strand". report:
    - forward-splice (not interesting)
    - back-splice and canonical splice-motif {[GT-AT],[GC-AG],[AT-AC]}, calculate strength using MaxEnt
    - back-splice and non-canonical splice-motif, find the maximal score using MaxEnt and the corresponding back-splice site
    - try to rescue circles using ">2-part-same-chromosome"
3. Check if both splice sites of the cirRNAs are on known exon-borders **_if_** an annotation is provided; otherwise all are CBR. report:
    - both known exon-border (termed : **_MEA_** .  **_M_** atch with **_E_** xisting **_A_** nnotation, which is somewhat more reliable)
    - at least one splice site sits on unknown exon-border (termed : **_CBR_** .   **_C_** hris **_B_** erge **_R_** escued)
4. Build gtf and pseudo-transcript for results from step3
5. Map all unmapped-reads to pseudo-transcipts and estimate the abundance of circRNAs
6. Generate bed track for visulization
7. Optional search of trans-splicing and fusion-circRNA events


# Before running ACFS, a few pre-process
0. Map the RNA-Seq reads to genome and transcriptom, and extract the unmapped reads. This is **_recommended_** as those mapped reads will NOT span the back-splice sites, and mapping them is a waste of time. 

1. Change fasta/fastq header format
    This is **_IMPORTANT_** ! ACFS expects a special header format so that multiple samples can be processed in one run. Do change the default header such as ">HWUSI-EAS100R:6:73:941:1973" into ">Truseq_sample1_HWUSI-EAS100R:6:73:941:1973", where the "sample1" is the name of your choice describing the sample. Do the conversion as:
    ```
    perl change_fastq_header.pl SRR650317_1.fasta SRR650317_1.fa Truseq_SRR650317left
    perl change_fastq_header.pl SRR650317_2.fasta SRR650317_2.fa Truseq_SRR650317right
    ```
    Make sure there is **_No underline_** within the sample name. e.g. ">Truseq_ctrl.1_Default_header" and ">Truseq_AGO2KO.1_Default_header" are OK; ">Truseq_ctrl_1_Default_header" is BAD because ACF will register "ctrl" as the sample name instead of "ctrl_1".

2. Merge sequences from multiple fasta/fastq files into one fasta file, which save time for mapping.
    ```
    perl Truseq_merge_unique_fa.pl UNMAP newid SRR650317_1.fa SRR650317_2.fa
    ```
    Alternatively, if there are many files to merge, generate a file (named filelist for example) contains the full-path of each file
    ```
    perl Truseq_merge_unique_filelist.pl UNMAP newid filelist
    ```
    UNMAP is the collasped fasta file which will be processed by ACFS, and UNMAP_expr contains the readcount per sequence in all the samples.
    
3. Build BWA index, using verion 0.73a (currently not support for other versions as the output format changes between versions)
    ```
    /bin/bwa073a/bwa index /data/iGenome/human/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa
    ```
    
4. Prepare for annotation (recommended)
    Download the gtf file from iGenome package    or    download ensembl gtf here : ftp://ftp.ensembl.org/pub/current_gtf/
    Then run:
    ```
    perl get_split_exon_border_biotype_genename.pl /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/genes.gtf /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf
    ```
    The first argument is the input gtf file, the second argument is the output

5. Strandedness is assumed as from the Truseq Stranded mRNA-Seq, so the reads are reverse-complementary to mRNAs.
    - If the reads are actually sense (the same 5'->3' direction as mRNA), please reverse-complement all reads.
    - If the reads are actually stransless or pair-ended, please run in parallel original reads and reverse-complemented reads.
    - No pair-end information is used, as the exact junction-site must be supported by a single read. Seeing is believing.
    

# Parameters
There are nine mandatory parameters to run ACFS in a basic mode, searching for fusion-circRNAs need to be enabled. Modify the config file "SPEC_example.txt" accordingly.

Mandatory paramters:

| Parameter | value | Note | 
| --------- | ----- | ---- | 
| BWA_folder | /home/bin/bwa037a/ | path of the folder of bwa | 
| BWA_genome_Index | /data/.../BWAIndex/genome.fa | full path to the index files | 
| BWA_genome_folder | /data/.../Chromosomes/ | full path to the folder containing **_separeted_** chromosome files | 
| ACF_folder | /home/bin/ACFS/ | path of the folder of ACFS | 
| CBR_folder | /home/bin/ACFS/CB_splice/ | path of the folder for MaxEnt | 
| Agtf | /data/.../Homo_sapiens.GRCh37.71_split_exon.gtf | processed annotation, see the previous section point-4 | 
| UNMAP | UNMAP | the collapsed fasta file | 
| UNMAP_expr | UNMAP_expr | the expression of the collasped reads | 
| Seq_len | 150 | length of sequencing reads | 

Optional parameters, the values in below are set as default:

| Parameter | value | Note |
| --------- | ----- | ---- |
| Thread | 16 | number of threads used in bwa |
| BWA_seed_length | 16 | bwa seed length  |
| BWA_min_score | 20 | bwa min score to triger report |
| minJump | 100 | the minimum distance of a back-splice. The smaller, the more likely you can find circles from short exons |
| maxJump | 1000000 | the maximum distance of a back-splice. The larger, the more likely you can find circles from long genes |
| minSplicingScore | 10 | the minimum score for the sum of splicing strength at both splice site, 10 corresponds to 97% of all human/mouse splice site pairs |
| minSampleCnt | 1 | the minimum number of samples that detect any given circle |
| minReadCnt | 1 | the minimum number of reads (from all samples) that detect any given circle |
| minMappingQuality | 20 | the minimum mapping quality of any given sequence |
| minSpanJunc | 6 | the minimum number of bases reach beyond the back-splice-site. The larger the less likely of false-positive |
| Coverage | 0.9 | the minimum percentage of any given read is aligned. The larger the better |
| ErrorRate | 0.05 | the maximum error rate for re-alignment. The smaller the better |
| Strandness | - | the strand information of sequencing, must be one of the {+, -, no}  |
| pre_defined_circle_bed | no | pre-defined circle annotation in bed6 or bed12 format (to increase sensitivity for lowly expressed circRNAs please include bed12 files of annotated one, e.g. merge the bed files from GSE61991) |
| Search_trans_splicing | no | set to "yes" to seach for trans-splicing reads |
| blat_search | yes | use blat to discard false positives results from gene duplication, turn off by "no" |
| blat_path | blat | full path to the executable, such as "/usr/bin/blat/blat", it is ignored if the blat_search option is set to no |
| trans_splicing_coverage | 0.9 | same as minMappingQuality |
| trans_splicing_minMappingQuality | 0 | the minimum mapping quality of any given sequence |
| trans_splicing_minSplicingScore | 10 | same as minSplicingScore |
| trans_splicing_maxSpan | 1000000 | the maximum distance to complete a fusion circRNA |


# run ACFS
1. make Pipeline Bash file
```
perl ACF_MAKE.pl SPEC_example.txt BASH_example.sh
```

2. find circles
```
bash BASH_example.sh
```


# Results
1. circRNAs are stored in two files, which can be visualzed using UCSC Genome Browser:
    - circle_candidates_MEA.bed12
    - circle_candidates_CBR.bed12  
    The higher the value in 5th column, the more "likely" that circle is true.  
    The name in 4th column can be seperated by "|" into four segments: circle-ID, sum-of-Splicing-Strength, number-of-supporting-samples, number-of-supporting-reads.  
    The sequence for the "middle exons" in almost every circle is hypothetically filled in using the annotations provided, true sequences could be determined by integrating mRNA-Seq data and validated by inward and outward PCRs.

2. The expression per sample is stored here:
    - circle_candidates_MEA.expr
    - circle_candidates_CBR.expr  
    The second column denote the name of the gene from which this circRNA is derived.

3. Fusion-circRNAs, if enabled:  
    - fusion_circRNAs  
    The Junctional sequences are stored in "unmap.trans.splicing.tsloci.fa".



# Tutorial
1. pre-processing for sequencing reads

    1A. **_for single-end RNA-Seq only_**. Note this sample is from mouse, therefore mouse annotation should be used. (Feasible but not a good idea in practice, since you don't want to map all reads again. Use unmapped reads instead.)
    ```
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX852%2FSRX852583/SRR1772422/SRR1772422.sra  
    fastq-dump.2 --fasta 0 --split-files SRR1772422.sra  
    perl change_fastq_header.pl SRR1772422.fasta SRR1772422.fa Truseq_HippoSyn  
    perl Truseq_merge_unique_fa.pl UNMAP newid SRR1772422.fa
    ```  
    1B. **_for paired-end RNA-Seq only_**. Note this sample is from human, therefore human annotation should be used. (Feasible but not a good idea in practice, since you don't want to map all reads again. Use unmapped reads instead.)
    ```
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX218%2FSRX218203/SRR650317/SRR650317.sra  
    fastq-dump.2 --fasta 0 --split-files SRR650317.sra  
    perl change_fastq_header.pl SRR650317_1.fasta SRR650317_1.fa Truseq_SRR650317left  
    perl change_fastq_header.pl SRR650317_2.fasta SRR650317_2.fa Truseq_SRR650317right  
    perl Truseq_merge_unique_fa.pl UNMAP newid SRR650317_1.fa SRR650317_2.fa  
    perl -ne 'chomp; if(m/^>/){s/^>//; print ">rc".$_,"\n";}else{my $t=$_; $t=~tr/[ATCG]/[TAGC]/; my $r=scalar reverse $t; print $r,"\n";}' UNMAP1 > UNMAP1.rc  
    cat UNMAP1 UNMAP1.rc > UNMAP  
    perl -e 'open IN,"UNMAP_expr1"; my $line=<IN>; print $line; while(<IN>){chomp; my @a=split("\t",$_); print join("\t",@a),"\n"; $a[0]="rc".$a[0]; print join("\t",@a),"\n";}' > UNMAP_expr
    ```
    1C. **_for single-end and paired-end RNA-Seq unmapped reads_**.  
    First align raw RNA-Seq reads to some reference with tools (such as Bowtie/BWA/STAR), obtain a SAM file containing all unmapped reads. For example if you choose Tophat, please convert <unmapped.bam> to <unmapped.sam>. Then run the following command:
    ```
    perl convert_unmapped_SAM_to_fa_for_acfs.pl unmapped.sam newName newNAME.fa
    ```  
    where <newName> is the new sample ID that is used to change the fasta/fastq header to ACF style. It could be "Ctrl12" or "TreatA.rep2", your choice. Just remember, there should be NO underscore as "_" in it.  
    <newNAME.fa> is the name of the converted file. Name it whatever you like.
    ```
    perl Truseq_merge_unique_fa.pl UNMAP newid newNAME.fa
    ```
2. prepare for annotation, if you haven't done so

    ```
    perl get_split_exon_border_biotype_genename.pl /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/genes.gtf /data/iGenome/human/Ensembl/GRCh37/Annotation/Genes/Homo_sapiens.GRCh37.71_split_exon.gtf
    ```
3. modify the config file SPEC_example.txt accordingly.
4. generate ACFS pipeline

    ```
    perl ACF_MAKE.pl SPEC_example.txt BASH_example.sh
    ```
5. run ACFS

    ```
    nohup bash BASH_example.sh &
    ```
6. take a look at the results when finished :

    - circle_candidates_MEA.bed12      # circRNAs back-splice at annotated boundarys of exon(s)
    - circle_candidates_CBR.bed12      # circRNAs back-splice at un-annotated boundarys of exon(s)  
    And the expression(readcounts) table for circRNAs, with circRNAs in rows and samples in columns
    - circle_candidates_MEA.expr       
    - circle_candidates_CBR.expr  
    And the potential fusion circRNAs  
    - fusion_circRNAs



# A few useful scripts for simulation  
To see the usage, simply run the perl scripts with no arguments.  

1. simulating SE reads from linear transcripts
    ```
    simulate_SE_reads_from_linear.pl
    ```

2. simulating PE reads from linear transcripts
    ```
    simulate_PE_reads_from_linear.pl
    ```

3. simulating circRNAs
    ```
    simulate_gtf_for_circRNA.pl
    get_split_exon_border_biotype_genename.pl
    get_seq_from_agtf.pl
    ```

4. simulating SE reads from circRNAs
    ```
    simulate_SE_reads_from_circRNA.pl
    ```

5. simulating PE reads from circRNAs
    ```
    simulate_PE_reads_from_circRNA.pl
    ```

6. simulating fusion-circRNAs
    ```
    simulate_gtf_for_fusion_circRNA.pl
    get_split_exon_border_biotype_genename.pl
    get_seq_from_agtf.pl
    ```

7. simulating SE reads from fusion-circRNAs
    ```
    simulate_reads_for_fusion_circRNA.pl
    ```



    
# Contact
This pipeline is developed and is maintained by Arthur Xintian You: arthur.yxt@gmail.com. I will do my best to respond in a timely manner.

# Cite
Nat Neurosci 2015 Apr;18(4):603-10 [PMID: 25714049]




