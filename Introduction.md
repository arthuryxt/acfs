# How to Find Your Circles #
> A pipeline to identify and quantify circRNA using single-ended stranded RNA sequencing data

# Introduction #
> Circular RNAs (cRNAs) are an old (first identified in 1993) class of cellular RNAs that possess neither 5' nor 3' ends. They are, as a group, lowly abundant and stable. As the majority of them arise from protein coding genes, presumably via "head-to-tail" splicing, they are thought to be biological relevant and one cRNA has been shown to be functionally important to brain. Yet many questions still lingers, like, How many cRNAs are out there? What's there abundance? How are they produced? Do they function, as a group, actively in the regulatory processes? The pipeline is designed to use high-throughput sequencing data to answer the first two questions.

> Also sync at https://github.com/arthuryxt/acfs

# Usage #
> Download the package

> Prepare according to "ACF\_readme.txt" : 1) fasta and expression files for unmapped reads; 2) non-overlapping exon-level gtf; 3) BWA genome index

> Generate parameter file (or modify "SPEC\_example.txt") and run bash like "bash BASH\_example.sh"

> Upload bed12 files to UCSC Genome Browser for visualisation