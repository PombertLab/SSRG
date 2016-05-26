# SNPs
A Simple Pipeline to Assess Genetic Diversity Between Bacterial Genomes

This simple tool calculates pairwise SNPs for every possible combination in the genomes provided as input.

-	SSRG.pl generates FASTQ datasets from fasta files at user specified read lengths, with coverages of 50X or 100X. Note that this works only for haploid genomes. This tool is especially useful to compare genomes in databases for which sequencing reads are unavailable.
-	get_SNPs.pl maps FASTQ files provided against references genomes using BWA, Bowtie2 or HISAT2, as specified by the user. SNPs and indels (optional) are then calculated with VarScan2.
-	SNPs reported in the VCF files can be summarized with count_SNPs.pl.
