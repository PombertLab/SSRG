# SSRG - A Simple Pipeline to Assess Genetic Diversity Between Bacterial Genomes

This simple pipeline can:
1)	Download genomes automatically from NCBI
2)	Calculate pairwise SNPs between FASTQ sequences and reference genomes using standard read mapping approaches.
3)	Run Mash (https://github.com/marbl/Mash; Ondov et al. 2016. DOI: 10.1186/s13059-016-0997-x) and plot the estimated genetic distances as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

Automatic downloads:
- queryNCBI.pl downloads genomes from the NCBI Genomes FTP site. The Tab/CSV-delimited input lists can be generated from NCBI's Genome Assembly and Annotation reports (e.g. http://www.ncbi.nlm.nih.gov/genome/genomes/186?#).

SSRG (read mapping/variant calling):
-	SSRG.pl (Synthetic Short Read Generator) generates FASTQ datasets from fasta files at user specified read lengths and desired sequencing depth. Note that this tool is valid only for haploid genomes. This tool is especially useful to compare genomes in databases for which sequencing reads are either unavailable or originating from different sequencing technologies (e.g. illumina, Nanopore and PacBio)
-	get_SNPs.pl maps FASTQ files provided against references genomes using BWA, Bowtie2 or HISAT2, as specified by the user. SNPs and indels (optional) are then calculated with Samtools + VarScan2, BCFtools, or FreeBayes. Both synthetic reads generated by SSRG.pl and/or regular FASTQ files using the Sanger Q33 scoring scheme can be used as input.
-	bam2fastq.pl extracts sequencing reads in FASTQ format from BAM alignment files
-	SNPs reported in the VCF files can be summarized with sort_stats.pl and/or count_SNPs.pl.
-	Synonymous and non-synonymous mutations can be summarized with synonymy.pl ## Valid only for intron-less genomes 

MASH (genetic distance estimations):
-	run_Mash.pl runs the genetic distance estimation tool from Ondov et al. 2016. DOI: 10.1186/s13059-016-0997-x.
-	MashToDistanceMatrix.pl converts the output of Mash to distance matrices.
-	MashR_plotter.pl plots Mash-derived distance matrices as heatmaps, phylogenetic trees, or clusters.

Requirements:
- Unix/Linux or MacOS X
- Perl 5, R, and Java (for VarScan)
- Samtools version 1.3.1+ - http://www.htslib.org/ (Li et al. 2009. DOI: 10.1093/bioinformatics/btp352)
- Mash - https://github.com/marbl/Mash (Ondov et al. 2016. DOI: 10.1186/s13059-016-0997-x)

Compatible read aligners
- BWA version 0.7.12 - http://bio-bwa.sourceforge.net/ (Li et al. 2010. DOI: 10.1093/bioinformatics/btp698)
- Bowtie2 version 2.2.9 or above - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml (Langmead et al. 2012. DOI: 10.1038/nmeth.1923)
- HISAT2 version 2.0 or above - https://ccb.jhu.edu/software/hisat2/index.shtml (Kim et al. 2015. DOI: 10.1038/nmeth.3317)
- Minimap2 version 2.14 or above - https://github.com/lh3/minimap2 (Li, H. 2018. DOI: 10.1093/bioinformatics/bty191)
- NGMLR version 0.2.7 or above (partial support) - https://github.com/philres/ngmlr (Sedlazeck et al. 2018. DOI: 10.1038/s41592-018-0001-7)

Compatible variant callers
- VarScan2 version 2.4.2 or above - http://dkoboldt.github.io/varscan/ (Koboldt et al. 2012. DOI: 10.1101/gr.129684.111)
- BCFtools version 1.3.1 or above - http://samtools.github.io/bcftools/
- FreeBayes - https://github.com/ekg/freebayes (Garrison et al. 2012. arXiv preprint arXiv:1207.3907 [q-bio.GN])
