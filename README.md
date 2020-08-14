# SSRG - A Simple Pipeline to Assess Genetic Diversity Between Bacterial Genomes
The SSRG pipeline was created as a simple, focused tool to investigate SNPs between prokaryote genomes. The pipeline uses the common SNP-calling approach of read-mapping against references, standardizes experimental conditions for more accurate SNP comparisons, and integrates ubiquitous methodologies for both analysis and visualization. The unique feature of the SSRG pipeline resides in the creation of synthetic short reads from complete or draft genomes, which can then be fed to the read mapping/variant calling tools. Note that this approach works only for haploid genomes. Alternatively, users can also select any compatible FASTQ datasets to use with the pipeline.

### Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
* [Workflows](#workflows)
* [Installation](#installation)
* [References](#references)

### Introduction
Assessing the genetic diversity between genomes often involves the calculation of single nucleotide polymorphisms (SNPs) and insertions/deletions (indels). This is usually done by mapping short accurate sequencing reads from one or more species against a reference genome, from which variants are called. This approach works well when short read data from published genomes are available in public repositories, which is not always the case, especially now that bacterial genome sequencing is shifting towards the use of long read technologies. While genomes and/or long reads can be aligned against each other, the results are often suboptimal when the investigated chromosomes are highly reorganized, which can cause the mapping to fail. A simple solution to this problem is to deconstruct the genomes or long reads into shorter fragments, a shotgun approach, and to use these smaller synthetic reads as input for mapping.

***Deconstructing genomes into synthetic reads has the following advantages:***
-	This approach allows the comparisons of genomes for which sequencing read datasets are not available in public repositories.
-	This approach helps standardize datasets by providing reads with the exact same parameters. For example, genomes generated from [Illumina](https://www.illumina.com/), [PacBio](https://www.pacb.com/) and/or [Oxford Nanopore](https://nanoporetech.com/) data can now be compared without fuss.
-	Because bases from complete or draft genomes have been queried multiple times by the sequencing depth, the underlying confidence in the base being called is thus higher than from a single sequencing read. This in turn leads to fewer false positives caused by sequencing errors.

***The SSRG pipeline currently can:***
1)	Download genomes automatically from NCBI using a CSV/Tab-delimited list of desired operational taxonomic units (OTU)
2)	Calculate pairwise SNPs between FASTQ sequences and reference genomes using standard read mapping approaches.
3)	Run [Mash](https://github.com/marbl/Mash) (Ondov et al. 2016. DOI: 10.1186/s13059-016-0997-x) and plot the estimated genetic distances as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

### Requirements
- Unix/Linux or MacOS X
- Perl 5, R, and Java (for VarScan)
- [Samtools](http://www.htslib.org/) 1.3.1+
- [Mash](https://github.com/marbl/Mash)

###### Compatible read aligners
- [BWA](http://bio-bwa.sourceforge.net/) 0.7.12+
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.2.9+
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) 2.0+
- [Minimap2](https://github.com/lh3/minimap2) 2.14+
- [NGMLR](https://github.com/philres/ngmlr) 0.2.7+ (partial support)

###### Compatible variant callers
- [VarScan2](http://dkoboldt.github.io/varscan/) 2.4.2+
- [BCFtools](http://samtools.github.io/bcftools/) 1.3.1+
- [FreeBayes](https://github.com/ekg/freebayes)

### Workflows
The SSRG pipeline features two independent workflows:
1.	Read-mapping/variant calling
2.	Genetic distances estimation

Users interested in point mutations should use the read-mapping/variant calling workflow. Users only interested in genetic distances should use the genetic distances estimation workflow. This workflow is based on [Mash](https://github.com/marbl/Mash), an excellent tool developed by [Ondov *et al.*](https://pubmed.ncbi.nlm.nih.gov/27323842/). The Mash workflow does not identify point mutations.

<p align="center"><img src="https://github.com/PombertLab/SNPs/blob/master/Manual/Workflow.png" alt="Workflow" width="1000"></p>
**FIGURE 1 - OVERVIEW OF THE SSRG PIPELINE**  ***I.*** Genomes can be downloaded automatically from NCBI using provided scripts and custom or NCBI-generated lists. ***II.*** SSRG.pl generates FASTQ datasets from FASTA files at user-specified read lengths and desired sequencing depth. Note that this approach should be used only for haploid genomes. SSRG.pl is especially useful to compare genomes in databases for which sequencing reads are unavailable. ***III.*** get_SNPs.pl maps FASTQ files against references genomes using BWA [2], Bowtie2 [3], HISAT2 [4], Minimap2 [5] or NGMLR [6] as specified by the user. SNPs and indels (optional) are then calculated with Samtools [7] + VarScan2 [8], BCFtools [9], or FreeBayes [10]. ***IV.*** sort_stats.pl generates a tab-delimited table of SNP metrics. ***V.*** run_Mash.pl can estimate genetic distances using the MinHash Reduction technique, as implemented in Mash [3]. ***VI.*** MashToDistanceCSV.pl converts the output of Mash to distance matrices. ***VII.*** MashR_plotter.pl can A) clusters operational taxonomic units (OTUs) according to their estimated genetic distances, using R and either MDS [11] or t-SNE [12,13] algorithms, B) plot Neighbor-joining or UPGMA trees from Mash distances, C) generate clustered heatmaps from these distances.

### Installation
###### On Fedora/Red Hat
```bash
sudo dnf install \
  perl \
  R \
  bwa \
  boost \
  boost-devel \
  zlib \
  zlib-devel \
  gsl \
  gsl-devel \
  autoconf \
  automake \
  java-1.?.?-openjdk \
  java-1.?.?-openjdk-devel
```

###### Downloading/installing from GitHub
```bash
git clone --recursive https://github.com/PombertLab/SNPs.git
chmod a+x SNPs/*.pl; chmod a+x SNPs/*/*.pl
```
Install the scripts in your $PATH (e.g. by adding to your ~/.bash_profile). To install for all users, you can create a shell script in /etc/profile.d/ on most Linux systems. For example:
```bash
sudo export PATH="$PATH:/path/to/SNPs" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/SSRG/" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/MASH/" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/Tools/NCBI/" >> /etc/profile.d/SSRG.sh;
```
In the above, replace **/path/to/SNPs** by your installation directory.

###### Installing R packages dependencies
To install R packages for all users, type the following:
```bash
sudo R
```
```R
install.packages c("Rcpp", "gplots", "ggplot2", "ggfortify", "RColorBrewer", "plotly", "ape", "Rtsne")
quit()
```

###### Optional
If desired, update the VarScan default location in the corresponding line in [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl):
```perl
my $varjar = '/opt/varscan/VarScan.v2.4.3.jar';
```
This can also be changed from the command line with the -var option

If dependencies are not installed in your $PATH, you can alternatively insert the installation directories at the top of  [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) to reflect your settings:
```perl
## User-defined environment variables. ## Leave these blank (e.g. $samtools = '';) if programs are set in $PATH.
## Alternatively, you can insert here the installation directories (e.g. $samtools = '/opt/samtools-1.3.1/bin/') to reflect your settings.
my $samtools = '';		## Path to samtools 1.3.1+ - http://www.htslib.org/
my $bcftools = '';		## Path to bcftools 1.3.1+ - http://www.htslib.org/
my $bwa = '';			## Path to BWA - http://bio-bwa.sourceforge.net/
my $bowtie2 = '';		## Path to Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
my $minimap2 = '';		## Path to Minimap2 - https://github.com/lh3/minimap2
my $ngmlr = '';			## Path to ngmlr https://github.com/philres/ngmlr
my $hisat2 = '';		## Path to HISAT2 - https://ccb.jhu.edu/software/hisat2/index.shtml
my $freebayes = '';		## Path to FreeBayes -  https://github.com/ekg/freebayes
```


### References
1.	[Ondov BD, Treangen TJ, Mallonee AB, Bergman NH, Koren S, Phillippy AM. Fast genome and metagenome distance estimation using MinHash](https://pubmed.ncbi.nlm.nih.gov/27323842/) 2016. Genome Biol. 17, 132. PMID: 27323842 PMCID: PMC4915045 DOI: 10.1186/s13059-016-0997-x
2.	[Li H, Durbin R. Fast and accurate long-read alignment with Burrows-Wheeler transform](https://pubmed.ncbi.nlm.nih.gov/20080505/) 2010. Bioinformatics 26, 589–595. PMID: 20080505 PMCID: PMC2828108 DOI: 10.1093/bioinformatics/btp698
3.	[Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2](https://pubmed.ncbi.nlm.nih.gov/22388286/) 2012. Nat Methods 9, 357–359. PMID: 22388286 PMCID: PMC3322381 DOI: 10.1038/nmeth.1923
4.	[Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements](https://pubmed.ncbi.nlm.nih.gov/25751142/) 2015. Nat. Methods 12, 357–360. PMID: 25751142 PMCID: PMC4655817 DOI: 10.1038/nmeth.3317
5.	[Li H. Minimap2: pairwise alignment for nucleotide sequences](https://pubmed.ncbi.nlm.nih.gov/29750242/) 2018. Bioinformatics 34, 3094–3100. PMID: 29750242 PMCID: PMC6137996 DOI: 10.1093/bioinformatics/bty191
6.	[Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A,  Schatz MC. Accurate detection of complex structural variations using single-molecule sequencing](https://pubmed.ncbi.nlm.nih.gov/29713083/) 2018. Nat. Methods 15, 461–468. PMID: 29713083 PMCID: PMC5990442 DOI: 10.1038/s41592-018-0001-7
7.	[Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. The Sequence Alignment/Map format and SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/) 2009. Bioinformatics 25, 2078–2079. PMID: 19505943 PMCID: PMC2723002 DOI: 10.1093/bioinformatics/btp352
8.	[Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan, MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing](https://pubmed.ncbi.nlm.nih.gov/22300766/) 2012. Genome Res. 22, 568–576. PMID: 22300766 PMCID: PMC3290792 DOI: 10.1101/gr.129684.111
9.	[Narasimhan V, Danecek P, Scally A, Xue Y, Tyler-Smith C, Durbin R. BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data](https://pubmed.ncbi.nlm.nih.gov/26826718/) 2016. Bioinformatics 32, 1749–1751. PMID: 26826718 PMCID: PMC4892413 DOI: 10.1093/bioinformatics/btw044
10.	[Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing](https://arxiv.org/abs/1207.3907) 2012. arXiv Prepr. arXiv1207.3907, 9.
11.	[R Core Team. R: A language and environment for statistical computing](https://www.r-project.org/). 2020.
12.	[van der Maaten LJP. Barnes-Hut-SNE](https://arxiv.org/abs/1301.3342) 2013. arxiv.org, 1301.3342v2.
13.	[van der Maaten LJP,  Hinton GE. Visualizing High-Dimensional Data Using t-SNE](https://jmlr.csail.mit.edu/papers/volume9/vandermaaten08a/vandermaaten08a.pdf) 2008. J. Mach. Learn. Res. 9, 2579–2605.



