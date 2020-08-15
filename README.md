<p align="center"><img src="https://github.com/PombertLab/SNPs/blob/master/Manual/logo.png" alt="SSRG - A simple pipeline to assess genetic diversity" width="900"></p>

The SSRG pipeline was created as a simple, focused tool to investigate genetic diversity between genomes. The pipeline features two independent workflows:
1.	Read-mapping/variant calling
2.	Genetic distance estimation

People interested in point mutations should use the read-mapping/variant calling workflow. The unique feature of the SSRG pipeline resides in the creation with [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl) of synthetic short reads  from complete or draft genomes, which can then be fed to the read mapping/variant calling tools. Note that this approach works only for haploid genomes. Alternatively, users can select any compatible FASTQ datasets to use with the pipeline.

People interested in genetic distances should use the genetic distance estimation workflow. This workflow is based on [Mash](https://github.com/marbl/Mash), an excellent tool developed by [Ondov *et al.*](https://pubmed.ncbi.nlm.nih.gov/27323842/). Genetic distances can be plotted as as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

### Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
  * [Requirements](#requirements)
  * [Compatible read aligners](#compatible-read-aligners)
  * [Compatible variant callers](#compatible-variant-callers)
  * [Downloading from GitHub](#downloading-from-github)
  * [Installing dependepencies](#installing-dependepencies)
* [Howto](#howto)
  * [Read mapping + variant calling](#read-mapping-and-variant-calling)
  * [Genetic distance estimation](#genetic-distance-estimation)
* [References](#references)

### Introduction
###### Why another read mapping pipeline?
Assessing the genetic diversity between genomes often involves the calculation of single nucleotide polymorphisms (SNPs) and insertions/deletions (indels). This is usually done by mapping short accurate sequencing reads from one or more species against a reference genome, from which variants are called. This approach works well when short read data from published genomes are available in public repositories but that is not always the case, especially now that genome sequencing is shifting towards the use of long read technologies. While genomes and/or long reads can be aligned against each other, the results are often suboptimal when the investigated chromosomes are highly reorganized, which can cause the mapping to fail. A simple solution to this problem is to deconstruct the genomes or long reads randomly into shorter fragments, much like DNA fragmentation protocols used in whole genome sequencing (WGS), and to use these smaller synthetic reads as input for mapping. We have implemented this approach in [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl). Note that this approach is only valid for haploid genomes and should not be used in presence of heterozygosity.

###### Deconstructing genomes into synthetic reads has the following advantages:
-	It enables comparisons between genomes for which sequencing datasets are not available in public repositories.
-	It helps standardize datasets by providing reads with the exact same parameters. For example, genomes generated from [Illumina](https://www.illumina.com/), [PacBio](https://www.pacb.com/) and/or [Oxford Nanopore](https://nanoporetech.com/) data can now be compared without fuss.
-	Because bases from complete or draft genomes have been queried multiple times by the sequencing depth, the underlying confidence in the base being called is thus higher than from a single sequencing read. This in turn should lead to fewer problems caused by sequencing errors.

###### The SSRG pipeline currently can:
1)	Download genomes automatically from NCBI using a CSV/Tab-delimited list of desired operational taxonomic units (OTU)
2)	Calculate pairwise SNPs between FASTQ sequences and reference genomes using standard read mapping approaches.
3)	Run [Mash](https://github.com/marbl/Mash) and plot the estimated genetic distances as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

### Installation
###### Requirements
- Unix/Linux, MacOS X, or WSL2
- Perl 5, R, and Java (for VarScan)
- [Samtools](http://www.htslib.org/) 1.3.1+
- [Mash](https://github.com/marbl/Mash) (for genetic distance estimations)

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

###### Downloading from GitHub
```bash
git clone --recursive https://github.com/PombertLab/SNPs.git
chmod a+x SNPs/*.pl; chmod a+x SNPs/*/*.pl
```
Install the scripts in your $PATH (e.g. by adding to your ~/.bash_profile).

To install for all users, a shell script can be created in /etc/profile.d/ on most Linux systems. For example:
```bash
sudo export PATH="$PATH:/path/to/SNPs" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/SSRG/" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/MASH/" >> /etc/profile.d/SSRG.sh; \
sudo export PATH="$PATH:/path/to/SNPs/Tools/NCBI/" >> /etc/profile.d/SSRG.sh;
```
In the above, replace **/path/to/SNPs** by your installation directory.

###### Installing dependepencies
On Fedora/RedHat:
```bash
sudo dnf install \
  perl \
  R \
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

To install R packages dependencies for all users, type the following:
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

### Howto
#### Read mapping and variant calling
Creating short read datasets with [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl) (paired ends; 250 bp; 50x sequencing depth) ## Optional
```Bash
SSRG.pl -f *.fasta -r 250 
```
Performing read mapping with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl)
- Read mapping only (using minimap2 on a paired end dataset with 16 threads)
```
get_SNPs.pl --threads 16 --fa *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --mapper minimap2 --rmo --bam
```
- Read mapping (minimap2; paired end dataset) + variant calling (varscan2; SNPs + indels)
```
get_SNPs.pl --threads 16 --fa *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --mapper minimap2 --caller varscan2 --type both --var ./VarScan.v2.4.3.jar
```

#### Genetic distance estimation
Runnning Mash with [run_Mash.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/run_Mash.pl) and converting output to distance matrices with [MashToDistanceMatrix.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashToDistanceMatrix.pl)
```
run_Mash.pl -f *.fasta -o Mash.txt
MashToDistanceMatrix.pl  -i Mash.txt -o Mash -f tsv
```
Plotting a quick Neighbor-joining tree with [MashR_plotter.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashR_plotter.pl)
```
MashR_plotter.pl -i Mash.tsv -if tsv -t tree -newick Mash.txt.tre
```
Plotting the same data using dimensionality reduction techniques
MashR_plotter.pl -t cluster -m tsne -i Mash.tsv -if tsv -o cluster_tsne --format pdf -fs 8 -lb -pe 10


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



