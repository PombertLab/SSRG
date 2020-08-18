<p align="center"><img src="https://github.com/PombertLab/SNPs/blob/master/Manual/logo.png" alt="SSRG - A simple pipeline to assess genetic diversity" width="900"></p>

The [SSRG pipeline](https://github.com/PombertLab/SNPs) was created as a simple and focused tool to investigate genetic diversity between genomes. The pipeline features two independent workflows:
1.	Read-mapping/variant calling
2.	Genetic distance estimation

People interested in ***point mutations*** should use the read-mapping/variant calling workflow. The unique feature of the SSRG pipeline resides in the creation with [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl) of synthetic short reads  from complete or draft genomes, which can then be fed to the read mapping/variant calling tools. Note that this approach works only for haploid genomes. Alternatively, users can select any FASTQ dataset to use with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl).

People interested in ***genetic distances*** should use the genetic distance estimation workflow. This workflow is based on [Mash](https://github.com/marbl/Mash), an excellent tool developed by [Ondov *et al.*](https://pubmed.ncbi.nlm.nih.gov/27323842/). Genetic distances can be plotted as as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

## Table of contents
* [Introduction](#introduction)
* [Installation](#installation)
  * [Requirements](#requirements)
  * [Compatible read aligners](#compatible-read-aligners)
  * [Compatible variant callers](#compatible-variant-callers)
  * [Downloading from GitHub](#downloading-from-github)
  * [Installing dependencies](#installing-dependencies)
  * [Testing the installation](#testing-the-installation)
* [Howto (Case examples)](#howto)
  * [Read mapping + variant calling](#read-mapping-and-variant-calling)
    * [Synthetic reads](#synthetic-reads)
    * [Read mapping](#read-mapping)
    * [Variant calling](#variant-calling)
  * [Genetic distance estimation with Mash](#genetic-distance-estimation-with-Mash)
* [References](#references)

## Introduction
##### Why another read mapping pipeline?
Assessing the genetic diversity between genomes often involves the calculation of ***[single nucleotide polymorphisms](https://ghr.nlm.nih.gov/primer/genomicresearch/snp) (SNPs)*** and ***[insertions/deletions](https://ghr.nlm.nih.gov/primer/mutationsanddisorders/possiblemutations) (indels)***. This is usually done by mapping short accurate sequencing reads from one or more species against a reference genome, from which variants are called. This approach works well when short read data from published genomes are available in public repositories such as [NCBI's SRA](https://www.ncbi.nlm.nih.gov/sra) but that is not always the case, especially now that genome sequencing is shifting towards the use of long read technologies. While genomes and/or long reads can be aligned against each other, the results are often suboptimal when the investigated chromosomes are highly reorganized, which can cause the mapping to fail. A simple solution to this problem is to deconstruct the genomes or long reads randomly into shorter fragments —much like DNA fragmentation protocols used in ***[whole genome sequencing](https://ghr.nlm.nih.gov/primer/testing/sequencing) (WGS)***— and to use these smaller synthetic reads as input for mapping. We have implemented this approach in [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl). Note that this approach is only valid for haploid genomes and should not be used in presence of heterozygosity.

##### Deconstructing genomes into synthetic reads has the following advantages:
-	It enables comparisons between genomes for which sequencing datasets are not available in public repositories.
-	It helps standardize datasets by providing reads with the exact same parameters. For example, genomes generated from [Illumina](https://www.illumina.com/), [PacBio](https://www.pacb.com/) and/or [Oxford Nanopore](https://nanoporetech.com/) data can now be compared without fuss.
-	Because bases from complete or draft genomes have been queried multiple times by the sequencing depth, the underlying confidence in the base being called is higher than from single sequencing reads (assuming of course that the corresponding loci in the assembled genomes are error-free). This in turn should lead to fewer problems caused by sequencing errors.

##### The SSRG pipeline currently can:
1)	Download genomes automatically from NCBI using a CSV/Tab-delimited list of desired operational taxonomic units (OTU)
2)	Calculate pairwise SNPs between FASTQ sequences and reference genomes using standard read mapping approaches.
3)	Run [Mash](https://github.com/marbl/Mash) and plot the estimated genetic distances as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

## Installation
##### Requirements
- Unix/Linux, MacOS X, or WSL2
- [Perl 5](https://www.perl.org/), [R](https://www.r-project.org/), and [Java](https://www.java.com/) (for VarScan2)
- [Samtools](http://www.htslib.org/) 1.3.1+
- [Mash](https://github.com/marbl/Mash) (for genetic distance estimations)

##### Compatible read aligners
- [BWA](http://bio-bwa.sourceforge.net/) 0.7.12+
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.2.9+
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) 2.0+
- [Minimap2](https://github.com/lh3/minimap2) 2.14+ (recommended)
- [NGMLR](https://github.com/philres/ngmlr) 0.2.7+ (partial support)

##### Compatible variant callers
- [VarScan2](http://dkoboldt.github.io/varscan/) 2.4.2+
- [BCFtools](http://samtools.github.io/bcftools/) 1.3.1+
- [FreeBayes](https://github.com/ekg/freebayes)

### Downloading from GitHub
```bash
git clone --recursive https://github.com/PombertLab/SNPs.git
chmod a+x `find ./SNPs/ -name *.pl`
```
For ease of use, the SSRG folders should be added to the $PATH variable. Bash users on Fedora/RedHat can type the following commands to add the SNPs folder and its subdirectories to their own .bash_profile:
```bash
cd ./SNPs; PWD=`pwd`;
echo 'PATH="$PATH:'"$PWD\"" >> ~/.bash_profile
echo 'PATH="$PATH:'"$PWD/SSRG\"" >> ~/.bash_profile
echo 'PATH="$PATH:'"$PWD/MASH\"" >> ~/.bash_profile
echo 'PATH="$PATH:'"$PWD/Tools/NCBI\"" >> ~/.bash_profile
echo 'export PATH' >> ~/.bash_profile
source ~/.bash_profile
```

### Installing dependencies
##### ***On Fedora/RedHat***
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
  libcurl-devel \
  openssl-devel \
  ncurses-devel \
  java-1.?.?-openjdk \
  java-1.?.?-openjdk-devel
```

##### ***R***

To install R dependencies, type:
```bash
R
```
Then in R, type:
```R
install.packages("ape")
install.packages("gplots")
install.packages("ggplot2")
install.packages("ggfortify")
install.packages("plotly")
install.packages("RColorBrewer")
install.packages("Rcpp")
install.packages("Rtsne")
quit()
```

##### ***Installing read mappers, variants callers and Mash***

Read mapping/variant calling with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) require [Samtools](http://www.htslib.org/) 1.3.1+. Read mappers and variant callers to be used with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) can be installed on a need-to basis. [Minimap2](https://github.com/lh3/minimap2) and [VarScan2](http://dkoboldt.github.io/varscan/) are recommended.

*Read mappers are available here:*

- [Minimap2](https://github.com/lh3/minimap2)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [NGMLR](https://github.com/philres/ngmlr)
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)
- [BWA](http://bio-bwa.sourceforge.net/). ## On Fedora, BWA can be installed from the DNF package manager:
```Bash
sudo dnf install bwa
```
*Variant callers are available here:*
- [BCFtools](http://samtools.github.io/bcftools/)
- [FreeBayes](https://github.com/ekg/freebayes)
- [VarScan2](http://dkoboldt.github.io/varscan/). ## The default location of the VarScan2 jar file can be updated in [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) by modifying the following line:
```perl
my $varjar = '/opt/varscan/VarScan.v2.4.4.jar';
```

To perform genetic distance estimations, [Mash](https://github.com/marbl/Mash) by [Ondov *et al.*](https://pubmed.ncbi.nlm.nih.gov/27323842/) is required. It is available [here](https://github.com/marbl/Mash/releases).

### Testing the installation
##### SSRG workflow (with Minimap2 and VarScan2)
1. Download a small CSV list of 3 *Streptococcus pneumoniae* genomes. ## This list was generated from the NCBI genome database (https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/176/)
```bash
wget https://raw.githubusercontent.com/PombertLab/SNPs/master/Misc/S_pneumoniae_3.csv
```
2. Download the corresponding FASTA and GenBank files from NCBI with [queryNCBI.pl](https://github.com/PombertLab/SNPs/blob/master/Tools/NCBI/queryNCBI.pl).
```bash
queryNCBI.pl -l S_pneumoniae_3.csv -fa -gb
```

3. Generate FASTQ datasets from the downloaded genomes with [SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl).
```bash
SSRG.pl -f *.fasta -r 250
```

4a. To test the read-mapping step with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) but skip variant calling, type:
```bash
get_SNPs.pl \
	-threads 2 \
	-mem 4 \
	-fa Streptococcus_pneumoniae_R6.fasta \
	-pe1 *R1.fastq \
	-pe2 *R2.fastq \
	-mapper minimap2 \
	-preset sr \
	-rmo \
	-bam
```
4b. To test the read-mapping and variant calling steps with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl), type:
```bash
## Running get_SNPs.pl with 8 threads, 16 Gb RAM (change according to your settings)
get_SNPs.pl \
	-threads 8 \
	-mem 16 \
	-fa Streptococcus_pneumoniae_R6.fasta \
	-pe1 *R1.fastq \
	-pe2 *R2.fastq \
	-mapper minimap2 \
	-preset sr \
	-bam \
	-caller varscan2 \
	-var ./VarScan.v2.4.4.jar
```

##### MASH workflow
1. Download a CSV list of 75 *Streptococcus pneumoniae* genomes. ## This list was generated from the NCBI genome database (https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/176/)
```bash
wget https://raw.githubusercontent.com/PombertLab/SNPs/master/Misc/S_pneumoniae_75.csv
```

2. Download the FASTA files automatically with [queryNCBI.pl](https://github.com/PombertLab/SNPs/blob/master/Tools/NCBI/queryNCBI.pl), then run Mash with [run_Mash.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/run_Mash.pl).
```bash
mkdir FASTA; cd FASTA/;
queryNCBI.pl -l ../S_pneumoniae_75.csv -fa
run_Mash.pl -f *.fasta -o ../S_pneumoniae_75.mash
cd ../
```
3. Convert the Mash output to a distance matrix with [MashToDistanceMatrix.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashToDistanceMatrix.pl).
```bash
MashToDistanceMatrix.pl  \
	-i S_pneumoniae_75.mash \
	-o S_pneumoniae_75 \
	-f tsv
```
4. Generate a quick neighbor-joining tree with [MashR_plotter.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashR_plotter.pl). The PDF generated should be similar to [S_pneumoniae_75_NJ_tree.pdf](https://github.com/PombertLab/SNPs/blob/master/Misc/S_pneumoniae_75_NJ_tree.pdf) from the Misc subdirectory.
```bash
MashR_plotter.pl \
	-i S_pneumoniae_75.tsv \
	-if tsv \
	-t tree \
	-newick S_pneumoniae_75.tre \
	-f pdf \
	-o S_pneumoniae_75_NJ_tree \
	-he 20
```
5. Generate a quick heatmap with [MashR_plotter.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashR_plotter.pl). The PDF generated should be similar to [S_pneumoniae_75_heatmap.pdf](https://github.com/PombertLab/SNPs/blob/master/Misc/S_pneumoniae_75_heatmap.pdf) from the Misc subdirectory.
```bash
MashR_plotter.pl \
	-i S_pneumoniae_75.tsv \
	-if tsv \
	-t heatmap \
	-f pdf \
	-o S_pneumoniae_75_heatmap \
	-colors white cyan magenta \
	-he 20 \
	-wd 20
```
6. Generate a quick t-SNE multidimensional reduction plot with [MashR_plotter.pl](https://github.com/PombertLab/SNPs/blob/master/MASH/MashR_plotter.pl). The PDF generated should be similar, ***but not identical***, to [S_pneumoniae_75_tSNE.pdf](https://github.com/PombertLab/SNPs/blob/master/Misc/S_pneumoniae_75_tSNE.pdf) from the Misc subdirectory. ## t-SNE graphs are generated using random seeds, which effect how the distances are represented in 2D.
```bash
MashR_plotter.pl \
	-i S_pneumoniae_75.tsv \
	-if tsv \
	-t cluster \
	-m tsne \
	-pe 15 \
	-cmode terrain \
	-f pdf \
	-o S_pneumoniae_75_tSNE \
	-he 20 \
	-wd 20 \
	-lb \
	-fs 25
```

## Howto
#### Read mapping and variant calling
##### Synthetic reads
[SSRG.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/SSRG.pl) can be used to generate FASTQ datasets from complete or draft genomes (ploidy = 1). For example, to generate FASTQ datasets (paired ends; 250 bp; 50x sequencing depth) from one or more genomes:
```Bash
SSRG.pl -f *.fasta -r 250 
```
Options for SSRG.pl are:
```
 -f (--fasta)		Fasta/multifasta file(s)
 -l (--list)		List of fasta file(s), one per line
 -r (--readsize)	Desired read size [default: 150]
 -t (--type)		Read type: single (SE) or paired ends (PE) [default: PE]
 -i (--insert)		PE insert size [default: 350]
 -s (--sdev)		PE insert size standard deviation (in percentage) [default: 10]
 -c (--coverage)	Desired sequencing depth [default: 50]
 -qs (--qscore)		Quality score associated with each base [default: 30]
 -q64          		Use the old Illumina Q64 FastQ format instead of the default Q33 Sanger/Illumina 1.8+ encoding
```


##### Read mapping
[get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) is the main script that handles read mapping and variant calling. The default read mapper in get_SNPs.pl is [Minimap2](https://github.com/lh3/minimap2), and can be changed by invoking the **--mapper** switch from the command line. Users can also change the default read mapper settings by modifying the following line in [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl):
```perl
my $mapper = 'minimap2';
```

[get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) can be used with single or paired ends datasets. BAM files generated during the alignments can be kept with the **--bam** command line switch (BAM files are discarded by default to save on disk space). The **--rmo** option (read mapping only) skips the variant calling process. When combined, the **--rmo** and **--bam** options generate BAM alignment files using the selected read mapper.

To generate BAM alignments with ***single end*** FASTQ datasets using the default read mapper, [Minimap2](https://github.com/lh3/minimap2), with its default sr (short reads) preset:
```bash
get_SNPs.pl -fa *.fasta -fq *.fastq -rmo -bam
```

To generate BAM alignments with ***paired end*** FASTQ datasets, this time using [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) as read mapper:
```Bash
get_SNPs.pl -fa *.fasta -pe1 *R1.fastq -pe2 *R2.fastq -mapper bowtie2 -rmo -bam
```

General and read mapping options for get_SNPs.pl are:
```
### Mapping options ###
-fa (--fasta)			Reference genome(s) in fasta file
-fq (--fastq)			Fastq reads (single ends) to be mapped against reference(s)
-pe1				Fastq reads #1 (paired ends) to be mapped against reference(s)
-pe2				Fastq reads #2 (paired ends) to be mapped against reference(s)
-mapper				Read mapping tool: bwa, bowtie2, minimap2, ngmlr or hisat2 [default: minimap2]
-threads			Number of processing threads [default: 16]
-mem				Max total memory for samtools (in Gb) [default: 16] ## mem/threads = memory per thread
-bam				Keeps BAM files generated
-sam				Keeps SAM files generated; SAM files can be quite large
-rmo (--read_mapping_only)	Do not perform variant calling; useful when only interested in bam/sam files and/or mapping stats
-ns (--no_stats)		Do not calculate stats; stats can take a while to compute for large eukaryote genomes

### Mapper-specific options ###
-X				BOWTIE2 - Maximum paired ends insert size [default: 750]
-preset				MINIMAP2 - Preset: sr, map-ont, map-pb or asm20 [default: sr]
-algo				BWA - Mapping algorithm:  bwasw, mem, samse [default: mem]
```

##### Variant calling
The default variant caller in get_SNPs.pl is [VarScan2](http://dkoboldt.github.io/varscan/). Other variant callers are implemented only partially. To perform variant calling with [get_SNPs.pl](https://github.com/PombertLab/SNPs/blob/master/SSRG/get_SNPs.pl) using default settings ([Minimap2](https://github.com/lh3/minimap2) + [VarScan2](http://dkoboldt.github.io/varscan/)) and single end datasets, simply type:
```bash
get_SNPs.pl -fa *.fasta -fq *.fastq
```

A more detailed command line using Minimap2, VarScan2, and paired end datasets should look like:
```bash
get_SNPs.pl \
   -threads 16 \
   -fa *.fasta \
   -pe1 *R1.fastq \
   -pe2 *R2.fastq \
   -mapper minimap2 \
   -caller varscan2 \
   -type both \
   -var ./VarScan.v2.4.3.jar
```

Variant calling options for get_SNPs.pl are:
```
### Variant calling options ###
-caller				[default: varscan2]	## Variant caller: varscan2, bcftools or freebayes
-type				[default: snp]		## snp, indel, or both
-ploidy				[default: 1]		## FreeBayes/BCFtools option; change ploidy (if needed)

### VarScan2 parameters ### see http://dkoboldt.github.io/varscan/using-varscan.html
-var				[default: /opt/varscan/VarScan.v2.4.3.jar]	## Which varscan jar file to use
-mc (--min-coverage)		[default: 15]		## Minimum read depth at a position to make a call
-mr (--min-reads2)		[default: 15]		## Minimum supporting reads at a position to call variants
-maq (--min-avg-qual)		[default: 28]		## Minimum base quality at a position to count a read
-mvf (--min-var-freq)		[default: 0.7]		## Minimum variant allele frequency threshold
-mhom (--min-freq-for-hom)	[default: 0.75]		## Minimum frequency to call homozygote
-pv (--p-value)			[default: 1e-02]	## P-value threshold for calling variants 
-sf (--strand-filter)		[default: 0]		## 0 or 1; 1 ignores variants with >90% support on one strand
```

#### Genetic distance estimation with Mash
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

## References
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



