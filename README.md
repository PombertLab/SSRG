<p align="center"><img src="https://github.com/PombertLab/SSRG/blob/master/Images/logo.png" alt="SSRG - A simple pipeline to assess genetic diversity" width="900"></p>

The [SSRG pipeline](https://github.com/PombertLab/SSRG) was created as a simple and focused tool to investigate genetic diversity between genomes. The pipeline features two independent workflows:
1.	Read-mapping/variant calling
2.	Genetic distance estimation

People interested in ***point mutations*** should use the read-mapping/variant calling workflow. The unique feature of the SSRG pipeline resides in the creation with [SSRG.pl](https://github.com/PombertLab/SSRG/blob/master/SSRG.pl) of synthetic short reads  from complete or draft genomes, which can then be fed to the read mapping/variant calling tools. Note that this approach works only for haploid genomes. Alternatively, users can select any FASTQ dataset to use with [get_SNPs.pl](https://github.com/PombertLab/SSRG/blob/master/get_SNPs.pl).

People interested in ***genetic distances*** should use the genetic distance estimation workflow. This workflow is based on [Mash](https://github.com/marbl/Mash), an excellent tool developed by [Ondov *et al.*](https://pubmed.ncbi.nlm.nih.gov/27323842/), and focuses on data visualization. Genetic distances can be plotted as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

## Table of contents
* [Introduction](#Introduction)
* [Requirements](#Requirements)
  * [Downloading from GitHub](#downloading-from-github)
  * [Installing dependencies](#installing-dependencies)
* [Examples](#Examples)
  * [Creating NCBI genome lists](#Creating-NCBI-genome-lists)
  * [SSRG workflow](#SSRG-workflow)
  * [Mash workflow](#Mash-workflow)
* [Misc scripts](#Misc-scripts)
* [Funding and acknowledgments](#Funding-and-acknowledgments)
* [References](#references)

## Introduction
##### Why another read mapping pipeline?
Assessing the genetic diversity between genomes often involves the calculation of ***[single nucleotide polymorphisms](https://ghr.nlm.nih.gov/primer/genomicresearch/snp) (SNPs)*** and ***[insertions/deletions](https://ghr.nlm.nih.gov/primer/mutationsanddisorders/possiblemutations) (indels)***. This is usually done by mapping short accurate sequencing reads from one or more species against a reference genome, from which variants are called. This approach works well when short read data from published genomes are available in public repositories such as [NCBI's SRA](https://www.ncbi.nlm.nih.gov/sra) but that is not always the case, especially now that genome sequencing is shifting towards the use of long read technologies. While genomes and/or long reads can be aligned against each other, the results are often suboptimal when the investigated chromosomes are highly reorganized, which can cause the mapping to fail. A simple solution to this problem is to deconstruct the genomes or long reads randomly into shorter fragments —much like DNA fragmentation protocols used in ***[whole genome sequencing](https://ghr.nlm.nih.gov/primer/testing/sequencing) (WGS)***— and to use these smaller synthetic reads as input for mapping. We have implemented this approach in [SSRG.pl](https://github.com/PombertLab/SSRG/blob/master/SSRG.pl). Note that this approach is only valid for haploid genomes and should not be used in presence of heterozygosity.

##### Deconstructing genomes into synthetic reads has the following advantages:
-	It enables comparisons between genomes for which sequencing datasets are not available in public repositories.
-	It helps standardize datasets by providing reads with the exact same parameters. For example, genomes generated from [Illumina](https://www.illumina.com/), [PacBio](https://www.pacb.com/) and/or [Oxford Nanopore](https://nanoporetech.com/) data can now be compared without fuss.
-	Because bases from complete or draft genomes have been queried multiple times by the sequencing depth, the underlying confidence in the base being called is higher than from single sequencing reads (assuming of course that the corresponding loci in the assembled genomes are error-free). This in turn should lead to fewer problems caused by sequencing errors.

##### The SSRG pipeline currently can:
1)	Download genomes automatically from [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/) using CSV/Tab-delimited lists of operational taxonomic units (OTU)
2)	Calculate pairwise SNPs between FASTQ sequences and reference genomes using standard read mapping approaches.
3)	Run [Mash](https://github.com/marbl/Mash) and plot the estimated genetic distances as heatmaps, neighbor-joining trees, or clusters (using dimensionality reduction techniques).

## Requirements
- [Perl 5](https://www.perl.org/)
#### For read mapping
- [Samtools](http://www.htslib.org/) 1.3.1+
- [Minimap2](https://github.com/lh3/minimap2) 2.14+
- [VarScan2](http://dkoboldt.github.io/varscan/) 2.4.2+
- [Java](https://www.java.com/) (for VarScan2)
#### For genetic distance estimations with Mash
- [Mash](https://github.com/marbl/Mash)
- [R](https://www.r-project.org/)

#### Optional
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 2.2.9+
- [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml) 2.0+
- [NGMLR](https://github.com/philres/ngmlr) 0.2.7+ (partial support)
- [BCFtools](http://samtools.github.io/bcftools/) 1.3.1+
- [FreeBayes](https://github.com/ekg/freebayes)

#### Downloading from GitHub
```bash
git clone https://github.com/PombertLab/SSRG.git
cd SSRG/
export PATH=$PATH:$(pwd)/
export PATH=$PATH:$(pwd)/Tools/MASH
export PATH=$PATH:$(pwd)/Tools/NCBI
```

To download the latest [VarScan2](http://dkoboldt.github.io/varscan/) with curl, type:
```Bash
curl \
   -o VarScan.v2.4.4.jar \
   -L https://github.com/dkoboldt/varscan/blob/master/VarScan.v2.4.4.jar?raw=true
```

If desired, the default location of the [VarScan2](http://dkoboldt.github.io/varscan/) jar file can be updated in [get_SNPs.pl](https://github.com/PombertLab/SSRG/blob/master/get_SNPs.pl) by modifying the following line:
```perl
my $varjar = 'VarScan.v2.4.4.jar';
```

#### Installing dependencies
##### ***On Fedora/RedHat***
On a Fedora distribution, the following packages can be installed with the DNF package manager to facililate the installation/compilation of read mappers and other tools from the source code.
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
##### ***R packages***
To install R data visualization packages, type:
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


## Examples
#### Creating NCBI genome lists
Comma-separated (.csv) lists of genomes for organisms of interest can be downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/genome/browse#!/overview/).
<p align="center"><img src="https://github.com/PombertLab/SSRG/blob/master/Images/NCBI_g1.png" alt="The NCBI Genome Information by Organism page" width="900"></p>

This can be done by entering a query in the search box, then by clicking on the lineage link (*e.g.* ***Prokaryotes***) next to ***Overview***. This links redirects to the corresponding list of genomes in NCBI. 
<p align="center"><img src="https://github.com/PombertLab/SSRG/blob/master/Images/NCBI_g2.png" alt="Searching the NCBI Genome Information page" width="900"></p>

Filters can be applied, then the list of corresponding genomes downloaded. Note that while the tooltip displays TSV Download, the file downloaded will be in CSV format.
<p align="center"><img src="https://github.com/PombertLab/SSRG/blob/master/Images/NCBI_g3.png" alt="Filtering the queries for complete genomes" width="900"></p>

#### SSRG workflow
###### A case example with Minimap2, VarScan2 and 3 *Streptococcus pneumoniae* genomes

1. To download FASTA and GenBank files from a small .csv list of 3 *Streptococcus pneumoniae* genomes from NCBI with [queryNCBI.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNCBI.pl), type:
```bash
queryNCBI.pl \
   -l Examples/S_pneumoniae_3.csv \
   -o DATASETS \
   -fa \
   -gb
```

Options for [queryNCBI.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNCBI.pl) are:
```
-l (--list)	TAB/CSV-delimited list from NCBI
-o (--outdir)	Output directory [Default: ./]
-fa (--fasta)	Retrieve fasta files
-gb (--genbank)	Retrieve GenBank annotation files (.gbk; if available)
-gff (--gff3)	Retrieve GFF3 annotation files (.gff; if available)
-p (--protein)	Retrieve protein sequences (.faa; if available)
-cds		Retrieve protein coding sequences (.fna; if available)
```

2. FASTQ datasets can be generated from the downloaded genomes with [SSRG.pl](https://github.com/PombertLab/SSRG/blob/master/SSRG.pl).
```bash
SSRG.pl \
   -f DATASETS/*.fasta \
   -o FASTQ \
   -r 250
```

Options for [SSRG.pl](https://github.com/PombertLab/SSRG/blob/master/SSRG.pl) are:
```
 -f (--fasta)		Fasta/multifasta file(s)
 -o (--outdir)		Output directory [Default: ./]
 -l (--list)		List of fasta file(s), one per line
 -r (--readsize)	Desired read size [default: 150]
 -t (--type)		Read type: single (SE) or paired ends (PE) [default: PE]
 -i (--insert)		PE insert size [default: 350]
 -s (--sdev)		PE insert size standard deviation (in percentage) [default: 10]
 -c (--coverage)	Desired sequencing depth [default: 50]
 -qs (--qscore)		Quality score associated with each base [default: 30]
 -q64          		Use the old Illumina Q64 FastQ format instead of the default Q33 Sanger/Illumina 1.8+ encoding
```

3. To test the read-mapping step with [get_SNPs.pl](https://github.com/PombertLab/SSRG/blob/master/get_SNPs.pl) but skip variant calling, type:
```bash
## Running get_SNPs.pl with 8 threads, 16 Gb RAM (change according to your settings)
get_SNPs.pl \
	-threads 8 \
	-mem 16 \
	-fa DATASETS/*.fasta \
	-pe1 FASTQ/*R1.fastq \
	-pe2 FASTQ/*R2.fastq \
	-o RESULTS \
	-mapper minimap2 \
	-preset sr \
	-rmo \
	-bam
```
4. To test the read-mapping and variant calling steps with [get_SNPs.pl](https://github.com/PombertLab/SSRG/blob/master/get_SNPs.pl), type:
```bash
## Running get_SNPs.pl with 8 threads, 16 Gb RAM (change according to your settings)
get_SNPs.pl \
	-threads 8 \
	-mem 16 \
	-fa DATASETS/*.fasta \
	-pe1 FASTQ/*R1.fastq \
	-pe2 FASTQ/*R2.fastq \
	-o RESULTS \
	-mapper minimap2 \
	-preset sr \
	-bam \
	-caller varscan2 \
	-var VarScan.v2.4.4.jar ## Replace jar file location accordingly
```

Options for [get_SNPs.pl](https://github.com/PombertLab/SSRG/blob/master/get_SNPs.pl) are:
```
-h (--help)	Display this list of options
-v (--version)	Display script version
-o (--outdir)	Output directory [Default: ./]

# Mapping options
-fa (--fasta)			Reference genome(s) in fasta file
-fq (--fastq)			Fastq reads (single ends) to be mapped against reference(s)
-pe1				Fastq reads #1 (paired ends) to be mapped against reference(s)
-pe2				Fastq reads #2 (paired ends) to be mapped against reference(s)
-mapper				Read mapping tool: bowtie2, minimap2, ngmlr or hisat2 [default: minimap2]
-threads			Number of processing threads [default: 16]
-mem				Max total memory for samtools (in Gb) [default: 16] ## mem/threads = memory per thread
-bam				Keeps BAM files generated
-idx (--index)			Type of bam index generated (bai or csi) [default = csi] ## .bai not compatible with -mem
-sam				Keeps SAM files generated; SAM files can be quite large
-rmo (--read_mapping_only)	Do not perform variant calling; useful when only interested in bam/sam files and/or mapping stats
-ns (--no_stats)		Do not calculate stats; stats can take a while to compute for large eukaryote genomes

# Mapper-specific options
-X				BOWTIE2 - Maximum paired ends insert size [default: 750]
-preset				MINIMAP2 - Preset: sr, map-ont, map-pb or asm20 [default: sr]
-algo				BWA - Mapping algorithm:  bwasw, mem, samse [default: mem]

### Variant calling options ###
-caller				[default: varscan2]	## Variant caller: varscan2, bcftools or freebayes
-type				[default: snp]		## snp, indel, or both
-ploidy				[default: 1]		## FreeBayes/BCFtools option; change ploidy (if needed)

# VarScan2 parameters - see http://dkoboldt.github.io/varscan/using-varscan.html
-var				[default: VarScan.v2.4.4.jar]	## Which varscan jar file to use
-mc (--min-coverage)		[default: 15]		## Minimum read depth at a position to make a call
-mr (--min-reads2)		[default: 15]		## Minimum supporting reads at a position to call variants
-maq (--min-avg-qual)		[default: 28]		## Minimum base quality at a position to count a read
-mvf (--min-var-freq)		[default: 0.7]		## Minimum variant allele frequency threshold
-mhom (--min-freq-for-hom)	[default: 0.75]		## Minimum frequency to call homozygote
-pv (--p-value)			[default: 1e-02]	## P-value threshold for calling variants 
-sf (--strand-filter)		[default: 0]		## 0 or 1; 1 ignores variants with >90% support on one strand
```

5. To checking for synonymous/non-synonymous SNPs against a reference genome (*e.g. Streptococcus pneumoniae* R6) with [synonymy.pl](https://github.com/PombertLab/SSRG/blob/master/synonymy.pl), type:
```bash
synonymy.pl \
	-gcode 11 \
	-fa DATASETS/Streptococcus_pneumoniae_R6.fasta \
	-ref DATASETS/Streptococcus_pneumoniae_R6.gbk \
	-format gb \
	-vcf RESULTS/minimap2.varscan2.VCFs/*.vcf \
	-o SYNONYMY \
	-v
```

Options for [synonymy.pl](https://github.com/PombertLab/SSRG/blob/master/synonymy.pl) are:
```
-fa (--fasta)	Reference genome in fasta format
-r (--ref)	Reference genome annotation in GBK or GFF format
-f (--format)	Reference file format; gb or gff [default: gff]
-vcf (--vcf)	SNPs in Variant Calling Format (VCF)
-o (--outdir)	Output directory [Default: ./]
-p (--prefix)	Table prefix [Default: synonymy]
-gc (--gcode)	NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		# For complete list, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-v (--verbose)	Adds verbosity
```

#### Mash workflow
###### Case example with Mash, R and 75 *Streptococcus pneumoniae* genomes
1. To download FASTA and GenBank files from a CSV list of 75 *Streptococcus pneumoniae* genomes, type:
```bash
queryNCBI.pl \
   -l Examples/S_pneumoniae_75.csv \
   -o DATASETS \
   -fa
```

Options for [queryNCBI.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNCBI.pl) are:
```
-l (--list)	TAB/CSV-delimited list from NCBI
-o (--outdir)	Output directory [Default: ./]
-fa (--fasta)	Retrieve fasta files
-gb (--genbank)	Retrieve GenBank annotation files (.gbk; if available)
-gff (--gff3)	Retrieve GFF3 annotation files (.gff; if available)
-p (--protein)	Retrieve protein sequences (.faa; if available)
-cds		Retrieve protein coding sequences (.fna; if available)
```

2. To run Mash on the downloaded FASTA files with [run_Mash.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/run_Mash.pl), type:
```bash
run_Mash.pl \
   -f DATASETS/*.fasta \
   -o MASH \
   -n S_pneumoniae_75.mash
```

Options for [run_Mash.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/run_Mash.pl) are:
```
-f (--fasta)	Reference genome(s) in fasta file
-o (--outdir)	Output directory [Default: ./]
-n (--name)	Output file name [Default: Mash.mash]
-s (--sort)	Sort Mash output by decreasing order of similarity
```

3. To convert the Mash output to a distance matrix with [MashToDistanceMatrix.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/MashToDistanceMatrix.pl), type:
```bash
MashToDistanceMatrix.pl  \
	-i MASH/S_pneumoniae_75.mash \
	-o MASH/ \
	-p S_pneumoniae_75 \
	-f tsv
```

Options for [MashToDistanceMatrix.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/MashToDistanceMatrix.pl) are:
```
-i (--input)	Mash output file
-o (--outdir)	Output folder [Default: ./]
-p (--prefix)	Output file prefix [Default: Mash]
-f (--format)	Output format; csv or tsv [Default: csv]
```

4. To generate a quick neighbor-joining tree with [MashR_plotter.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/MashR_plotter.pl), type:
```bash
MashR_plotter.pl \
	-i MASH/S_pneumoniae_75.tsv \
	-if tsv \
	-t tree \
	-newick S_pneumoniae_75.tre \
	-f pdf \
	-o PLOTS \
	-n S_pneumoniae_75_NJ_tree \
	-he 20
```
The PDF generated should be similar to [S_pneumoniae_75_NJ_tree.pdf](https://github.com/PombertLab/SSRG/blob/master/Images/S_pneumoniae_75_NJ_tree.pdf) from the [Images/](https://github.com/PombertLab/SSRG/tree/master/Images) directory.

5. To generate a quick heatmap with [MashR_plotter.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/MashR_plotter.pl), type:
```bash
MashR_plotter.pl \
	-i MASH/S_pneumoniae_75.tsv \
	-if tsv \
	-t heatmap \
	-f pdf \
	-o PLOTS \
	-n S_pneumoniae_75_heatmap \
	-colors white cyan magenta \
	-he 20 \
	-wd 20
```
The PDF generated should be similar to [S_pneumoniae_75_heatmap.pdf](https://github.com/PombertLab/SSRG/blob/master/Images/S_pneumoniae_75_heatmap.pdf) from the [Images/](https://github.com/PombertLab/SSRG/tree/master/Images) directory.

6. To generate a quick t-SNE multidimensional reduction plot with [MashR_plotter.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/MASH/MashR_plotter.pl), type:
```bash
MashR_plotter.pl \
	-i MASH/S_pneumoniae_75.tsv \
	-if tsv \
	-t cluster \
	-m tsne \
	-pe 15 \
	-cmode terrain \
	-f pdf \
	-o PLOTS \
	-n S_pneumoniae_75_tSNE \
	-he 20 \
	-wd 20 \
	-lb \
	-fs 25
```
The PDF generated should be similar, ***but not identical***, to [S_pneumoniae_75_tSNE.pdf](https://github.com/PombertLab/SSRG/blob/master/Images/S_pneumoniae_75_tSNE.pdf) from the [Images/](https://github.com/PombertLab/SSRG/tree/master/Images) directory. t-SNE graphs are generated using random seeds, which affect how the distances are represented in 2D.

## Misc scripts
The script [queryNuccore.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNuccore.pl) retrieves genomes/proteins from the NCBI Nucleotide database using a list of accession numbers (one per line). This script uses [efetch](https://dataguide.nlm.nih.gov/edirect/efetch.html) from NCBI's [E-utilities](http://www.ncbi.nlm.nih.gov/books/NBK25499/) suite.

To download FASTA and GenBank files with [queryNuccore.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNuccore.pl), type:
```Bash
queryNuccore.pl \
  -l Examples/list_example_queryNuccore.txt
  -o DATASETS \
  -fa \
  -gb
```

Options for [queryNuccore.pl](https://github.com/PombertLab/SSRG/blob/master/Tools/NCBI/queryNuccore.pl) are:
```
-l (--list)		Accession numbers list, one accession per line
-o (--outdir)		Output folder [Default: ./]
-db (--database)	NCBI database to be queried [Default: nuccore]
-fa (--fasta)		Reference genome(s) in fasta format
-gb (--genbank)		Reference genome(s) in GenBank format
-sqn (--sequin)		Reference genome(s) in Sequin ASN format
-p (--proteins)		Protein sequences (amino acids)
-c (--cds)		Protein sequences (nucleotides)
```

The script [get_SRA.pl](https://github.com/PombertLab/SSRG/tree/master/Tools/NCBI) downloads data from the [NCBI Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) and converts it to FASTQ format using fasterq-dump from the [NCBI SRA toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software).

To download NCBI SRA datasets with [get_SRA.pl](https://github.com/PombertLab/SSRG/tree/master/Tools/NCBI) using 10 threads (-t 10), type:
```Bash
get_SRA.pl \
   -t 10 \
   -o FASTQ \
   -l Examples/list_example_get_SRA.txt \
   -p \
   -v
```

Note that this process can take quite a while depending on the size of the requested datasets and the available bandwidth.

Options for [get_SRA.pl](https://github.com/PombertLab/SSRG/tree/master/Tools/NCBI) are:
```
-t (--threads)	Number of CPU threads to use [Default: 10]
-o (--outdir)	Output directory [Default: ./]
-l (--list)	List(s) of SRA accesssion numbers, one accession number per line
-p (--progess)	Show progess
-v (--verbose)	Add verbosity
-f (--force)	Force to overwrite existing file(s)
```

## Funding and acknowledgments
This work was supported in part by the National Institute of Allergy and Infectious Diseases of the National Institutes of Health (award number R15AI128627) to Jean-Francois Pombert. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## References

1.	Li H. **Minimap2: pairwise alignment for nucleotide sequences.** *Bioinformatics.* 2018. 34, 3094–3100. PMID: [29750242](https://pubmed.ncbi.nlm.nih.gov/29750242) PMCID: [PMC6137996](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6137996) DOI: [10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)

2.	Langmead B, Salzberg SL. **Fast gapped-read alignment with Bowtie 2.** *Nat Methods.* 2012. 9, 357–359. PMID: [22388286](https://pubmed.ncbi.nlm.nih.gov/22388286) PMCID: [PMC3322381](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3322381) DOI: [10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)

3.	Kim D, Langmead B, Salzberg SL. **HISAT: a fast spliced aligner with low memory requirements.** *Nat. Methods.* 2015. 12, 357–360. PMID: [25751142](https://pubmed.ncbi.nlm.nih.gov/25751142) PMCID: [PMC4655817](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4655817) DOI: [10.1038/nmeth.3317](https://doi.org/10.1038/nmeth.3317)

4.	Sedlazeck FJ, Rescheneder P, Smolka M, Fang H, Nattestad M, von Haeseler A,  Schatz MC. **Accurate detection of complex structural variations using single-molecule sequencing.** *Nat. Methods.* 2018. 15, 461–468. PMID: [29713083](https://pubmed.ncbi.nlm.nih.gov/29713083) PMCID: [PMC5990442](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5990442) DOI: [10.1038/s41592-018-0001-7](https://doi.org/10.1038/s41592-018-0001-7)

5.	Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R. **The Sequence Alignment/Map format and SAMtools.** *Bioinformatics.* 2009. 25, 2078–2079. PMID: [19505943](https://pubmed.ncbi.nlm.nih.gov/19505943) PMCID: [PMC2723002](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002) DOI: [10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

6.	Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan, MD, Lin L, Miller CA, Mardis ER, Ding L, Wilson RK. **VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing.** *Genome Res.* 2012. 22, 568–576. PMID: [22300766](https://pubmed.ncbi.nlm.nih.gov/22300766) PMCID: [PMC3290792](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3290792) DOI: [10.1101/gr.129684.111](https://doi.org/10.1101/gr.129684.111)

7.	Narasimhan V, Danecek P, Scally A, Xue Y, Tyler-Smith C, Durbin R. **BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data.** *Bioinformatics.* 2016. 32, 1749–1751. PMID: [26826718](https://pubmed.ncbi.nlm.nih.gov/26826718) PMCID: [PMC4892413](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4892413) DOI: [10.1093/bioinformatics/btw044](https://doi.org/10.1093/bioinformatics/btw044)

8.	Ondov BD, Treangen TJ, Mallonee AB, Bergman NH, Koren S, Phillippy AM. **Fast genome and metagenome distance estimation using MinHash.** 2016. *Genome Biol.* 17, 132. PMID: [27323842](https://pubmed.ncbi.nlm.nih.gov/27323842) PMCID: [PMC4915045](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4915045) DOI: [10.1186/s13059-016-0997-x](https://doi.org/10.1186/s13059-016-0997-x)

9.	Garrison E, Marth G. **Haplotype-based variant detection from short-read sequencing.** *arXiv* Prepr. 2012. 9, [arXiv1207.3907](https://arxiv.org/abs/1207.3907).

10.	R Core Team. **R: A language and environment for statistical computing.** 2021. [URL](https://www.r-project.org/).

11.	van der Maaten LJP. **Barnes-Hut-SNE.** *arxiv.org* 2013. [1301.3342v2](https://arxiv.org/abs/1301.3342)

12.	van der Maaten LJP,  Hinton GE. **Visualizing High-Dimensional Data Using t-SNE.** *J. Mach. Learn. Res.* 2008. 9, 2579–2605. [URL](https://jmlr.csail.mit.edu/papers/volume9/vandermaaten08a/vandermaaten08a.pdf)
