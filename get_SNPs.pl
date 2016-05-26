#!/usr/bin/perl
## Pombert JF, Illinois tech - 2016
## Version 1.2

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $bwa = '/usr/bin/'; ## Path to BWA
my $bowtie2 = '/opt/bowtie2-2.2.9/'; ## Path to Bowtie2
my $hisat2 = '/opt/hisat2-2.0.1-beta/'; ## Path to HISAT2
my $varscan = '/opt/VarScan/VarScan.v2.3.7.jar'; ## Define which VarScan2 jar file to use

my $usage = "
USAGE = perl get_SNPs.pl [options]

OPTIONS:
--fasta		reference genome(s) in fasta file
--fastq		fastq reads to be mapped against reference(s)
--mapper	bwa, bowtie2 or hisat2 (not recommended) [default: bwa]
--algo		bwa mapping algorithm [default: bwasw]
--threads	number of processing threads 2,4,8... [default: 16]
--indel		Type 'yes' to calculate Indels [default: no]
--temp		Type 'yes' to keep temporary SAM/BAM files [default: no], WARNING temporary files generated can be quite large
";

die "$usage\n" unless@ARGV;

## Defining options
my $mapper = 'bwa';
my $algo = 'bwasw';
my $threads = 16;
my $temp = 'no';
my $indel = 'no';
my @fasta;
my @fastq;

GetOptions(
	'mapper=s' => \$mapper,
	'algo=s' => \$algo,
	'threads=i' => \$threads,
	'temp=s' => \$temp,
	'indel=s' => \$indel,
	'fasta=s@{1,}' => \@fasta,
	'fastq=s@{1,}' => \@fastq,
);

## Running BWA mapping
if ($mapper eq 'bwa'){
	## Creating Burrows-Wheeler BWA indexes
	foreach my $fasta (@fasta){system "$bwa"."bwa index -a is $fasta";}
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			system "$bwa"."bwa $algo -t $threads $fasta $fastq -f $fastq.$fasta.sam";
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.BWA.SNP.vcf";
			if ($indel eq 'yes'){
				system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.BWA.indel.vcf";
			}
			if ($temp eq 'no'){system "rm *.bam *.sam";}
		}
	}
	## Cleaning up
	system "mkdir BWA_indexes; mv *.amb *.ann *.bwt *.fai *.pac *.sa BWA_indexes/";
	system "mkdir BWA_VCFs; mv *.vcf BWA_VCFs/";
}
## Running Bowtie2 mapping
elsif ($mapper eq 'bowtie2'){
	## Creating Burrows-Wheeler HISAT2 indexes
	foreach my $fasta (@fasta){system "$bowtie2"."bowtie2-build $fasta $fasta.bt2";}
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			system "$bowtie2"."bowtie2 -p $threads -x $fasta.bt2 -U $fastq -S $fastq.$fasta.sam";
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.Bowtie2.SNP.vcf";
			if ($indel eq 'yes'){
				system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.Bowtie2.indel.vcf";
			}
			if ($temp eq 'no'){system "rm *.bam *.sam";}
		}
	}
	## Cleaning up
	system "mkdir Bowtie2_indexes; mv *.bt2 *.fai Bowtie2_indexes/";
	system "mkdir Bowtie2_VCFs; mv *.vcf Bowtie2_VCFs/";
}
## Running HISAT2 mapping
elsif ($mapper eq 'hisat2'){
	## Creating Burrows-Wheeler HISAT2 indexes
	foreach my $fasta (@fasta){system "$hisat2"."hisat2-build $fasta $fasta.ht";}
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			system "$hisat2"."hisat2 -p $threads --phred33 -x $fasta.ht -U $fastq --no-spliced-alignment -S $fastq.$fasta.sam";
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.HISAT2.SNP.vcf";
			if ($indel eq 'yes'){
				system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.HISAT2.indel.vcf";
			}
			if ($temp eq 'no'){system "rm *.bam *.sam";}
		}
	}
	## Cleaning up
	system "mkdir HISAT2_indexes; mv *.ht2 *.fai HISAT2_indexes/";
	system "mkdir HISAT2_VCFs; mv *.vcf HISAT2_VCFs/";
}