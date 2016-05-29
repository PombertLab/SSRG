#!/usr/bin/perl
## Pombert JF, Illinois tech - 2016
## Version 1.2a

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## User-defined environment variables. Change these to reflect your settings.
my $bwa = '/usr/bin/'; ## Path to BWA - http://bio-bwa.sourceforge.net/
my $bowtie2 = '/opt/bowtie2-2.2.9/'; ## Path to Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
my $hisat2 = '/opt/hisat2-2.0.1-beta/'; ## Path to HISAT2 - https://ccb.jhu.edu/software/hisat2/index.shtml
my $varscan = '/opt/VarScan/VarScan.v2.4.2.jar'; ## Define which VarScan2 jar file to use -  https://github.com/dkoboldt/varscan

## Defining options
my $usage = "
USAGE = perl get_SNPs.pl [options]

EXAMPLE: get_SNPs.pl --fasta *.fasta --fastq *.fastq --mapper bowtie2 --threads 16 --indel --bam

OPTIONS:
--fasta		Reference genome(s) in fasta file
--fastq		Fastq reads to be mapped against reference(s)
--mapper	Read mapping tool: bwa, bowtie2 or hisat2 (not recommended) [default: bowtie2]
--algo		BWA mapping algorithm:  bwasw, mem, samse [default: bwasw]
--threads	Number of processing threads: 2,4,8 ... [default: 16]
--indel		Calculates indels
--bam		Keeps BAM files generated
--sam		Keeps SAM files generated; SAM files can be quite large
";

die "$usage\n" unless@ARGV;

my $mapper = 'bowtie2';
my $algo = 'bwasw';
my $threads = 16;
my $bam = '';
my $sam = '';
my $indel = '';
my @fasta;
my @fastq;

GetOptions(
	'mapper=s' => \$mapper,
	'algo=s' => \$algo,
	'threads=i' => \$threads,
	'bam' => \$bam,
	'sam' => \$sam,
	'indel' => \$indel,
	'fasta=s@{1,}' => \@fasta,
	'fastq=s@{1,}' => \@fastq,
);

## Timestamps
my $start = localtime();
my $tstart = time;
my $todo = scalar(@fastq)*scalar(@fasta);
open LOG, ">time.$mapper.log"; ## Keep track of running time
print LOG "Mapping/SNP calling started on: $start\n";
print LOG "A total of $todo pairwise comparisons will be performed\n";
open MAP, ">>mapping.$mapper.log"; ## Keep track of read mapper STDERR messages
my $comparison = 0;

## Running BWA mapping
if ($mapper eq 'bwa'){
	## Creating Burrows-Wheeler BWA indexes
	foreach my $fasta (@fasta){system "$bwa"."bwa index -a is $fasta";}
	my $index_time = time - $tstart;
	print LOG "Time required to create all indexes: $index_time seconds\n";
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";
			system "$bwa"."bwa $algo -t $threads $fasta $fastq -f $fastq.$fasta.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.BWA.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.BWA.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam";
			unless ($bam) {system "rm $fastq.$fasta.sorted.bam";}
			unless ($sam) {system "rm $fastq.$fasta.sam";}
			## Logs
			my $run_time = time - $tstart; my $mend = localtime();
			$comparison++;
			print LOG "Comparison # $comparison : $fastq vs. $fasta - cumulative time elapsed: $run_time seconds\n";
			print MAP "\n".'###'." Mapping ended on $mend\n\n";
		}
	}
	## Cleaning up
	system "mkdir BWA_indexes; mv *.amb *.ann *.bwt *.fai *.pac *.sa BWA_indexes/";
	system "mkdir BWA_VCFs; mv *.vcf BWA_VCFs/";
}
## Running Bowtie2 mapping
elsif ($mapper eq 'bowtie2'){
	## Creating Burrows-Wheeler HISAT2 indexes
	foreach my $fasta (@fasta){system "$bowtie2"."bowtie2-build --threads $threads $fasta $fasta.bt2";}
	my $index_time = time - $tstart;
	print LOG "Time required to create all indexes: $index_time seconds\n";
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";
			system "$bowtie2"."bowtie2 -p $threads -x $fasta.bt2 -U $fastq -S $fastq.$fasta.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.Bowtie2.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.Bowtie2.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam";
			unless ($bam) {system "rm $fastq.$fasta.sorted.bam";}
			unless ($sam) {system "rm $fastq.$fasta.sam";}
			## Logs
			my $run_time = time - $tstart; my $mend = localtime();
			$comparison++;
			print LOG "Comparison # $comparison : $fastq vs. $fasta - cumulative time elapsed: $run_time seconds\n";
			print MAP "\n".'###'." Mapping ended on $mend\n\n";
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
	my $index_time = time - $tstart;
	print LOG "Time required to create all indexes: $index_time seconds\n";
	## Running the read mapping
	foreach my $fastq (@fastq){
		foreach my $fasta (@fasta){
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";
			system "$hisat2"."hisat2 -p $threads --phred33 -x $fasta.ht -U $fastq --no-spliced-alignment -S $fastq.$fasta.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.sorted";
			system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2snp --output-vcf --strand-filter 0 > $fastq.$fasta.HISAT2.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.sorted.bam | java -jar $varscan mpileup2indel --output-vcf --strand-filter 0 > $fastq.$fasta.HISAT2.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam";
			unless ($bam) {system "rm $fastq.$fasta.sorted.bam";}
			unless ($sam) {system "rm $fastq.$fasta.sam";}
			## Logs
			my $run_time = time - $tstart; my $mend = localtime();
			$comparison++;
			print LOG "Comparison # $comparison : $fastq vs. $fasta - cumulative time elapsed: $run_time seconds\n";
			print MAP "\n".'###'." Mapping ended on $mend\n\n";
		}
	}
	## Cleaning up
	system "mkdir HISAT2_indexes; mv *.ht2 *.fai HISAT2_indexes/";
	system "mkdir HISAT2_VCFs; mv *.vcf HISAT2_VCFs/";
}

my $end = localtime();
my $time_taken = time - $tstart;
my $average_time = sprintf("%.1f", $time_taken/$comparison);

print "\nMapping/SNP calling started on: $start\n";
print "Mapping/SNP calling ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";
print LOG "Mapping/SNP calling ended on: $end\n";
print LOG "Time elapsed: $time_taken seconds\n";
print LOG "Average time per pairwise comparison: $average_time seconds\n";