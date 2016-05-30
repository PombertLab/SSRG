#!/usr/bin/perl
## Pombert JF, Illinois tech - 2016
## Version 1.2b

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## User-defined environment variables. Change these to reflect your settings.
my $bwa = '/usr/bin/'; ## Path to BWA - http://bio-bwa.sourceforge.net/
my $bowtie2 = '/opt/bowtie2-2.2.9/'; ## Path to Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
my $hisat2 = '/opt/hisat2-2.0.4/'; ## Path to HISAT2 - https://ccb.jhu.edu/software/hisat2/index.shtml
my $varscan = '/opt/VarScan/VarScan.v2.4.2.jar'; ## Define which VarScan2 jar file to use -  https://github.com/dkoboldt/varscan

## Defining options
my $usage = "
USAGE = perl get_SNPs.pl [options]

EXAMPLE: get_SNPs.pl --fasta *.fasta --fastq *.fastq --mapper bowtie2 --threads 16 --indel --sf 1 --mc 20

OPTIONS:

## Mapping options
--fasta		Reference genome(s) in fasta file
--fastq		Fastq reads to be mapped against reference(s)
--mapper	Read mapping tool: bwa, bowtie2 or hisat2 [default: bowtie2]
--algo		BWA mapping algorithm:  bwasw, mem, samse [default: bwasw]
--threads	Number of processing threads [default: 16]
--bam		Keeps BAM files generated
--sam		Keeps SAM files generated; SAM files can be quite large

## VarScan2 parameters (see http://dkoboldt.github.io/varscan/using-varscan.html)
--indel				Calculates indels 	## Runs mpileup2indel
--mc (--min-coverage)		[default: 8]		## Minimum read depth at a position to make a call
--mr (--min-reads2)		[default: 2]		## Minimum supporting reads at a position to call variants
--maq (--min-avg-qual)		[default: 15]		## Minimum base quality at a position to count a read
--mvf (--min-var-freq)		[default: 0.01]		## Minimum variant allele frequency threshold
--mhom (--min-freq-for-hom)	[default: 0.75]		## Minimum frequency to call homozygote
--pv (--p-value)		[default: 99e-02]	## P-value threshold for calling variants 
--sf (--strand-filter)		[default: 0]		## 0 or 1; 1 ignores variants with >90% support on one strand
";

die "$usage\n" unless@ARGV;

## Mapping
my $mapper = 'bowtie2';
my $algo = 'bwasw';
my $threads = 16;
my $bam = '';
my $sam = '';
my @fasta;
my @fastq;
## VarScan
my $indel = '';
my $mc = 8;
my $mr = 2;
my $maq = 15;
my $mvf = '0.01';
my $mhom = '0.75';
my $pv = '99e-02';
my $sf = 0;

GetOptions(
	'mapper=s' => \$mapper,
	'algo=s' => \$algo,
	'threads=i' => \$threads,
	'bam' => \$bam,
	'sam' => \$sam,
	'fasta=s@{1,}' => \@fasta,
	'fastq=s@{1,}' => \@fastq,
	'indel' => \$indel,
	'mc|min-coverage=i' => \$mc,
	'mr|min-reads2=i' => \$mr,
	'maq|min-avg-qual=i' => \$maq,
	'mvf|min-var-freq=s' => \$mvf,
	'mhom|min-freq-for-hom=s' => \$mhom,
	'pv|p-value=s' => \$pv,
	'sf|strand-filter=i' => \$sf,
);

## Timestamps and logs
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
			system "$bwa"."bwa $algo -t $threads $fasta $fastq -f $fastq.$fasta.$mapper.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.$mapper.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.$mapper";
			system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2snp --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2indel --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam"; ## Discarding unsorted BAM file
			unless ($bam) {system "rm $fastq.$fasta.$mapper.bam";}
			unless ($sam) {system "rm $fastq.$fasta.$mapper.sam";}
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
			system "$bowtie2"."bowtie2 -p $threads -x $fasta.bt2 -U $fastq -S $fastq.$fasta.$mapper.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.$mapper.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.$mapper";
			system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2snp --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2indel --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam"; ## Discarding unsorted BAM file
			unless ($bam) {system "rm $fastq.$fasta.$mapper.bam";}
			unless ($sam) {system "rm $fastq.$fasta.$mapper.sam";}
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
			system "$hisat2"."hisat2 -p $threads --phred33 -x $fasta.ht -U $fastq --no-spliced-alignment -S $fastq.$fasta.$mapper.sam 2>> mapping.$mapper.log"; ## Appending STDERR to mapping.$mapper.log"
			system "samtools view -bS $fastq.$fasta.$mapper.sam -o $fastq.$fasta.bam";
			system "samtools sort $fastq.$fasta.bam $fastq.$fasta.$mapper";
			system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2snp --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.SNP.vcf";
			if ($indel){system "samtools mpileup -f $fasta $fastq.$fasta.$mapper.bam | java -jar $varscan mpileup2indel --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $fastq.$fasta.$mapper.indel.vcf";}
			## Cleaning temp files
			system "rm $fastq.$fasta.bam"; ## Discarding unsorted BAM file
			unless ($bam) {system "rm $fastq.$fasta.$mapper.bam";}
			unless ($sam) {system "rm $fastq.$fasta.$mapper.sam";}
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