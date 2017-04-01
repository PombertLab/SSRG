#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2016
## Version 1.4

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## User-defined environment variables.
## Change these to reflect your settings or leave blank (e.g. $bwa = '';) if set in $PATH
my $samtools = '/opt/samtools-1.3.1/bin/';	## Path to samtools 1.3.1+ - http://www.htslib.org/
my $bcftools = '/opt/bcftools-1.3.1/';		## Path to bcftools 1.3.1+ - http://www.htslib.org/
my $bwa = '/usr/bin/';				## Path to BWA - http://bio-bwa.sourceforge.net/
my $bowtie2 = '/opt/bowtie2-2.3.1/';		## Path to Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
my $hisat2 = '/opt/hisat2-2.0.4/';		## Path to HISAT2 - https://ccb.jhu.edu/software/hisat2/index.shtml
my $freebayes = '/opt/freebayes/bin/';		## Path to FreeBayes -  https://github.com/ekg/freebayes
my $varscan = '/opt/varscan/';			## Path to VarScan2 jar file to use -  https://github.com/dkoboldt/varscan
my $mash = '/opt/Mash/bin/';			## Path to Mash - https://github.com/marbl/Mash
my $varjar = 'VarScan.v2.4.3.jar';		## Define which VarScan2 jar file to use

## Usage definition
my $usage = "\nUSAGE = perl get_SNPs.pl [options]\n
EXAMPLE (simple): get_SNPs.pl -fa *.fasta -fq *.fastq
EXAMPLE (advanced): get_SNPs.pl --fasta *.fasta --fastq *.fastq --mapper bowtie2 --caller freebayes --threads 16
EXAMPLE (paired ends): get_SNPs.pl --fasta *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --X 1000 --mapper bowtie2 --caller freebayes --threads 16\n";
my $hint = "Type get_SNPs.pl -h (--help) for list of options\n";
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:

-h (--help)	Display this list of options

## Genetic distances
-mh	Evaluate genetic distances using Mash (Ondov et al. DOI: 10.1186/s13059-016-0997-x)
-out	Output file name [default: Mash.txt]
-sort	Sort Mash output by decreasing order of similarity

## Mapping options
-fa (--fasta)	Reference genome(s) in fasta file
-fq (--fastq)	Fastq reads (single ends) to be mapped against reference(s)
-pe1		Fastq reads #1 (paired ends) to be mapped against reference(s)
-pe2		Fastq reads #2 (paired ends) to be mapped against reference(s)
-X		Maximum paired ends insert size (for bowtie2) [default: 750]
-mapper		Read mapping tool: bwa, bowtie2 or hisat2 [default: bowtie2]
-caller		Variant caller: varscan2, bcftools or freebayes [default: varscan2]
-algo		BWA mapping algorithm:  bwasw, mem, samse [default: bwasw]
-threads	Number of processing threads [default: 16]
-bam		Keeps BAM files generated
-sam		Keeps SAM files generated; SAM files can be quite large

## VarScan2 parameters (see http://dkoboldt.github.io/varscan/using-varscan.html)
-indel				Calculates indels 	## Runs mpileup2indel
-mc (--min-coverage)		[default: 15]		## Minimum read depth at a position to make a call
-mr (--min-reads2)		[default: 5]		## Minimum supporting reads at a position to call variants
-maq (--min-avg-qual)		[default: 28]		## Minimum base quality at a position to count a read
-mvf (--min-var-freq)		[default: 0.2]		## Minimum variant allele frequency threshold
-mhom (--min-freq-for-hom)	[default: 0.75]		## Minimum frequency to call homozygote
-pv (--p-value)			[default: 1e-02]	## P-value threshold for calling variants 
-sf (--strand-filter)		[default: 0]		## 0 or 1; 1 ignores variants with >90% support on one strand

## FreeBayes/BCFtools (see https://github.com/ekg/freebayes/; https://samtools.github.io/bcftools/bcftools.html)
-ploidy			[default: 1]		## Change ploidy (if needed)

END_OPTIONS

my $help ='';
## Genetic distances
my $mh = '';
my $out = 'Mash.txt';
my $sort = '';
## Mapping
my $mapper = 'bowtie2';
my $caller = 'varscan2';
my $algo = 'bwasw';
my $threads = 16;
my $bam = '';
my $sam = '';
my @fasta;
my @fastq;
my @pe1;
my @pe2;
my $maxins = '750'; 
## VarScan
my $indel = '';
my $mc = 15;
my $mr = 5;
my $maq = 28;
my $mvf = '0.2';
my $mhom = '0.75';
my $pv = '1e-02';
my $sf = 0;
## FreeBayes/BCFtools
my $ploidy = 1;

GetOptions(
	'h|help' => \$help,
	'mh' => \$mh,
	'out=s' => \$out,
	'sort' => \$sort,
	'mapper=s' => \$mapper,
	'caller=s' => \$caller,
	'algo=s' => \$algo,
	'threads=i' => \$threads,
	'bam' => \$bam,
	'sam' => \$sam,
	'fa|fasta=s@{1,}' => \@fasta,
	'fq|fastq=s@{1,}' => \@fastq,
	'pe1=s@{1,}' => \@pe1,
	'pe2=s@{1,}' => \@pe2,
	'X=i' => \$maxins,
	'indel' => \$indel,
	'mc|min-coverage=i' => \$mc,
	'mr|min-reads2=i' => \$mr,
	'maq|min-avg-qual=i' => \$maq,
	'mvf|min-var-freq=s' => \$mvf,
	'mhom|min-freq-for-hom=s' => \$mhom,
	'pv|p-value=s' => \$pv,
	'sf|strand-filter=i' => \$sf,
	'ploidy=i' => \$ploidy,
);

if ($help){die "$usage\n$options";}

## Timestamps and logs
my $start = localtime();
my $tstart = time;
my $todo = 0;
if (@fastq){$todo += scalar(@fastq)*scalar(@fasta);}
if (@pe1 && @pe2){$todo += scalar(@pe1)*scalar(@fasta);}
open LOG, ">time.$mapper.$caller.log"; ## Keep track of running time
print LOG "Mapping/SNP calling started on: $start\n";
print LOG "A total of $todo pairwise comparisons will be performed\n";
open MAP, ">>mapping.$mapper.log"; ## Keep track of read mapper STDERR messages
my $comparison = 0;

## Genetic distances with MASH
if ($mh){
	system "echo Running Mash genetic distance analysis...";
	system "$mash"."mash sketch -o reference @fasta";
	system "$mash"."mash dist reference.msh @fasta > $out";
	if ($sort){
		system "echo Sorting out Mash results - See $out.sorted";
		system "sort -gk3 $out > $out.sorted";
	}
 }

## Running read mapping/SNP calling
my $fasta = undef; my $fastq = undef; my $file = undef;
foreach (@fasta){ ## Creating indexes
	if ($mapper eq 'bwa'){system "$bwa"."bwa index -a is $_";}
	elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2-build --threads $threads $_ $_.bt2";}
	elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2-build $_ $_.ht";}
}
my $index_time = time - $tstart;
print LOG "Time required to create all indexes: $index_time seconds\n";

if (@fastq){ ## Single ends mode
	foreach (@fastq){
		$fastq = $_;
		$file = $fastq;
		foreach (@fasta){
			$fasta = $_;
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";
			if ($mapper eq 'bwa'){system "$bwa"."bwa $algo -t $threads $fasta $fastq -f $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";} ## Appending STDERR to mapping.$mapper.log"
			elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2 -p $threads -x $fasta.bt2 -U $fastq -S $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";} 
			elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2 -p $threads --phred33 -x $fasta.ht -U $fastq --no-spliced-alignment -S $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";}
			samtools();
			variant();
			stats(); logs(); ## Printing stats
			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $file.$fasta.$mapper.bam";}
			unless ($sam) {system "rm $file.$fasta.$mapper.sam";}
		}
	}
}

if (@pe1 && @pe2){ ## Paired ends mode
	while (my $pe1 = shift@pe1){
		my $pe2 = shift@pe2;
		$file = $pe1;
		foreach (@fasta){
			$fasta = $_;
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$pe1 + $pe2 vs. $fasta\n";
			if ($mapper eq 'bwa'){system "$bwa"."bwa $algo -t $threads $fasta $pe1 $pe2 -f $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";} ## Appending STDERR to mapping.$mapper.log"
			elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2 -p $threads -x $fasta.bt2 -X $maxins -1 $pe1 -2 $pe2 -S $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";} 
			elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2 -p $threads --phred33 -x $fasta.ht -1 $pe1 -2 $pe2 --no-spliced-alignment -S $file.$fasta.$mapper.sam 2>> mapping.$mapper.log";}
			samtools();
			variant();
			stats(); logs(); ## Printing stats
			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $file.$fasta.$mapper.bam";}
			unless ($sam) {system "rm $file.$fasta.$mapper.sam";}
		}
	}
}

## Cleaning up
if ($mapper eq 'bwa'){system "mkdir $mapper.indexes; mv *.amb *.ann *.bwt *.fai *.pac *.sa $mapper.indexes/";}
elsif ($mapper eq 'bowtie2'){system "mkdir $mapper.indexes; mv *.bt2 *.fai $mapper.indexes/";}
elsif ($mapper eq 'hisat2'){system "mkdir $mapper.indexes; mv *.ht2 *.fai $mapper.indexes/";}
system "mkdir $mapper.$caller.VCFs; mv *.vcf $mapper.$caller.VCFs/";

## Printing timestamps and logs
my $end = localtime();
my $time_taken = time - $tstart;
my $average_time = sprintf("%.1f", $time_taken/$comparison);
print "\nMapping/SNP calling started on: $start\n";
print "Mapping/SNP calling ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";
print LOG "Mapping/SNP calling ended on: $end\n";
print LOG "Time elapsed: $time_taken seconds\n";
print LOG "Average time per pairwise comparison: $average_time seconds\n";

### Subroutines ###
sub samtools{
	system "$samtools"."samtools view -@ $threads -bS $file.$fasta.$mapper.sam -o $file.$fasta.bam";
	system "$samtools"."samtools sort -@ $threads -o $file.$fasta.$mapper.bam $file.$fasta.bam";
	system "$samtools"."samtools depth -aa $file.$fasta.$mapper.bam > $file.$fasta.$mapper.coverage"; ## Printing per base coverage information
	system "rm $file.$fasta.bam"; ## Discarding unsorted BAM file
}

sub variant{
	if ($caller eq 'varscan2'){
		system "$samtools"."samtools mpileup -f $fasta $file.$fasta.$mapper.bam | java -jar $varscan$varjar mpileup2snp --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $file.$fasta.$mapper.SNP.vcf";
		if ($indel){system "$samtools"."samtools mpileup -f $fasta $file.$fasta.$mapper.bam | java -jar $varscan$varjar mpileup2indel --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $file.$fasta.$mapper.indel.vcf";}
	}
	elsif ($caller eq 'bcftools'){
		system "$samtools"."samtools mpileup -ugf $fasta $file.$fasta.$mapper.bam | $bcftools"."bcftools call -vmO v -V indels --ploidy $ploidy --output $file.$fasta.$mapper.SNP.vcf";
		if ($indel){system "$samtools"."samtools mpileup -ugf $fasta $file.$fasta.$mapper.bam | $bcftools"."bcftools call -vmO v -V snps --ploidy $ploidy --output $file.$fasta.$mapper.indel.vcf";}
	}
	elsif ($caller eq 'freebayes'){ ## single thread only, parallel version behaving wonky
		system "$samtools"."samtools index $file.$fasta.$mapper.bam";
		system "$freebayes"."freebayes -f $fasta -p $ploidy $file.$fasta.$mapper.bam > $file.$fasta.$mapper.$caller.SNP.vcf";
		system "rm $file.$fasta.$mapper.bam.bai";
	}
}

sub stats{
	open COV, "<$file.$fasta.$mapper.coverage";
	open SN, "<$file.$fasta.$mapper.SNP.vcf";
	open STATS, ">$file.$fasta.$mapper.stats";
	my $total = 0; my $covered = 0; my $nocov = 0; my $max = 0; my $sumcov; my $sn =0;
	while (my $line = <SN>){
		if ($line =~ /^#/){next;}
		else {$sn++;}
	}
	while (my $line = <COV>){
		chomp $line;
		$total++; 
		if ($line =~ /^\S+\s+\d+\s+(\d+)/){
			my $coverage = $1;
			$sumcov += $coverage;
			if ($coverage >= 1) {
				$covered++;
				if ($coverage > $max){$max = $coverage;}
			}
			else {$nocov++;}
		}
	}
	print STATS "FASTQ file(s) used: $file (and mate, if PE)\n";
	print STATS "Reference fasta file used: $fasta\n";
	print STATS "\nTotal number of bases in the reference genome\t$total bp\n"."Number of bases covered by at least one read\t$covered\n". "Number of bases without coverage\t$nocov\n";
	print STATS "Maximum sequencing depth\t$max"."X\n";
	my $avc = sprintf("%.2f", ($sumcov/$total));
	print STATS "Average sequencing depth\t$avc"."X\n";
	if ($total == $covered){print STATS "Sequencing breadth (percentage of bases covered by at least one read)\t100%\n";}
	if ($total != $covered){my $av_cov = sprintf("%.2f%%", ($covered/$total)*100); print STATS "Sequencing breadth (percentage of bases covered by at least one read)\t$av_cov\n";}
	my $snkb = sprintf("%.2f", ($sn/$covered)*1000);
	print STATS "Total number of SNPs found: $sn\n";
	print STATS "Average number of SNPs per Kb: $snkb\n";
	print STATS "\n## SAMTOOLS flagstat metrics\n";
	close STATS;
	system "$samtools"."samtools flagstat $file.$fasta.$mapper.bam >> $file.$fasta.$mapper.stats";
}

sub logs{
	my $run_time = time - $tstart; my $mend = localtime();
	$comparison++;
	print LOG "Comparison # $comparison : $file (and mate, if PE) vs. $fasta - cumulative time elapsed: $run_time seconds\n";
	print MAP "\n".'###'." Mapping ended on $mend\n\n";
}
