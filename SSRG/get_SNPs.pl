#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2018
my $version = '1.7';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

## User-defined environment variables. ## Leave these blank (e.g. $samtools = '';) if programs are set in $PATH.
## Alternatively, you can insert here the installation directories (e.g. $samtools = '/opt/samtools-1.3.1/bin/') to reflect your settings.
my $samtools = '';		## Path to samtools 1.3.1+ - http://www.htslib.org/
my $bcftools = '';		## Path to bcftools 1.3.1+ - http://www.htslib.org/
my $bwa = '';			## Path to BWA - http://bio-bwa.sourceforge.net/
my $bowtie2 = '';		## Path to Bowtie2 - http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
my $hisat2 = '';		## Path to HISAT2 - https://ccb.jhu.edu/software/hisat2/index.shtml
my $freebayes = '';		## Path to FreeBayes -  https://github.com/ekg/freebayes
my $mash = '';			## Path to Mash - https://github.com/marbl/Mash

## Usage definition
my $hint = "Type get_SNPs.pl -h (--help) for list of options\n";
my $usage = <<'USAGE';

USAGE = perl get_SNPs.pl [options]

EXAMPLE (simple): get_SNPs.pl -fa *.fasta -fq *.fastq
EXAMPLE (advanced): get_SNPs.pl --fasta *.fasta --fastq *.fastq --mapper bowtie2 --caller varscan2 --type both --var ./VarScan.v2.4.3.jar --threads 16
EXAMPLE (paired ends): get_SNPs.pl --fasta *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --X 1000 --mapper bowtie2 --caller freebayes --threads 16
EXAMPLE (MASH): get_SNPs.pl -fa *.fasta -out Mash.txt -sort
USAGE
die "$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:
-h (--help)	Display this list of options
-v (--version)	Display script version

### Mapping options ###
-fa (--fasta)	Reference genome(s) in fasta file
-fq (--fastq)	Fastq reads (single ends) to be mapped against reference(s)
-pe1		Fastq reads #1 (paired ends) to be mapped against reference(s)
-pe2		Fastq reads #2 (paired ends) to be mapped against reference(s)
-X		Maximum paired ends insert size (for bowtie2) [default: 750]
-mapper		Read mapping tool: bwa, bowtie2 or hisat2 [default: bowtie2]
-algo		BWA mapping algorithm:  bwasw, mem, samse [default: bwasw]
-threads	Number of processing threads [default: 16]
-bam		Keeps BAM files generated
-sam		Keeps SAM files generated; SAM files can be quite large

### Variant caller options ###
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

### MASH Genetic distances ### OPTIONAL - see Ondov et al. DOI: 10.1186/s13059-016-0997-x
-mh		Evaluate genetic distances using Mash ## Quick, does not require read mapping...
-out		Output file name [default: Mash.txt]
-sort		Sort Mash output by decreasing order of similarity

END_OPTIONS

my $help =''; my $vn;
## Mapping
my $mapper = 'bowtie2';
my $algo = 'bwasw';
my $threads = 16;
my $sam = ''; my $bam = '';
my @fasta; my @fastq;
my @pe1; my @pe2;
my $maxins = '750';
## Variant caller
my $caller = 'varscan2';
my $type = 'snp';
my $ploidy = 1;
## VarScan-specific
my $varjar = '/opt/varscan/VarScan.v2.4.3.jar';	## Define which VarScan2 jar file to use; or use -var from the CMD line
my $mc = 15;
my $mr = 15;
my $maq = 28;
my $mvf = '0.7';
my $mhom = '0.75';
my $pv = '1e-02';
my $sf = 0;
## MASH
my $out = 'Mash.txt'; my $mh = ''; my $sort = ''; 

open MAP, ">>mapping.$mapper.log"; ## Keep track of read mapper STDERR messages
print MAP "COMMAND LINE:\nget_SNPs.pl @ARGV\n";

GetOptions(
	'h|help' => \$help, 'v|version', => \$vn,
	## Mapping
	'mapper=s' => \$mapper,
	'algo=s' => \$algo,
	'threads=i' => \$threads,
	'sam' => \$sam, 'bam' => \$bam,
	'fa|fasta=s@{1,}' => \@fasta, 'fq|fastq=s@{1,}' => \@fastq,
	'pe1=s@{1,}' => \@pe1, 'pe2=s@{1,}' => \@pe2,
	'X=i' => \$maxins,
	## Variant caller
	'caller=s' => \$caller, 'type=s' => \$type, 'ploidy=i' => \$ploidy,
	## VarScan-specific
	'var=s' => \$varjar,
	'mc|min-coverage=i' => \$mc,
	'mr|min-reads2=i' => \$mr,
	'maq|min-avg-qual=i' => \$maq,
	'mvf|min-var-freq=s' => \$mvf,
	'mhom|min-freq-for-hom=s' => \$mhom,
	'pv|p-value=s' => \$pv,
	'sf|strand-filter=i' => \$sf,
	## MASH
	'mh' => \$mh, 'out=s' => \$out,	'sort' => \$sort,
);

if ($help){die "$usage\n$options";} if ($vn){die "\nversion $version\n\n";}

## Timestamps and logs
my $start = localtime();
my $tstart = time;
my $todo = 0;
if (@fastq){$todo += scalar(@fastq)*scalar(@fasta);}
if (@pe1 && @pe2){$todo += scalar(@pe1)*scalar(@fasta);}
open LOG, ">time.$mapper.$caller.log"; ## Keep track of running time
print LOG "Mapping/SNP calling started on: $start\n";
print LOG "A total of $todo pairwise comparisons will be performed\n";
print_options();
my $comparison = 0;

## Genetic distances with MASH
if ($mh){
	system "echo Running Mash genetic distance analysis...";
	system "$mash"."mash sketch -o reference @fasta";
	system "$mash"."mash dist reference.msh @fasta > $out";
	if ($sort){system "echo Sorting out Mash results - See $out.sorted"; system "sort -gk3 $out > $out.sorted";}
	exit;
 }

## Running read mapping/SNP calling
my $fasta = undef; my $fastq = undef; my $file = undef; my $fa = undef; my $dir; my $qdir;
foreach (@fasta){ ## Creating indexes
	$fasta = $_; ($fa, $dir) = fileparse($fasta);
	if ($mapper eq 'bwa'){system "$bwa"."bwa index -a is $_";}
	elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2-build --threads $threads $_ $fa.bt2";}
	elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2-build $_ $fa.ht";}
}
my $index_time = time - $tstart;
print LOG "Time required to create all indexes: $index_time seconds\n";

if (@fastq){ ## Single ends mode
	foreach (@fastq){
		$fastq = $_;
		print "\n## FASTQ information:\n";
		($file, $qdir) = fileparse($fastq); print "FASTQ parsed as: $file\n"; print "FASTQ input DIR parsed as: $qdir\n"; ## Debugging print statement
		foreach (@fasta){
			print "\n## FASTA information:\n";
			$fasta = $_; ($fa, $dir) = fileparse($fasta); print "FASTA parsed as: $fa\n"; print "FASTA input DIR parsed as: $dir\n\n"; ## Debugging print statement
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";
			print "Mapping $fastq on $fasta with $mapper...\n";
			if ($mapper eq 'bwa'){system "$bwa"."bwa $algo -t $threads $fasta $fastq -f $file.$fa.$mapper.sam 2>> mapping.$mapper.log";} ## Appending STDERR to mapping.$mapper.log"
			elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2 -p $threads -x $fa.bt2 -U $fastq -S $file.$fa.$mapper.sam 2>> mapping.$mapper.log";} 
			elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2 -p $threads --phred33 -x $fa.ht -U $fastq --no-spliced-alignment -S $file.$fa.$mapper.sam 2>> mapping.$mapper.log";}
			samtools(); variant(); stats(); logs(); ## Printing stats
			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $file.$fa.$mapper.bam";}
			unless ($sam) {system "rm $file.$fa.$mapper.sam";}
		}
	}
}

if (@pe1 && @pe2){ ## Paired ends mode
	while (my $pe1 = shift@pe1){
		my $pe2 = shift@pe2;
		print "\n## FASTQ information:\n";
		($file, $qdir) = fileparse($pe1); print "R1 FASTQ parsed as: $file\n"; ## Debugging print statement
		my $r2 = fileparse($pe2); print "R2 FASTQ parsed as: $r2\n"; print "FASTQ input DIR parsed as: $qdir\n";## Debugging print statement
		foreach (@fasta){
			print "\n## FASTA information:\n";
			$fasta = $_; ($fa, $dir) = fileparse($fasta); print "FASTA parsed as: $fa\n"; print "FASTA input DIR parsed as: $dir\n\n"; ## Debugging print statement
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$pe1 + $pe2 vs. $fasta\n";
			print "Mapping $pe1 and $pe2 on $fasta with $mapper...\n";
			if ($mapper eq 'bwa'){system "$bwa"."bwa $algo -t $threads $fasta $pe1 $pe2 -f $file.$fa.$mapper.sam 2>> mapping.$mapper.log";} ## Appending STDERR to mapping.$mapper.log"
			elsif ($mapper eq 'bowtie2'){system "$bowtie2"."bowtie2 -p $threads -x $fa.bt2 -X $maxins -1 $pe1 -2 $pe2 -S $file.$fa.$mapper.sam 2>> mapping.$mapper.log";} 
			elsif ($mapper eq 'hisat2'){system "$hisat2"."hisat2 -p $threads --phred33 -x $fa.ht -1 $pe1 -2 $pe2 --no-spliced-alignment -S $file.$fa.$mapper.sam 2>> mapping.$mapper.log";}
			samtools(); variant(); stats(); logs(); ## Printing stats
			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $file.$fa.$mapper.bam";}
			unless ($sam) {system "rm $file.$fa.$mapper.sam";}
		}
	}
}

## Cleaning up
if ($mapper eq 'bwa'){system "mkdir $mapper.indexes; mv ${dir}*.amb ${dir}*.ann ${dir}*.bwt ${dir}*.fai ${dir}*.pac ${dir}*.sa $mapper.indexes/";}
elsif ($mapper eq 'bowtie2'){system "mkdir $mapper.indexes; mv *.bt2 ${dir}*.fai $mapper.indexes/";}
elsif ($mapper eq 'hisat2'){system "mkdir $mapper.indexes; mv *.ht2 ${dir}*.fai $mapper.indexes/";}
system "mkdir $mapper.$caller.VCFs; mv *.vcf $mapper.$caller.VCFs/";
system "mkdir $mapper.$caller.coverage; mv *.$mapper.coverage $mapper.$caller.coverage/";
system "mkdir $mapper.$caller.stats; mv *.$mapper.$type.stats $mapper.$caller.stats/";

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
	print "Running samtools on $file.$fa.$mapper.sam...\n";
	system "$samtools"."samtools view -@ $threads -bS $file.$fa.$mapper.sam -o $file.$fa.bam";
	system "$samtools"."samtools sort -@ $threads -o $file.$fa.$mapper.bam $file.$fa.bam";
	system "$samtools"."samtools depth -aa $file.$fa.$mapper.bam > $file.$fa.$mapper.coverage"; ## Printing per base coverage information
	system "rm $file.$fa.bam"; ## Discarding unsorted BAM file
}

sub variant{
	print "Calling variants with $caller on $fasta...\n\n";
	if ($caller eq 'varscan2'){
		if (($type eq 'snp')||($type eq 'indel')){system "$samtools"."samtools mpileup -f $fasta $file.$fa.$mapper.bam | java -jar $varjar mpileup2$type --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --output-vcf > $file.$fa.$mapper.$type.vcf";}
		elsif ($type eq 'both'){system "$samtools"."samtools mpileup -f $fasta $file.$fa.$mapper.bam | java -jar $varjar mpileup2cns --min-coverage $mc --min-reads2 $mr --min-avg-qual $maq --min-var-freq $mvf --min-freq-for-hom $mhom --p-value $pv --strand-filter $sf --variants --output-vcf > $file.$fa.$mapper.$type.vcf";}
		else {print "\nERROR: Unrecognized variant type. Please use: snp, indel, or both\n\n"; exit;}
	}
	elsif ($caller eq 'bcftools'){
		if ($type eq 'snp'){system "$samtools"."samtools mpileup -ugf $fasta $file.$fa.$mapper.bam | $bcftools"."bcftools call -vmO v -V indels --ploidy $ploidy --output $file.$fa.$mapper.$type.vcf";}
		elsif ($type eq 'indel'){system "$samtools"."samtools mpileup -ugf $fasta $file.$fa.$mapper.bam | $bcftools"."bcftools call -vmO v -V snps --ploidy $ploidy --output $file.$fa.$mapper.$type.vcf";}
		elsif ($type eq 'both'){system "$samtools"."samtools mpileup -ugf $fasta $file.$fa.$mapper.bam | $bcftools"."bcftools call -vmO v --ploidy $ploidy --output $file.$fa.$mapper.$type.vcf";}
		else {print "\nERROR: Unrecognized variant type. Please use: snp, indel, or both\n\n"; exit;}
	}
	elsif ($caller eq 'freebayes'){ ## single thread only, parallel version behaving wonky
		system "$samtools"."samtools index $file.$fa.$mapper.bam";
		system "$freebayes"."freebayes -f $fasta -p $ploidy $file.$fa.$mapper.bam > $file.$fa.$mapper.$caller.SNP.vcf";
		system "rm $file.$fa.$mapper.bam.bai";
	}
}

sub stats{
	print "\nCalculating stats...\n";
	open COV, "<$file.$fa.$mapper.coverage";
	open SN, "<$file.$fa.$mapper.$type.vcf";
	open STATS, ">$file.$fa.$mapper.$type.stats";
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
	if ($type eq 'both'){print STATS "Total number of SNPs + indels found: $sn\n"; print STATS "Average number of SNPs + indels per Kb: $snkb\n";}
	elsif ($type eq 'snp') {print STATS "Total number of SNPs found: $sn\n"; print STATS "Average number of SNPs per Kb: $snkb\n";}
	elsif ($type eq 'indel') {print STATS "Total number of indels found: $sn\n"; print STATS "Average number of indels per Kb: $snkb\n";}
	print STATS "\n## SAMTOOLS flagstat metrics\n";
	close STATS;
	system "$samtools"."samtools flagstat $file.$fa.$mapper.bam >> $file.$fa.$mapper.$type.stats";
}

sub logs{
	my $run_time = time - $tstart; my $mend = localtime();
	$comparison++;
	print LOG "Comparison # $comparison : $file (and mate, if PE) vs. $fasta - cumulative time elapsed: $run_time seconds\n";
	print MAP "\n".'###'." Mapping ended on $mend\n\n";
}

sub print_options{
	print MAP "\nOPTIONS:\n\n";
	print MAP "get_SNP.pl version: $version\n";
	print MAP "Number of threads: $threads\n";
	print MAP "Read mapper: $mapper\n";
	if ($mapper eq 'bwa'){print MAP "BWA alogrithm: $algo\n";}
	if ($mapper eq 'bowtie'){print MAP "Max insert size for bowtie: $maxins nt\n";}
	print MAP "Variant caller: $caller\n";
	if ($caller eq 'varscan2'){
		print MAP "   Varscan jar file used: $varjar\n";
		print MAP "   Minimum read depth at a position to make a call: $mc\n";
		print MAP "   Minimum supporting reads at a position to call variants: $mr\n";
		print MAP "   Minimum base quality at a position to count a read: $maq\n";
		print MAP "   Minimum variant allele frequency threshold: $mvf\n";
		print MAP "   Minimum frequency to call homozygote: $mhom\n";
		print MAP "   P-value threshold for calling variants: $pv\n";
		print MAP "   Strand filter: $sf\t\#\# 1 ignores variants with >90% support on one strand\n";
	}
	if ($caller eq "bcftools|freebayes"){print MAP "Setting ploidy to: $ploidy\n";}
	if ($type eq 'snp'){print MAP "Searching for SNPs...\n";}
	elsif ($type eq 'indel') {print MAP "Searching for indels...\n";}
	elsif ($type eq 'both') {print MAP "Searching for SNPs and indels...\n";}
}
