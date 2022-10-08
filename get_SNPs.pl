#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2020
my $version = '2.0d';
my $name = 'get_SNPs.pl';
my $updated = '2022-10-08';

use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

## Usage definition
my $hint = "Type get_SNPs.pl -h (--help) for list of options\n";
my $usage = <<"USAGE";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Automates read-mapping of FASTQ datasets against reference sequence(s) and
		performs variant calling (if desired)
		
USAGE		${name} [options]

EXAMPLES:
simple			${name} -fa *.fasta -fq *.fastq -o RESULTS
advanced		${name} --fasta *.fasta --fastq *.fastq --mapper bowtie2 --caller varscan2 --type both --var ./VarScan.v2.4.3.jar --threads 16
paired ends		${name} --fasta *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --X 1000 --mapper bowtie2 --caller bcftools --threads 16
read-mapping only	${name} --fasta *.fasta --pe1 *R1.fastq --pe2 *R2.fastq --X 1000 --mapper bowtie2 --rmo --bam --threads 16
USAGE
die "\n$usage\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';
OPTIONS:
-h (--help)	Display this list of options
-v (--version)	Display script version
-o (--outdir)	Output directory [Default: ./]

# Mapping options
-fa (--fasta)			Reference genome(s) in fasta file
-fq (--fastq)			Fastq reads (single ends) to be mapped against reference(s)
-pe1				Fastq reads #1 (paired ends) to be mapped against reference(s)
-pe2				Fastq reads #2 (paired ends) to be mapped against reference(s)
-mapper				Read mapping tool: bowtie2, minimap2, winnowmap, ngmlr or hisat2 [default: minimap2]
-threads			Number of processing threads [default: 16]
-mem				Max total memory for samtools (in Gb) [default: 16] ## mem/threads = memory per thread
-bam				Keeps BAM files generated
-idx (--index)			Type of bam index generated (bai or csi) [default = csi] ## .bai not compatible with -mem
-sam				Keeps SAM files generated; SAM files can be quite large
-rmo (--read_mapping_only)	Do not perform variant calling; useful when only interested in bam/sam files and/or mapping stats
-ns (--no_stats)		Do not calculate stats; stats can take a while to compute for large eukaryote genomes

# Mapper-specific options
-preset				MINIMAP2/Winnowmap - Preset: sr, map-ont, map-pb or asm20 [default: sr]
-X				BOWTIE2 - Maximum paired ends insert size [default: 750]

# Variant calling options
-caller				[default: varscan2]	## Variant caller: varscan2 or bcftools
-type				[default: snp]		## snp, indel, or both
-ploidy				[default: 1]		## BCFtools option; change ploidy (if needed)

# VarScan2 parameters - see http://dkoboldt.github.io/varscan/using-varscan.html
-var				[default: VarScan.v2.4.4.jar]	## Which varscan jar file to use
-mc (--min-coverage)		[default: 15]		## Minimum read depth at a position to make a call
-mr (--min-reads2)		[default: 15]		## Minimum supporting reads at a position to call variants
-maq (--min-avg-qual)		[default: 28]		## Minimum base quality at a position to count a read
-mvf (--min-var-freq)		[default: 0.7]		## Minimum variant allele frequency threshold
-mhom (--min-freq-for-hom)	[default: 0.75]		## Minimum frequency to call homozygote
-pv (--p-value)			[default: 1e-02]	## P-value threshold for calling variants 
-sf (--strand-filter)		[default: 0]		## 0 or 1; 1 ignores variants with >90% support on one strand
END_OPTIONS
my @command = @ARGV;

my $help ='';
my $vn;
my $outdir = './';

## Mapping
my $mapper = 'minimap2';
my $threads = 16;
my $mem = 16;
my $index = 'csi';
my $sam = '';
my $bam = '';
my @fasta;
my @fastq;
my @pe1;
my @pe2;
my $rmo;
my $nostat;

## Mapper-specific
my $preset = 'sr';
my $maxins = '750';

## Variant calling
my $caller = 'varscan2';
my $type = 'snp';
my $ploidy = 1;

## VarScan-specific
my $varjar = 'VarScan.v2.4.4.jar';
my $mc = 15;
my $mr = 15;
my $maq = 28;
my $mvf = '0.7';
my $mhom = '0.75';
my $pv = '1e-02';
my $sf = 0;

GetOptions(
	'h|help' => \$help,
	'v|version', => \$vn,
	'o|outdir=s' => \$outdir,
	
	## Mapping
	'mapper=s' => \$mapper,
	'threads=i' => \$threads,
	'mem=i' => \$mem,
	'idx|index=s' => \$index,
	'sam' => \$sam,
	'bam' => \$bam,
	'fa|fasta=s@{1,}' => \@fasta,
	'fq|fastq=s@{1,}' => \@fastq,
	'pe1=s@{1,}' => \@pe1,
	'pe2=s@{1,}' => \@pe2,
	'rmo|read_mapping_only' => \$rmo,
	'ns|no_stats' => \$nostat,
	
	## Mapper-specific
	'preset=s' => \$preset,
	'X=i' => \$maxins,

	## Variant calling
	'caller=s' => \$caller,
	'type=s' => \$type,
	'ploidy=i' => \$ploidy,
	
	## VarScan-specific
	'var=s' => \$varjar,
	'mc|min-coverage=i' => \$mc,
	'mr|min-reads2=i' => \$mr,
	'maq|min-avg-qual=i' => \$maq,
	'mvf|min-var-freq=s' => \$mvf,
	'mhom|min-freq-for-hom=s' => \$mhom,
	'pv|p-value=s' => \$pv,
	'sf|strand-filter=i' => \$sf,
);

if ($help){die "\n$usage\n$options\n";}
if ($vn){die "\nversion $version\n\n";}

### Program checks for samtools, read mappers and variant callers
chkinstall('samtools');
if ($mapper =~ /bowtie2|hisat2|minimap2|winnowmap|ngmlr/){ chkinstall($mapper);}
else {die "\nMapper option $mapper is unrecognized. Please use bowtie2, hisat2, minimap2, winnowmap or ngmlr...\n\n";}
unless ($rmo){
	if ($caller =~ /bcftools/){ chkinstall($caller);}
	elsif ($caller eq 'varscan2'){ unless (-f $varjar){ die "Cannot find varscan jar file: $varjar\n";} }
	else {die "\nVariant caller option $caller is unrecognized. Please use varscan2 or bcftools...\n\n";}
}

## Option checks
if (scalar@fasta == 0) { die "\nPlease enter at least one FASTA reference before proceeding...\n\n"; }
if ((scalar@fastq == 0) && (scalar@pe1 == 0)) { die "\nPlease enter at least one FASTQ dataset before proceeding...\n\n"; }

## Setting the caller variable to rmo (read-mapping only); doesn't make sense to
## keep the default variant caller label for file names if -rmo is enabled...
if ($rmo){$caller = 'rmo';} 

## Creating output directory
unless (-d $outdir){
	make_path( $outdir, { mode => 0755 } ) or die "Can't create output directory $outdir: $!\n";
}

## Timestamps and logs
my $start = localtime();
my $tstart = time;
my $todo = 0;
if (@fastq){$todo += scalar(@fastq)*scalar(@fasta);}
if (@pe1 && @pe2){$todo += scalar(@pe1)*scalar(@fasta);}

my $map_log = "${outdir}/mapping.$mapper.log";
my $time_log = "${outdir}/time.$mapper.$caller.log";
open MAP, ">>", "$map_log" or die "Can't create file $map_log: $!\n";
open LOG, ">", "$time_log" or die "Can't create file $time_log.log: $!\n";

print MAP "COMMAND LINE:\n${name} @command\n";
print LOG "Mapping/SNP calling started on: $start\n";
print LOG "A total of $todo pairwise comparisons will be performed\n";
print_options();
my $comparison = 0;

##### Running read mapping/SNP calling
my $fasta; my $fastq; my $file; my $fa; my $dir; my $qdir; my $flagstat;

## Creating indexes and/or kmer frequencies 
foreach (@fasta){
	$fasta = $_;
	($fa, $dir) = fileparse($fasta);
	unless (-f $fasta) { 
		die "\nFASTA file named $fasta not found. Please check your command line...\n\n";
	}
	if ($mapper eq 'bowtie2'){
		system ("bowtie2-build --threads $threads $_ $fa.bt2") == 0 or checksig();
	}
	elsif ($mapper eq 'hisat2'){
		system ("hisat2-build $_ $fa.ht") == 0 or checksig();
	}
	elsif ($mapper eq 'winnowmap'){
		system ("meryl count k=15 output merylDB $fasta") == 0 or checksig();
		system ("meryl print greater-than distinct=0.9998 merylDB > $fa.repetitive_k15.txt") == 0 or checksig();
	}
}
my $index_time = time - $tstart;
print LOG "Time required to create all indexes: $index_time seconds\n";

## Read mapping files
my $log_file = "${outdir}/mapping.$mapper.log";
my $sam_file;
my $bam_file;

## Single end read-mapping
if (@fastq){
	foreach (@fastq){
		$fastq = $_;
		($file, $qdir) = fileparse($fastq);
		unless (-f $fastq) { 
			die "\nFASTQ file named $fastq not found. Please check your command line...\n\n";
		}
		print "\n## FASTQ information:\n";
		print "FASTQ parsed as: $file\n";
		print "FASTQ input DIR parsed as: $qdir\n";
		
		foreach (@fasta){
			$fasta = $_;
			($fa, $dir) = fileparse($fasta);
			$sam_file = "${outdir}/$file.$fa.$mapper.sam";
			$bam_file = "${outdir}/$file.$fa.$mapper.bam";

			print "\n## FASTA information:\n";
			print "FASTA parsed as: $fa\n";
			print "FASTA input DIR parsed as: $dir\n\n";
			print "Mapping $fastq on $fasta with $mapper...\n";
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$fastq vs. $fasta\n";

			## Read-mapping
			if ($mapper eq 'bowtie2'){
				system ("bowtie2 \\
				  --rg-id $fastq \\
				  --rg SM:$fasta \\
				  -p $threads \\
				  -x $fa.bt2 \\
				  -U $fastq \\
				  -S $sam_file \\
				  2>> $log_file") == 0 or checksig();
			} 
			elsif ($mapper eq 'minimap2'){ #-R \@RG\\\\tID:$fastq\\\\tSM:$fasta
				system ("minimap2 \\
				  -t $threads \\
				  --MD \\
				  -R \@RG\\\\tID:$fastq\\\\tSM:$fasta \\
				  -L \\
				  -ax $preset \\
				  $fasta \\
				  $fastq \\
				  1> $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}
			elsif ($mapper eq 'winnowmap'){
				system ("winnowmap \\
				  -t $threads \\
				  -W $fa.repetitive_k15.txt \\
				  --MD \\
				  -R \@RG\\\\tID:$fastq\\\\tSM:$fasta \\
				  -L \\
				  -ax $preset \\
				  $fasta \\
				  $fastq \\
				  1> $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}
			elsif ($mapper eq 'ngmlr'){
				system ("ngmlr \\
				  -t $threads \\
				  --rg-id $fastq \\
				  --rg-sm $fasta \\
				  -r $fasta \\
				  -q $fastq \\
				  -o $sam_file \\
				  2>&1 | tee $log_file") == 0 or checksig();
			}
			elsif ($mapper eq 'hisat2'){
				system ("hisat2 \\
				  -p $threads \\
				  --phred33 \\
				  --rg-id $fastq \\
				  --rg SM:$fasta \\
				  -x $fa.ht \\
				  -U $fastq \\
				  --no-spliced-alignment \\
				  -S $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}

			samtools(); ## Converting to BAM
			unless ($rmo){variant();} ## Calling variants
			unless ($nostat){stats();} ## Printing stats
			logs();

			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $bam_file $bam_file.$index";}
			unless ($sam) {system "rm $sam_file";}
		}
	}
}

## Paired ends read mapping
if (@pe1 && @pe2){
	while (my $pe1 = shift@pe1){
		my $pe2 = shift@pe2;
		unless (-f $pe1) { die "\nFASTQ PE1 file named $pe1 not found. Please check your command line...\n\n"; }
		unless (-f $pe2) { die "\nFASTQ PE2 file named $pe2 not found. Please check your command line...\n\n"; }
		($file, $qdir) = fileparse($pe1);
		my $r2 = fileparse($pe2);
		print "\n## FASTQ information:\n";
		print "R1 FASTQ parsed as: $file\n";
		print "R2 FASTQ parsed as: $r2\n";
		print "FASTQ input DIR parsed as: $qdir\n";

		foreach (@fasta){
			$fasta = $_;
			($fa, $dir) = fileparse($fasta);
			$sam_file = "${outdir}/$file.$fa.$mapper.sam";
			$bam_file = "${outdir}/$file.$fa.$mapper.bam";

			print "\n## FASTA information:\n";
			print "FASTA parsed as: $fa\n";
			print "FASTA input DIR parsed as: $dir\n\n";
			print "Mapping $pe1 and $pe2 on $fasta with $mapper...\n";
			my $mstart = localtime();
			print MAP "\n".'###'." Mapping started on $mstart\n";
			print MAP "\n$pe1 + $pe2 vs. $fasta\n";

			## Read mapping
			if ($mapper eq 'bowtie2'){
				system ("bowtie2 \\
				  --rg-id $pe1 \\
				  --rg SM:$fasta \\
				  -p $threads \\
				  -x $fa.bt2 \\
				  -X $maxins \\
				  -1 $pe1 \\
				  -2 $pe2 \\
				  -S $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}
			elsif ($mapper eq 'hisat2'){
				system ("hisat2 \\
				  -p $threads \\
				  --phred33 \\
				  --rg-id $pe1 \\
				  --rg SM:$fasta \\
				  -x $fa.ht \\
				  -1 $pe1 \\
				  -2 $pe2 \\
				  --no-spliced-alignment \\
				  -S $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}
			elsif ($mapper eq 'minimap2'){
				system ("minimap2 \\
				  -t $threads \\
				  -R \@RG\\\\tID:$pe1\\\\tSM:$fasta \\
				  -ax $preset \\
				  $fasta \\
				  $pe1 \\
				  $pe2 \\
				  1> $sam_file \\
				  2>> $log_file") == 0 or checksig();
			}

			samtools(); ## Converting to BAM 
			unless ($rmo){variant();} ## Calling variants
			unless ($nostat){stats();} ## Printing stats
			logs();
			
			## Cleaning SAM/BAM files
			unless ($bam) {system "rm $bam_file $bam_file.$index";}
			unless ($sam) {system "rm $sam_file";}
		}
	}
}

## Cleaning up
if ($mapper =~ /bowtie2|hisat2|ngmlr/){	mkdir ("${outdir}/$mapper.indexes",0755);}

if ($mapper eq 'bowtie2'){ system "mv ${outdir}/*.bt2 ${outdir}/*.fai ${outdir}/$mapper.indexes/"; }
elsif ($mapper eq 'hisat2'){ system "mv ${outdir}/*.ht2 ${outdir}/*.fai ${outdir}/$mapper.indexes/"; }
elsif ($mapper eq 'ngmlr'){ system "mv ${outdir}/*.ngm ${outdir}/$mapper.indexes/"; }

if ($bam){
	mkdir ("${outdir}/$mapper.BAM",0755);
	system "mv ${outdir}/*.bam ${outdir}/*.$index ${outdir}/$mapper.BAM/";
}

if ($sam){
	mkdir ("${outdir}/$mapper.SAM",0755);
	system "mv ${outdir}/*.sam ${outdir}/$mapper.SAM/";
}

unless ($rmo) {
	mkdir ("${outdir}/$mapper.$caller.VCFs",0755);
	system "mv ${outdir}/*.vcf ${outdir}/$mapper.$caller.VCFs/";
}

unless ($nostat){
	mkdir ("${outdir}/$mapper.$caller.depth",0755);
	mkdir ("${outdir}/$mapper.$caller.stats",0755);
	system "mv ${outdir}/*.$mapper.depth ${outdir}/$mapper.$caller.depth/";
	system "mv ${outdir}/*.$mapper.$type.stats ${outdir}/$mapper.$caller.stats/";
}

mkdir ("${outdir}/$mapper.$caller.coverage",0755);
system "mv ${outdir}/*.$mapper.coverage ${outdir}/$mapper.$caller.coverage/";

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


##### Subroutines #####

# sub to check if programs are installed
sub chkinstall {
	my $exe = $_[0];
	my $prog = `echo \$(command -v $exe)`;
	chomp $prog;
	if ($prog eq ''){die "\nERROR: Cannot find $exe. Please install $exe in your \$PATH\n\n";}
}

# sub to check for SIGINT and SIGTERM
sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	print "\nExit code = $exit_code; modulo = $modulo \n";

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n";
		exit(131);
	}

}

# sub to run convert SAM files to binary (BAM) format with samtools
sub samtools {
	my $coverage_file = "${outdir}/$file.$fa.$mapper.coverage";

	my $sammem = int(($mem/$threads)*1024);
	print "Running samtools on $sam_file...\n";
	print "Using $sammem Mb per thread for samtools\n";
	system "samtools view -@ $threads -bS $sam_file -o ${outdir}/tmp.bam";
	system "samtools sort -@ $threads -m ${sammem}M -o $bam_file ${outdir}/tmp.bam";
	if($index eq 'bai') { ## Can't use the memory -m CMD line switch with -b / .bai indexes
		print "samtools index (-b) not compatible with (-m); ";
		print "running samtools index with default memory settings...\n";
		system "samtools index -b -@ $threads $bam_file";
	}
	else { ## .csi indexes; work with -m
		system "samtools index -c -@ $threads -m ${sammem}M $bam_file";
	}
	system "samtools depth -aa $bam_file > $coverage_file"; ## Printing per base coverage information
	system "rm ${outdir}/tmp.bam"; ## Discarding unsorted temporary BAM file
	$flagstat = `samtools flagstat $bam_file`;
}

# sub to run the variant calling process
sub variant {
	(my $passQC) = ($flagstat =~ /(\d+)\s+\+\s+\d+ in total/s);
	print "\nQC-passed reads mapping to the genome = $passQC\n\n";

	my $vcf_file = "${outdir}/$file.$fa.$mapper.$type.vcf";

	## Checking if BAM file is empty
	if ($passQC == 0){
		print "No QC-passed reads mapped to the genome; skipping variants calling\n\n";
		open VCF, ">", "$vcf_file" or die "Can't create file $vcf_file: $!\n";
		print VCF '## No QC-passed reads mapped to the genome'."\n";
		close VCF;
	}

	## If not empty, proceed; empty BAM files create isssues when piping mpileup
	## to some variant callers (e.g. VarScan2).
	else{
		print "Calling variants with $caller on $fasta...\n\n";
		if ($caller eq 'varscan2'){
			my $vflag = '';
			my $cns;
			if (($type eq 'snp')||($type eq 'indel')){ $cns = $type; }
			elsif ($type eq 'both'){
				$cns = 'cns';
				$vflag = '--variants';
			}
			else {die "\nERROR: Unrecognized variant type. Please use: snp, indel, or both\n\n";}
			system "samtools \\
			  mpileup \\
			  -f $fasta \\
			  $bam_file \\
			  | \\
			  java -jar $varjar \\
			  mpileup2$cns \\
			  --min-coverage $mc \\
			  --min-reads2 $mr \\
			  --min-avg-qual $maq \\
			  --min-var-freq $mvf \\
			  --min-freq-for-hom $mhom \\
			  --p-value $pv \\
			  --strand-filter $sf \\
			  $vflag \\
			  --output-vcf \\
			  > $vcf_file";
		}
		elsif ($caller eq 'bcftools'){ ## tested with 1.10.2
			unless ($type =~ /snp|indel|both/){ die "\nERROR: Unrecognized variant type. Please use: snp, indel, or both\n\n"; }
			my $varV = '';
			if ($type eq 'snp'){ $varV = '-V indels'; }
			elsif ($type eq 'indel') { $varV = '-V snps'; }
			system "bcftools \\
			  mpileup \\
			  -f $fasta \\
			  $bam_file \\
			  | \\
			  bcftools \\
			  call \\
			  -vmO v \\
			  $varV \\
			  --ploidy $ploidy \\
			  --output $vcf_file";
		}
	}
}

## Sub to calculate read mapping stats/metrics
sub stats {
	my $run_time = time - $tstart; my $mend = localtime();

	my $coverage_file = "${outdir}/$file.$fa.$mapper.coverage";
	my $stats_file = "${outdir}/$file.$fa.$mapper.$type.stats";
	my $depth_file = "${outdir}/$file.$fa.$mapper.depth";
	my $vcf_file = "${outdir}/$file.$fa.$mapper.$type.vcf";

	open COV, "<", "$coverage_file" or die "Can't read file $coverage_file: $!\n";
	open STATS, ">", "$stats_file" or die "Can't create file $stats_file: $!\n";
	open DEPTH, ">", "$depth_file" or die "Can't create file $depth_file: $!\n";
	unless ($rmo){open SN, "<", "$vcf_file" or die "Can't read file $vcf_file.vcf: $!\n";}

	print "\nCalculating stats...\n";
	print DEPTH "Contig\tAverage depth\tAverage (all contigs)\tDifference\tRatio Contig\/Average\n";

	my $total = 0; my $covered = 0; my $nocov = 0; my $max = 0; my $sumcov; my $sn =0;
	
	## Run only if variant calling has been performed
	unless ($rmo){
		while (my $line = <SN>){
			if ($line =~ /^#/){next;}
			else {$sn++;}
		}
	}

	## Checks for sequencing depth
	my %data = (); my @contigs = ();
	while (my $line = <COV>){
		chomp $line;
		$total++; 
		if ($line =~ /^(\S+)\s+(\d+)\s+(\d+)/){
			my $contig = $1; my $position = $2; my $coverage = $3;
			$sumcov += $coverage;
			if ($coverage >= 1) {
				$covered++;
				if ($coverage > $max){$max = $coverage;}
			}
			else {$nocov++;}
			if (exists $data{$contig}){
				$data{$contig}[0] += 1; ## Keeping track of contig size 
				$data{$contig}[1] += $coverage; ## Keeping track of cumulative sequencing depth 
			}
			else{
				$data{$contig}[0] = 1; ## Initializing new contig 
				$data{$contig}[1] += $coverage; ## Keeping track of cumulative sequencing depth 
				push (@contigs, $contig);
			}
		}
	}

	## In case the output of samtools depth -aa generates a blank file (if no read maps to the reference)
	my $avc; if ($total > 0){$avc = sprintf("%.2f", ($sumcov/$total));}	else{$avc = 0;}

	## Printing sequencing depths per contig
	while (my $tmp = shift@contigs){
		my $average = ($data{$tmp}[1]/$data{$tmp}[0]);
		$average = sprintf("%.2f", $average);
		my $diff = $average - $avc; $diff = sprintf("%.2f", $diff);
		my $ratio;

		## Preventing possible division by zero
		if ($avc > 0){$ratio = $average/$avc; $ratio = sprintf("%.2f", $ratio);}
		else{$ratio = 0;}
		
		print DEPTH "$tmp\t$average\t$avc\t$diff\t$ratio\n";
	}
	print STATS "FASTQ file(s) used: $file (and mate, if PE)\n";
	print STATS "Reference fasta file used: $fasta\n";
	
	if ($total > 0){
		print STATS "\nTotal number of bases in the reference genome\t$total bp\n";
		print STATS "Number of bases covered by at least one read\t$covered\n";
		print STATS "Number of bases without coverage\t$nocov\n";
		print STATS "Maximum sequencing depth\t$max"."X\n";
		print STATS "Average sequencing depth\t$avc"."X\n";
		if ($total == $covered){print STATS "Sequencing breadth (percentage of bases covered by at least one read)\t100%\n";}
		if ($total != $covered){my $av_cov = sprintf("%.2f%%", ($covered/$total)*100);
		print STATS "Sequencing breadth (percentage of bases covered by at least one read)\t$av_cov\n";}
	}
	else{
		print STATS "\nNo read was found to map to the reference: ";
		print STATS "the output of samtools depth -aa ($coverage_file) is blank.\n";
	}
	
	## Execute only if variant calling has been performed
	unless ($rmo){
		my $snkb = sprintf("%.2f", ($sn/$covered)*1000);
		if ($type eq 'both'){
			print STATS "Total number of SNPs + indels found: $sn\n";
			print STATS "Average number of SNPs + indels per Kb: $snkb\n";
		}
		elsif ($type eq 'snp') {
			print STATS "Total number of SNPs found: $sn\n";
			print STATS "Average number of SNPs per Kb: $snkb\n";
		}
		elsif ($type eq 'indel') {
			print STATS "Total number of indels found: $sn\n";
			print STATS "Average number of indels per Kb: $snkb\n";
		}
	}
	print STATS "\n## SAMTOOLS flagstat metrics\n";
	print STATS "$flagstat\n";
	close STATS;
	print "Time to calculate stats: $run_time seconds\n";
}

## LOG subroutines
sub logs {
	my $run_time = time - $tstart; my $mend = localtime();
	$comparison++;
	print LOG "Comparison # $comparison : $file (and mate, if PE) vs. $fasta - cumulative time elapsed: $run_time seconds\n";
	print MAP "\n".'###'." Mapping ended on $mend\n\n";
}

sub print_options {
	print MAP "\nOPTIONS:\n\n";
	print MAP "get_SNP.pl version: $version\n";
	print MAP "Number of threads: $threads\n";
	print MAP "Read mapper: $mapper\n";
	if ($mapper eq 'bowtie'){print MAP "Max insert size for bowtie: $maxins nt\n";}
	if ($rmo){print MAP "Read-mapping only requested, skipping variant calling.\n";}
	if ($bam){print MAP "Keeping BAM files are requested.\n";}
	if ($sam){print MAP "Keeping SAM files are requested.\n";}
	if ($nostat){print MAP "Skipping stats calculations.\n";}
	unless ($rmo){
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
		if ($caller eq "bcftools"){print MAP "Setting ploidy to: $ploidy\n";}
		if ($type eq 'snp'){print MAP "Searching for SNPs...\n";}
		elsif ($type eq 'indel') {print MAP "Searching for indels...\n";}
		elsif ($type eq 'both') {print MAP "Searching for SNPs and indels...\n";}
	}
}
