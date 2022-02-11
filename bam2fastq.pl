#!/usr/bin/perl
## Pombert Lab, IIT, 2019
my $name = 'bam2fastq.pl';
my $version = '0.4';
my $updated = '2022-02-11';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Extract sequencing reads in FASTQ format from BAM alignment files
REQUIREMENTS	Samtools 1.3.1+

USAGE 1		${name} -b file.bam -o FASTQ -t pe -auto map -p reads -s fastq
USAGE 2		${name} -b file.bam -o FASTQ -t pe -f 1 -g 12 -p reads -s fastq

OPTIONS:
-b (--bam)	BAM alignment file

## Output options
-o (--outdir)	Output directory [Default: ./]
-t (--type)	Alignment type: pe (paired ends) or single [Default: pe]
-p (--prefix)	Output file(s) prefix ## Defaults to bam file name if not set
-s (--suffix)	Output file(s) suffix [Default: fastq]

## SAM/BAM flag options
-a (--auto)	Extract reads automatically: all, map or unmap
-f SAM/BAM (-f) flag # https://www.samformat.info/sam-format-flag
-g SAM/BAM (-F) flag # https://www.samformat.info/sam-format-flag
OPTIONS
die "\n$usage\n" unless @ARGV;
my @command_line = @ARGV;

my $bam;
my $outdir = './';
my $type = 'pe';
my $auto;
my $F1;
my $F2;
my $prefix;
my $suffix = 'fastq';
GetOptions(
	'b|bam=s' => \$bam,
	'o|outdir=s' => \$outdir,
	't|type=s'	=> \$type,
	'a|auto=s' => \$auto,
	'f=i' => \$F1,
	'g=i' => \$F2,
	'p|prefix=s' => \$prefix,
	's|suffix=s' => \$suffix
);

## Program + option check
my $samtools = `echo \$(command -v samtools)`;
chomp $samtools;
if ($samtools eq ''){
	die "\nERROR: Cannot find Samtools. Please install Samtools in your path\n\n";
}
unless (($type eq 'pe') || ($type eq 'se')) {
	die "\nUnrecognized type $type. Please use 'pe' for paired-ends or 'se' for single ends\n";
}
if ($auto){
	$auto = lc($auto);
	unless (($auto eq 'all') || ($auto eq 'map') || ($auto eq 'unmap')) {
		die "\nUnrecognized mode: $auto. Please enter 'all', 'map' or 'unmap' to extract reads that map or do not map to the reference(s), respectively\n";
	}
}

## Creating output directory
unless (-d $outdir){
	make_path( $outdir, { mode => 0755 } ) or die "Can't create folder $outdir: $!\n";
}

## Creating log file
my $log_file = "${outdir}/b2fq.log";
my $date = localtime();
my $time_start = time;
open LOG, ">", "$log_file" or die "Can't create log file $log_file: $!\n";
print LOG "Command executed on $date:\n\n";
print LOG "${name} @command_line\n";

## Autopopulating prefix from basename if not defined from command line
unless ($prefix){
	my $basename = fileparse($bam);
	$basename =~ s/\.\w+$//;
	$prefix = $basename;
}

## Extracting reads from BAM files with samtools bam2fq
my $flags;

if ($type eq 'se'){ ## Single ends

	my $filename_SE = "${outdir}/$prefix.$suffix";

	## Setting flags, if any
	if ($auto){
		if ($auto eq 'all'){ $flags = ''; } ## All reads, no flag needed
		elsif ($auto eq 'map'){ $flags = '-F 4'; }
		elsif ($auto eq 'unmap'){ $flags = '-f 4';} 
	}
	elsif ($F1){ $flags = "-f $F1";}
	elsif ($F2){ $flags = "-F $F2";}

	if ($auto){
		if ($auto eq 'all'){
			print "\nExtracting all single reads from $bam\n";
		}
	}
	else {
		print "\nExtracting single reads [flag: $flags] from $bam\n";
	}

	## Saving to file
	print "Saving to $filename_SE\n\n";
	system "samtools \\
		bam2fq \\
		$flags \\
		$bam \\
		> $filename_SE";
}

elsif ($type eq 'pe'){ ## Paired ends (illumina)

	my $filename_R1 = "${outdir}/${prefix}_R1.$suffix";
	my $filename_R2 = "${outdir}/${prefix}_R2.$suffix";

	## Setting flags, if any
	if ($auto){
		if ($auto eq 'all'){ $flags = ''; } ## All reads, no flag needed
		elsif ($auto eq 'map'){ $flags = '-f 1 -F 12'; } ## R1 + R2 mapped
		elsif ($auto eq 'unmap'){  $flags = '-f 12 -F 256'; } ## R1 + R2 didn’t map
	}
	else {
		my $f1 = "-f $F1";
		my $f2 = "-F $F2";
		$flags = "$f1".' '."$f2";
	}

	if ($auto){
		if ($auto eq 'all'){
			print "\nExtracting all paired reads from $bam\n";
		}
	}
	else {
		print "\nExtracting paired-end reads [flags: $flags] from $bam\n";
	}

	## Saving to file
	print "Saving to $filename_R1 and $filename_R2\n\n";
	system "samtools \\
		bam2fq \\
		$flags \\
		-1 $filename_R1 \\
		-2 $filename_R2 \\
		$bam";
}

my $time_end = time;
my $run_time = $time_end - $time_start;
print LOG "\nRuntime: $run_time seconds\n";