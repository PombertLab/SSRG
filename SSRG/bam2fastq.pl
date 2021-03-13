#!/usr/bin/perl
## Pombert Lab, IIT, 2019
my $name = 'bam2fastq.pl';
my $version = '0.2a';
my $updated = '13/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Extract sequencing reads in FASTQ format from BAM alignment files
REQUIREMENTS	Samtools 1.3.1+

USAGE		${name} \\
		  -b file.bam \\
		  -t pe \\
		  -e map \\
		  -p reads \\
		  -s fastq

OPTIONS:
-b (--bam)	BAM alignment file
-t (--type)	Alignment type: pe (paired ends) or single [Default: pe]
-e (--extract)	Reads to extract: map, unmap [Default: map]
-p (--prefix)	Output file(s) prefix [Default: reads]
-s (--suffix)	Output file(s) suffix [Default: fastq]
OPTIONS
die "\n$usage\n" unless @ARGV;

my $bam;
my $type = 'pe';
my $extract = 'map';
my $prefix = 'reads';
my $suffix = 'fastq';
GetOptions(
	'b|bam=s' => \$bam,
	't|type=s'	=> \$type,
	'e|extract=s' => \$extract,
	'p|prefix=s' => \$prefix,
	's|suffix=s' => \$suffix
);

## Program + option check
my $samtools = `command -v samtools`;
chomp $samtools;
if ($samtools eq ''){
	die "\nERROR: Cannot find Samtools. Please install Samtools in your path\n\n";
}
unless (($type eq 'pe') || ($type eq 'se')) {
	die "\nUnrecognized type $type. Please use 'pe' for paired-ends or 'se' for single ends\n";
}
unless (($extract eq 'map') || ($extract eq 'unmap')) {
	die "\nUnrecognized reads to extact: $extract. Please enter 'map' or 'unmap' to extract reads that map or do not map to the reference(s), respectively\n";
}

## Extracting reads from BAM files with samtools bam2fq
if ($type eq 'se'){ ## Single ends
	if ($extract eq 'map'){
		print "\nExtracting single reads mapping to reference(s) [flag: -F 4] from $bam\nSaving to $prefix.$suffix\n\n";
		system "samtools bam2fq -F 4 $bam > $prefix.$suffix";
	}
	elsif ($extract eq 'unmap'){
		print "\nExtracting single reads not mapping to reference(s) [flag: -f 4] from $bam\nSaving $prefix.$suffix\n\n";
		system "samtools bam2fq -f 4 $bam > $prefix.$suffix";
	}
}
elsif ($type eq 'pe'){ ## Paired ends (illumina)
	if ($extract eq 'map'){ ## R1 + R2 mapped
		print "\nExtracting paired-end reads mapping to reference(s) [flags: -f 1 -F 12] from $bam\nSaving to ${prefix}_R1.$suffix and ${prefix}_R2.$suffix\n\n";
		system "samtools bam2fq -f 1 -F 12 -1 ${prefix}_R1.$suffix -2 ${prefix}_R2.$suffix $bam";
	}
	elsif ($extract eq 'unmap'){ ## R1 + R2 didn’t map
		print "\nExtracting paired-end reads not mapping to reference(s) [flags: -f 12 -F 256] from $bam\nSaving to ${prefix}_R1.$suffix and ${prefix}_R2.$suffix\n\n";
		system "samtools bam2fq -f 12 -F 256 -1 ${prefix}_R1.$suffix -2 ${prefix}_R2.$suffix $bam";	
	}
}
