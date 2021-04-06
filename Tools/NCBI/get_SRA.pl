#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '0.5';
my $name = 'get_SRA.pl';
my $updated = '2021-03-15';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Downloads data from NCBI SRA files from the NCBI FTP site and converts to FASTQ format using fasterq-dump
		from the NCBI SRA toolkit
REQUIREMENTS	NCBI SRA toolkit:
		https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

EXAMPLE		${name} \\
		  -t 10 \\
		  -o FASTQ \\
		  -l accession_list(s) \\
		  -p \\
		  -v

OPTIONS:
-t (--threads)	Number of CPU threads to use [Default: 10]
-o (--outdir)	Output directory [Default: ./]
-l (--list)	List(s) of SRA accesssion numbers, one accession number per line
-p (--progess)	Show progess
-v (--verbose)	Add verbosity
-f (--force)	Force to overwrite existing file(s)
OPTIONS
die "\n$usage\n" unless @ARGV;

## Defining options
my $threads = 8;
my $outdir = './';
my @list;
my $verbose;
my $progess;
my $force;
GetOptions(
	't|threads=i' => \$threads,
	'o|outdir=s' => \$outdir,
	'l|list=s' => \@list,
	'v|verbose' => \$verbose,
	'p|progress' => \$progess,
	'f|force' => \$force
);

## Creating flags
my $verbosity = '';
my $overwrite = '';
my $pgress;
if ($verbose){ $verbosity = '--verbose'; }
if ($force){ $overwrite = '--force'; }
if ($progess){ $pgress = '--progress'; }

## Working on list file
while (my $list = shift@list){ 
	open IN, "<", "$list" or die "Can't read file $list: $!\n";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^\w+/){
			my $sra = $line;
			print "\nDownloading/converting $sra to FASTQ format with fasterq-dump ...\n\n";
			system "fasterq-dump \\
			  --threads $threads \\
			  $verbosity \\
			  $overwrite \\
			  $pgress \\
			  --outdir $outdir \\
			  --force \\
			  -3 $sra";
		}
	}
}