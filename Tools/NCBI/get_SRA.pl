#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '0.3; No longer works, needs an overhaul';
my $name = 'get_SRA.pl';
my $updated = '15/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Downloads data from NCBI SRA files from the NCBI FTP site and converts to FASTQ format
REQUIREMENTS	NCBI SRA toolkit:
		https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software

EXAMPLE		${name} \\
		  -l accession_list(s) \\
		  -fq pe \\
		  -Q 33

OPTIONS:
-l (--list)	List(s) of SRA accesssion numbers, one accession number per line
-fq (--fastq)	Converts SRA files to FASTQ paired-end (pe) or single (se) format [default: pe] ## Non-pe files will default to se
-Q		FASTQ quality score format; 33 (Sanger) or 64 (illumina) [default: 33]
OPTIONS
die "\n$usage\n" unless @ARGV;

## Defining options
my @list;
my $fastq = 'pe';
my $qscore = 33;
GetOptions(
	'l|list=s' => \@list,
	'fq|fastq=s' => \$fastq,
	'Q=i' => \$qscore
);

## Downloading/converting SRA data
my $ftp = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/';
my $prefix = undef;
my $run = undef;
my $id = undef;
my $sra = undef;
while (my $list = shift@list){ 
	open IN, "<", "$list";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^(\w{3})(\w{3})(\w{3,})/){
			$prefix = $1;
			$run = $2;
			$id = $3;
			$sra = "$line.sra";
			my $url = "$ftp$prefix".'/'."$prefix$run".'/'."$line".'/'."$sra";
			print "Downloading $sra...";
			system "wget $url";
		}
		if ($fastq eq 'pe'){
			print "Converting $sra to FASTQ format...";
			system "fastq-dump -Q $qscore --split-3 $sra";
		}
		elsif ($fastq eq 'se'){
			print "Converting $sra to FASTQ format...";
			system "fastq-dump -Q $qscore $sra";
		}
	}
}