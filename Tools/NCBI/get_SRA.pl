#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
## Downloads SRA files from the NCBI FTP site 
my $version = '0.2';
my $name = 'get_SRA.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);  

## Usage definition
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Downloads data from NCBI SRA and converts to FASTQ format
REQUIREMENTS	NCBI SRA toolkit - https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
EXAMPLE		get_SRA.pl -l accession_list(s) -fq pe -Q 33 -sra /opt/sratoolkit.2.9.6/bin/

OPTIONS:
-l (--list)	List(s) of SRA accesssion numbers, one accession number per line
-fq (--fastq)	Converts SRA files to FASTQ paired-end (pe) or single (se) format [default: pe] ## Non-pe files will default to se
-Q		FASTQ quality score format; 33 (Sanger) or 64 (illumina) [default: 33]
-sra		Path to NCBI sratoolkit binaries
OPTIONS
die "$usage\n" unless @ARGV;

## Defining options
my @list;
my $fastq = 'pe';
my $qscore = 33;
my $sratoolkit = ''; ## Path to NCBI SRA toolkit binaries
GetOptions(
	'l|list=s' => \@list,
	'fq|fastq=s' => \$fastq,
	'Q=i' => \$qscore,
	'sra=s' => \$sratoolkit
);

## Downloading/converting SRA data
my $ftp = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/';
my $prefix = undef;
my $run = undef;
my $id = undef;
my $sra = undef;
while (my $list = shift@list){ 
	open IN, "<$list";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^(\w{3})(\w{3})(\w{3,})/){
			$prefix = $1;
			$run = $2;
			$id = $3;
			$sra = "$line.sra";
			my $url = "$ftp$prefix".'/'."$prefix$run".'/'."$line".'/'."$sra";
			system "echo Downloading $sra...";
			system "wget $url";
		}
		if ($fastq eq 'pe'){
			system "echo Converting $sra to FASTQ format...";
			system "$sratoolkit"."fastq-dump -Q $qscore --split-3 $sra";
		}
		elsif ($fastq eq 'se'){
			system "echo Converting $sra to FASTQ format...";
			system "$sratoolkit"."fastq-dump -Q $qscore $sra";
		}
	}
}