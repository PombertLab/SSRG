#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2016
## Downloads SRA files from the NCBI FTP site 

use strict;
use warnings;
use Getopt::Long qw(GetOptions);  

## Environment variables. Leave blank if set in $PATH or replace with your settings
my $sratoolkit = '/opt/sratoolkit.2.6.3/bin/'; ## Path to NCBI SRA toolkit binaries

my $usage = '
USAGE = get_SRA.pl [options]

EXAMPLE = get_SRA.pl -l accession_list(s) -fq pe -Q 33 ## -fq converts SRA to fastq format

Options:
-l (--list)	List(s) of SRA accesssion numbers, one accession number per line
-fq (--fastq)	Converts SRA files to FASTQ paired-end (pe) or single (se) format ## Requires NCBI SRA toolkit; non-pe files will default to se
-Q		FASTQ quality score format; 33 (Sanger) or 64 (illumina) [default: 33]
';
die "$usage\n" unless @ARGV;

my @list;
my $fastq = '';
my $qscore = 33;
GetOptions(
'l|list=s' => \@list,
'fq|fastq=s' => \$fastq,
'Q=i' => \$qscore,
);

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