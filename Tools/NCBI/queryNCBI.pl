#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2020
my $version = '1.3';
my $name = 'queryNCBI.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Retrieve multifasta files from NCBI's FTP
		TAB/CSV-delimited lists can be generated from NCBI's Genome Assembly and Annotation reports, e.g.:
		1) goto http://www.ncbi.nlm.nih.gov/genome/genomes/159?
		2) click on the Download Table link in the upper right corner 
		
EXAMPLE		queryNCBI.pl -l genome_list.csv -fa -gb -p

OPTIONS:
-l (--list)	TAB/CSV-delimited list from NCBI
-fa (--fasta)	Retrieve fasta files
-gb (--genbank)	Retrieve GenBank annotation files (if available)
-p (--protein)	Retrieve protein sequences (if available)
-cds		Retrieve protein coding sequences (if available)
OPTIONS
die "$usage\n" unless @ARGV;

## Defining options
my $fasta;
my $gbk;
my $list;
my $protein;
my $cds;
GetOptions(
	'fa|fasta' => \$fasta,
	'gb|genbank' => \$gbk,
	'p|protein' => \$protein,
	'cds' => \$cds,
	'l|list=s' => \$list
);

## Downloading data from NCBI
my $start = localtime(); my $tstart = time;
system "dos2unix $list"; ## making sure that the line breaks are in UNIX format
open IN, "<$list";
my $url = undef;
my $accession = undef;
my $fa = '_genomic.fna.gz';
my $gb = '_genomic.gbff.gz';
my $prot = '_protein.faa.gz';
my $cd = '_cds_from_genomic.fna.gz';

while (my $line = <IN>){
	chomp $line;
	$line =~ s/\"//g; ## Removing quotes
	if ($line =~ /^#/){next;} ## Discarding comments
	else {
		my @genome = split(/\t|,/, $line);
		## Organism info
		my $organism = $genome[0]; my $genus = undef; my $species = undef;
		if ($organism =~ /^(\S+)\s+(\S+)/){$genus = $1; $species = $2;} ## Getting rid of random strain tags
		# Strain info
		my $strain = $genome[2]; $strain =~ s/ /_/g; $strain =~ s/\//_/g; $strain =~ s/\|/_/g; $strain =~ s/\./_/g; $strain =~ s/\(//g; $strain =~ s/\)//g; ## Deleting or substituting problematic characters with underscores
		# Accession info
		my $genbankFTP = $genome[$#genome]; if ($genbankFTP =~ /.*\/(\S+)$/){$accession = $1;}
		## Downloading files
		if ($fasta){
			my $filename = $accession."$fa";
			$url = $genbankFTP.'/'.$accession.$fa;
			system "wget $url"; ## Fetching compressed file
			system "mv $filename  $genus".'_'."$species".'_'."$strain".'.fasta.gz'; ## renaming file according to human readable species/strains name
			system "gunzip $genus".'_'."$species".'_'."$strain".'.fasta.gz'; ## decompressing file
		}
		if ($gbk){
			my $filename = $accession."$gb";
			$url = $genbankFTP.'/'.$accession.$gb;
			system "wget $url"; ## Fetching compressed file
			system "mv $filename  $genus".'_'."$species".'_'."$strain".'.gbk.gz'; ## renaming file according to human readable species/strains name
			system "gunzip $genus".'_'."$species".'_'."$strain".'.gbk.gz'; ## decompressing file
		}
		if ($protein){
			my $filename = $accession."$prot";
			$url = $genbankFTP.'/'.$accession.$prot;
			system "wget $url"; ## Fetching compressed file
			system "mv $filename  $genus".'_'."$species".'_'."$strain".'.faa.gz'; ## renaming file according to human readable species/strains name
			system "gunzip $genus".'_'."$species".'_'."$strain".'.faa.gz'; ## decompressing file
		}
		if ($cds){
			my $filename = $accession."$cd";
			$url = $genbankFTP.'/'.$accession.$cd;
			system "wget $url"; ## Fetching compressed file
			system "mv $filename  $genus".'_'."$species".'_'."$strain".'.cds.gz'; ## renaming file according to human readable species/strains name
			system "gunzip $genus".'_'."$species".'_'."$strain".'.cds.gz'; ## decompressing file
		}
	}
}
my $end = localtime(); my $time_taken = time - $tstart;
print "\nDownloads started on: $start\n";
print "Downloads completed on: $end\n";
print "Time elapsed: $time_taken seconds\n";
