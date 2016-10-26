#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2016
## Retrieve multifasta files from NCBI's FTP
## TAB/CSV-delimited lists can be generated from NCBI's Genome Assembly and Annotation reports

use strict;
use warnings;

my $usage = 'USAGE = queryNCBI.pl genome_list';
die "\n$usage\n
list =	TAB/CSV-delimited list generated from NCBI's Genome Assembly and Annotation reports, e.g.:
	1) goto http://www.ncbi.nlm.nih.gov/genome/genomes/159?
	2) click on the Download Table link in the upper right corner 
\n" unless@ARGV;

my $start = localtime();
my $tstart = time;

system "dos2unix $ARGV[0]"; ## making sure that the line breaks are in UNIX format
open IN, "<$ARGV[0]";

my $url = undef;
my $accession = undef;
my $suffix = '_genomic.fna.gz';

while (my $line = <IN>){
	chomp $line;
	$line =~ s/\"//g; ## Removing quotes
	if ($line =~ /^#/){next;} ## Discarding comments
	else {
		my @genome = split(/\t|,/, $line);
		## Organism info
		my $organism = $genome[0];
		my $genus = undef;
		my $species = undef;
		if ($organism =~ /^(\w+)\s+(\w+)/){$genus = $1; $species = $2;} ## Getting rid of random strain tags
		# Strain info
		my $strain = $genome[1]; $strain =~ s/ /_/g; $strain =~ s/\//_/g; $strain =~ s/\|/_/g; ## Substituting problematic characters with underscores
		# Accession info
		my $genbankFTP = $genome[$#genome];
		if ($genbankFTP =~ /.*\/(\S+)$/){$accession = $1;}
		my $filename = $accession."_genomic.fna.gz";
		## Downloading files
		$url = $genbankFTP.'/'.$accession.$suffix;
		system "wget $url"; ## Fetching compressed file
		system "mv $filename  $genus".'_'."$species".'_'."$strain".'.fasta.gz'; ## renaming file according to human readable species/strains name
		system "gunzip $genus".'_'."$species".'_'."$strain".'.fasta.gz'; ## decompressing file
	}
}
	
my $end = localtime();
my $time_taken = time -$tstart;

print "\nDownloads started on: $start\n";
print "Downloads completed on: $end\n";
print "Time elapsed: $time_taken seconds\n";