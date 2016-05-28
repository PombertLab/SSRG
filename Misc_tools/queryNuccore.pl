#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2016
## Retrieve fasta files from the NCBI Nucleotide database using a list of accession numbers (one per line)

use LWP::Simple;
use strict;
use warnings;

my $usage = "USAGE = queryNuccore.pl Acession_list";
die "$usage\n" unless@ARGV;

my $start = localtime();
my $tstart = time;

open IN, "<$ARGV[0]";

my $url1 = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=';
my $url2 = '&rettype=fasta&retmode=text';

while (my $accession = <IN>){
	chomp $accession;
	if ($accession =~ /^#/){next;} ## Discarding comments
	else{
		my $filename = "$accession".'.fasta';
		my $URL = $url1.$accession.$url2;
		system "echo Downloading $accession.fasta ...";
		getstore($URL, $filename);
	}
}

my $time_taken = time - $tstart;
my $end = localtime();

print "\nDownloading started on: $start\n";
print "Downloading ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";
