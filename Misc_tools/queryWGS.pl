#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2016
## Retrieve multifasta files from WGS accession number lists
## Lists can be generated from http://www.ncbi.nlm.nih.gov/Traces/wgs/

use strict;
use warnings;

my $usage = 'USAGE = queryWGS.pl WGS_list.txt';
die "\n$usage\n
list = TAB-delimited list generated from NCBI: http://www.ncbi.nlm.nih.gov/Traces/wgs/
" unless@ARGV;

my $start = localtime();
my $tstart = time;

open IN, "<$ARGV[0]";

my $taxa = 0;

while (my $line = <IN>){
	chomp $line;
	if ($line =~ /^#/){next;} ## Discarding comments
	elsif ($line =~ /^Prefix/){next;}
	elsif ($line =~ /^NZ_/){next;}
	elsif ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(.*)\t(\d+)\t(\d+)\t(no|yes)/){
		my $prefix = $1;
		my $div = $2;
		my $bioproject = $3;
		my $biosample = $4;
		my $organism = $5; $organism =~ s/ /_/g;
		my $contigs = $6;
		my $length = $7;
		$taxa++;
		system 'wget http://www.ncbi.nlm.nih.gov/Traces/wgs/?download='."$prefix".'.1.fsa_nt.gz';
		system "mv index.html\?download=$prefix.1.fsa_nt.gz $organism.$prefix.fasta.gz";
		system" gunzip $organism.$prefix.fasta.gz";
	}
}

my $end = localtime();
my $time_taken = time -$tstart;

print "Genomes found: $taxa\n";

print "\nDownloads started on: $start\n";
print "Downloads completed on: $end\n";
print "Time elapsed: $time_taken seconds\n";