#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2019
my $version = '0.2';
my $name = 'queryWGS.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Retrieves multifasta files from the NCBI Sequence Set Browser
		Lists can be generated from http://www.ncbi.nlm.nih.gov/Traces/wgs/
		CSV lists must include all columns; from 1 (type) to 12 (date)
		
EXAMPLE		queryWGS.pl -l wgs_selector.csv -fa -gb

OPTIONS:
-l (--list)	CSV list from the NCBI Sequence Set Browser [Default: wgs_selector.csv]
-fa (--fasta)	Reference genome(s) in fasta format
-gb (--genbank)	Reference genome(s) in GenBank file format (GBFF; if available)
OPTIONS
die "$usage\n" unless @ARGV;

## Defining options
my $list = 'wgs_selector.csv';
my $fasta;
my $genbank;
GetOptions(
	'l|list=s' => \$list,
	'fa|fasta' => \$fasta,
	'gb|genbank' => \$genbank,
);

## Downloading data from NCBI
my $start = localtime(); my $tstart = time;
open LIST, "<$list" or die "Error. Can\'t open CSV list $list...\n"; 
my $taxa = 0;
my $ftp = 'ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/';
while (my $line = <LIST>){
	chomp $line;
	if ($line =~ /^#/){next;} ## Discarding comments
	elsif ($line =~ /^Prefix/i){next;}
	else {
		my @columns = split(",", $line);
		my $prefix = $columns[0];		## Prefix
		my $type = $columns[1];			## Type
		my $db = $columns[2];			## Source database
		my $locus = $columns[3];		## Targeted Locus Name
		my $div = $columns[4];			## Division; Bacteria (BCT), Eukaryotes...
		my $organism = $columns[5];		## Organism
		$organism =~ s/ /_/g;
		my $bioproject = $columns[6];	## Bioproject
		my $biosample = $columns[7];	## Biosample
		my $strain = $columns[8];		## Infraspecific Name
		my $source = $columns[9];		## Other Source
		my $contig_len = $columns[10];	## Contigs length sum
		my $contig_num = $columns[11];	## Contigs number
		my $contig_prot = $columns[12];	## Contigs proteins
		my $contig_anno = $columns[13];	## Contigs annotations
		my $scaff_len = $columns[14];	## Scaffolds length sum
		my $scaff_num = $columns[15];	## Scaffolds number
		my $scaff_prot = $columns[16];	## Scaffolds proteins
		my $scaff_anno = $columns[17];	## Scaffolds annotations
		my $update = $columns[19];		## Date - update
		my $date_creatd = $columns[19];	## Date - created
		my $length = $7;
		if ($db eq 'RefSeq'){next;} ## RefSeq data is redundant with GenBank/ENA entries; skipping
		else{
			$taxa++;
			my $path = '';
			my @prefixes = unpack ("(A2)*", $prefix);
			for (0..$#prefixes-1){$path .= $prefixes[$_]; $path .= '/';}
			$path .= $prefix; $path .= '/';
			if ($fasta){
				system "wget ${ftp}${path}${prefix}".'.1.fsa_nt.gz';
				system "mv $prefix.1.fsa_nt.gz $organism.$prefix.fasta.gz";
				system" gunzip $organism.$prefix.fasta.gz";
			}
			if ($genbank){
				system "wget ${ftp}${path}${prefix}".'.1.gbff.gz';
				system "mv $prefix.1.gbff.gz $organism.$prefix.gbff.gz";
				system" gunzip $organism.$prefix.gbff.gz";
			}
		}
	}
}

my $end = localtime();
my $time_taken = time -$tstart;

print "Genomes found: $taxa\n";

print "\nDownloads started on: $start\n";
print "Downloads completed on: $end\n";
print "Time elapsed: $time_taken seconds\n";