#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '1.1';
my $name = 'run_Mash.pl';
my $updated = '14/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs Mash on the provided fasta files
REQUIREMENTS	Mash - https://github.com/marbl/Mash (Ondov et al. DOI: 10.1186/s13059-016-0997-x)

USAGE		${name} \
		  -f DATASETS/*.fasta \
		  -d MASH \
		  -o S_pneumoniae_75.mash 

OPTIONS:
-f (--fasta)		Reference genome(s) in fasta file
-d (--dir)		Output directory [Default: ./]
-o (--out)		Output file name [Default: Mash.mash]
-s (--sort)		Sort Mash output by decreasing order of similarity
OPTIONS
die "\n$usage\n" unless@ARGV;

my @fasta;
my $outdir = './';
my $out = 'Mash.mash';
my $sort;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'd|dir=s' => \$outdir,
	'o|out=s' => \$out,
	's|sort' => \$sort,
);

## Checking if Mash is installed
my $prog = `command -v mash`;
chomp $prog;
if ($prog eq ''){
	die "\nERROR: Cannot find Mash. Please install Mash in your \$PATH\n\n";
}

## Creating output directory
unless (-d $outdir){
	mkdir ($outdir,0755) or die "Can't created folder $outdir: $!\n";
}

## Running Mash
print "Running Mash genetic distance analysis...";
system "mash sketch -o reference @fasta";
system "mash dist reference.msh @fasta > ${outdir}/$out";
if ($sort){
	print "Sorting out Mash results -- See ${outdir}/$out.sorted";
	system "sort -gk3 ${outdir}/$out > ${outdir}/$out.sorted";
}
