#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '1.1b';
my $name = 'run_Mash.pl';
my $updated = '2022-01-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs Mash on the provided fasta files
REQUIREMENTS	Mash - https://github.com/marbl/Mash (Ondov et al. DOI: 10.1186/s13059-016-0997-x)

USAGE		${name} \
		  -f DATASETS/*.fasta \
		  -o MASH \
		  -n S_pneumoniae_75.mash 

OPTIONS:
-f (--fasta)	Reference genome(s) in fasta file
-o (--outdir)	Output directory [Default: ./]
-n (--name)	Output file name [Default: Mash.mash]
-s (--sort)	Sort Mash output by decreasing order of similarity
OPTIONS
die "\n$usage\n" unless@ARGV;

my @fasta;
my $outdir = './';
my $filename = 'Mash.mash';
my $sort;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|outdir=s' => \$outdir,
	'n|name=s' => \$filename,
	's|sort' => \$sort,
);

## Checking if Mash is installed
my $prog = `echo \$(command -v mash)`;
chomp $prog;
if ($prog eq ''){
	die "\nERROR: Cannot find Mash. Please install Mash in your \$PATH\n\n";
}

## Creating output directory
unless (-d $outdir){
	make_path( $outdir, { mode => 0755 } ) or die "Can't created folder $outdir: $!\n";
}

## Running Mash
print "Running Mash genetic distance analysis...";
system "mash sketch -o ${outdir}/reference @fasta";
system "mash dist ${outdir}/reference.msh @fasta > ${outdir}/$filename";
if ($sort){
	print "Sorting out Mash results -- See ${outdir}/$filename.sorted";
	system "sort -gk3 ${outdir}/$filename > ${outdir}/$filename.sorted";
}
