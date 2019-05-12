#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2019
my $version = '1.0';
my $name = 'run_Mash.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);
my $mash = '';	## Path to Mash; leave blank (e.g. $mash = '';) if set in $PATH

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Runs Mash on the provided fasta files
REQUIREMENTS	Mash - https://github.com/marbl/Mash (Ondov et al. DOI: 10.1186/s13059-016-0997-x)

USAGE		run_Mash.pl -f *.fasta -o Mash.txt 

OPTIONS:
-f (--fasta)		Reference genome(s) in fasta file
-o (--out)		Output file name [Default: Mash.txt]
-s (--sort)		Sort Mash output by decreasing order of similarity
OPTIONS
die "$options\n" unless@ARGV;

my @fasta;
my $out = 'Mash.txt';
my $sort;
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|out=s' => \$out,
	's|sort' => \$sort,
);

## Checking if Mash is installed
my $prog = `command -v ${mash}mash`;
chomp $prog;
if ($prog eq ''){
	print "\nERROR: Cannot find Mash. Please install Mash in your \$PATH or include its installation directory at the top of run_Mash.pl\n\n";
	exit;
}

## Running Mash
system "echo Running Mash genetic distance analysis...";
system "$mash"."mash sketch -o reference @fasta";
system "$mash"."mash dist reference.msh @fasta > $out";
if ($sort){
	system "echo Sorting out Mash results -- See $out.sorted";
	system "sort -gk3 $out > $out.sorted";
}
