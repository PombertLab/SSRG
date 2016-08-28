#!/usr/bin/perl
## Pombert JF, Illinois Tech - 2016

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Path to Mash; leave blank (e.g. $mash = '';) if set in $PATH
my $mash = '/opt/Mash/bin/';

## Defining options
my $usage = "
Mash - https://github.com/marbl/Mash (Ondov et al. DOI: 10.1186/s13059-016-0997-x)

USAGE = perl run_Mash.pl [options]
EXAMPLE = run_Mash.pl --fasta *.fasta --out Mash.txt --sort 
OPTIONS:
--fasta		Reference genome(s) in fasta file
--out		Output file name
--sort		Sort Mash output by decreasing order of similarity
";

die "$usage\n" unless@ARGV;

my @fasta;
my $out = '';
my $sort = '';

GetOptions(
	'fasta=s@{1,}' => \@fasta,
	'out=s' => \$out,
	'sort' => \$sort,
);

system "echo Running Mash genetic distance analysis...";
system "$mash"."mash sketch -o reference @fasta";
system "$mash"."mash dist reference.msh @fasta > $out";
if ($sort){
	system "echo Sorting out Mash results -- See $out.sorted";
	system "sort -gk3 $out > $out.sorted";
}
