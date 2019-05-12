#!/usr/bin/perl
## Pombert lab, 2019
my $version = '1.1d';
my $name = 'MashToDistanceMatrix.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Converts the Mash output to a CSV or TSV matrix for easy plotting with R
USAGE		MashToDistanceMatrix.pl -i Mash.txt -o Mash -t csv

OPTIONS:
-i (--input)		Mash output file [Default: Mash.txt]
-o (--out)		Output file prefix [Default: Mash]
-f (--format)	Output format; csv or tsv [Default: csv]
OPTIONS
die "$options\n" unless @ARGV;

my $mash = 'Mash.txt';
my $out = 'Mash';
my $format = 'csv';
GetOptions(
	'i|input=s' => \$mash,
	'o|out=s' => \$out,
	'f|format=s' => \$format
);

## Reformatting Mash output to a CSV or TSV matrix
$format = lc($format);
my $sep;
if ($format eq 'csv'){$sep = ",";}
elsif ($format eq 'tsv'){$sep = "\t";}
else {die "Format $format in not recognized. Please specify csv or tsv...\n";}
open IN, "<$mash";
open MASH, ">$out.$format"; open CORR, ">$out.alternatedist.$format";
print MASH "OTU"; print CORR "OTU";
my %hash = (); my %dist = (); my @OTU = ();
while (my $line = <IN>){
	chomp $line;
	my @columns = split ("\t", $line); ## https://github.com/marbl/Mash/blob/master/doc/sphinx/tutorials.rst
	my $reference = $columns[0];
	my $query = $columns[1];
	my $mashdist = $columns[2];
	my $pval = $columns[3];
	my ($numerator, $denominator) = $columns[4] =~ /(\d+)\/(\d+)/;
	my $similarity = $numerator/$denominator;
	my $distance = 1 - $similarity; ## Applying an alternate method to calculate genetic distances, i.e. 1 - proportion of kmers shared; e.g. 1 - 857/1000 = 0.143.
	if (exists $hash{$reference}){
		$hash{$reference} .= "${sep}$mashdist";
		$dist{$reference} .= "${sep}$distance";
	}
	else{
		$hash{$reference} .= "${sep}$mashdist";
		$dist{$reference} .= "${sep}$distance";
		push (@OTU, $reference);
	}
}
foreach (@OTU) {my $sp = $_;  $sp =~ s/.fasta$//; print MASH "${sep}$sp"; print CORR "${sep}$sp";}
print MASH "\n"; print CORR "\n";
while (my $taxa = shift@OTU){
	my $spp = $taxa; $spp =~ s/.fasta$//;
	print MASH "$spp"."$hash{$taxa}\n";
	print CORR "$spp"."$dist{$taxa}\n";
}
