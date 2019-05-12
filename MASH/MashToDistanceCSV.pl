#!/usr/bin/perl
## Pombert lab, 2019
my $version = '1.1c';
my $name = 'MashToDistanceCSV.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Converts the Mash output to a CSV matrix for easy plotting with R
USAGE		MashToDistanceCSV.pl -i Mash.txt -o Mash

OPTIONS:
-i (--input)		Mash output file (unsorted)
-o (--out)		Output file prefix [Default: Mash]
OPTIONS
die "$options\n" unless @ARGV;

my $mash;
my $out = 'Mash';
GetOptions(
	'i|input=s' => \$mash,
	'o|out=s' => \$out,
);

## Reformatting Mash output to a CSV matrix
open IN, "<$mash";
open MASH, ">$out.csv"; open CORR, ">$out.alternatedist.csv";
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
		$hash{$reference} .= "\,$mashdist";
		$dist{$reference} .= "\,$distance";
	}
	else{
		$hash{$reference} .= "\,$mashdist";
		$dist{$reference} .= "\,$distance";
		push (@OTU, $reference);
	}
}
foreach (@OTU) {my $sp = $_;  $sp =~ s/.fasta$//; print MASH "\,$sp"; print CORR "\,$sp";}
print MASH "\n"; print CORR "\n";
while (my $taxa = shift@OTU){
	my $spp = $taxa; $spp =~ s/.fasta$//;
	print MASH "$spp"."$hash{$taxa}\n";
	print CORR "$spp"."$dist{$taxa}\n";
}
