#!/usr/bin/perl
## Pombert lab, 2019
my $version = '1.2';
my $name = 'MashToDistanceMatrix.pl';
my $updated = '14/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts the Mash output to a CSV or TSV matrix for easy plotting with R

USAGE		${name} \\
		  -i MASH/S_pneumoniae_75.mash \\
		  -o MASH \\
		  -p S_pneumoniae_75 \\
		  -f csv

OPTIONS:
-i (--input)		Mash output file
-o (--outdir)	Output folder [Default: ./]
-p (--prefix)	Output file prefix [Default: Mash]
-f (--format)	Output format; csv or tsv [Default: csv]
OPTIONS
die "$usage\n" unless @ARGV;

my $mash;
my $outdir = './';
my $prefix = 'Mash';
my $format = 'csv';
GetOptions(
	'i|input=s' => \$mash,
	'o|outdir=s' => \$outdir,
	'p|prefix=s' => \$prefix,
	'f|format=s' => \$format
);

## Reformatting Mash output to a CSV or TSV matrix
$format = lc($format);
my $sep;
if ($format eq 'csv'){$sep = ",";}
elsif ($format eq 'tsv'){$sep = "\t";}
else {die "Format $format in not recognized. Please specify csv or tsv...\n";}

open IN, "<", "$mash" or die "Can't read file $mash: $!\n";
my $outfile = "${outdir}/$prefix.$format";
my $alternate = "${outdir}/$prefix.alternatedist.$format";
open MASH, ">", "$outfile" or die "Can't create file $outfile: $!\n";
open CORR, ">", "$alternate" or die "Can't create file $outfile: $!\n";
print MASH "OTU";
print CORR "OTU";

my %hash = (); my %dist = (); my @OTU = ();
while (my $line = <IN>){
	chomp $line;
	my @columns = split ("\t", $line); ## https://github.com/marbl/Mash/blob/master/doc/sphinx/tutorials.rst
	my $reference = fileparse($columns[0]);
	my $query = fileparse($columns[1]);
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
