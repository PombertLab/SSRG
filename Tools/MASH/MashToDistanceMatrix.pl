#!/usr/bin/perl
## Pombert lab, 2019
my $version = '1.3';
my $name = 'MashToDistanceMatrix.pl';
my $updated = '2022-01-22';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);

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
-i (--input)	Mash output file
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

## Outdir
unless (-d $outdir) { 
	make_path( $outdir, { mode => 0755 }	) or die "Can't create $outdir: $!\n";
}

## Choosing CSV or TSV
$format = lc($format);
my $sep;
if ($format eq 'csv'){
	$sep = ",";
}
elsif ($format eq 'tsv'){
	$sep = "\t";
}
else {
	die "Format $format in not recognized. Please specify csv or tsv...\n";
}

## Working on Mash input file
## https://github.com/marbl/Mash/blob/master/doc/sphinx/tutorials.rst

my %mash = ();
my %dist = ();
my @OTU = ();

open IN, "<", "$mash" or die "Can't read file $mash: $!\n";
while (my $line = <IN>){

	chomp $line;
	my @columns = split ("\t", $line); 

	my $reference = fileparse($columns[0]);
	$reference =~ s/.fasta$//;

	my $query = fileparse($columns[1]);
	$query =~ s/.fasta$//;
	
	my $mashdist = $columns[2];
	my $pval = $columns[3];
	
	my ($numerator, $denominator) = $columns[4] =~ /(\d+)\/(\d+)/;
	
	my $similarity = $numerator/$denominator;

	## Applying an alternate method to calculate genetic distances,
	## i.e. 1 - proportion of kmers shared; e.g. 1 - 857/1000 = 0.143.
	my $distance = 1 - $similarity;
	$distance = sprintf("%.3f", $distance);
	
	unless (exists $mash{$reference}){
		push (@OTU, $reference);
	}

	$mash{$reference}{$query} = $mashdist;
	$dist{$reference}{$query} = $distance;

}

## Creating matrix header
my $outfile = "${outdir}/$prefix.$format";
my $alternate = "${outdir}/$prefix.alternatedist.$format";
open MASH, ">", "$outfile" or die "Can't create file $outfile: $!\n";
open CORR, ">", "$alternate" or die "Can't create file $outfile: $!\n";

print MASH "OTU";
print CORR "OTU";
for (0..$#OTU) {
	my $sp = $OTU[$_];
	print MASH "$sep"."$sp";
	print CORR "$sep"."$sp";
}
print MASH "\n";
print CORR "\n";

## Creating the matrix core
for (0..$#OTU){

	my $spp = $OTU[$_];
	print MASH "$spp";
	print CORR "$spp";

	for (0..$#OTU){
		my $otu = $OTU[$_];
		my $mvalue = $mash{$spp}{$otu};
		my $dvalue = $dist{$spp}{$otu};
		print MASH "$sep"."$mvalue";
		print CORR "$sep"."$dvalue";
	}
	print MASH "\n";
	print CORR "\n";
}
