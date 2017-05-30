#!/usr/bin/perl

## Pombert lab, 2016
## Convert MASH output to CSV format
## Version 1.1b

use strict;
use warnings;

my $usage = "perl MashToDistanceCSV.pl Mash.txt\n";
die $usage unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.txt//;
	open MASH, ">$file.mashdist.csv";
	open CORR, ">$file.alternatedist.csv";
	
	my %hash = ();
	my %dist = ();
	my @OTU = ();
	
	print MASH "OTU";
	print CORR "OTU";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\d+)\/(\d+)$/){ ## https://github.com/marbl/Mash/blob/master/doc/sphinx/tutorials.rst
			my $reference = $1;
			my $query = $2;
			my $mashdist= $3;
			my $pval = $4;
			my $similarity = "$5".'/'."$6";
			my $distance = 1-($5/$6); ## Applying an alternate method to calculate genetic distances, i.e. 1 - proportion of kmers shared; e.g. 1 - 857/1000 = 0.143.
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
	}
	foreach (@OTU) {my $sp = $_;  $sp =~ s/.fasta$//; print MASH "\,$sp"; print CORR "\,$sp";}
	print MASH "\n";
	print CORR "\n";
	while (my $taxa = shift@OTU){
		my $spp = $taxa; $spp =~ s/.fasta$//;
		print MASH "$spp"."$hash{$taxa}\n";
		print CORR "$spp"."$dist{$taxa}\n";
	}
}