#!/usr/bin/perl

## Pombert lab, 2016
## Convert MASH output to CSV format

use strict;
use warnings;

my $usage = "perl MashToDistanceCSV.pl Mash.txt\n";
die $usage unless @ARGV;

while (my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.txt//;
	open OUT, ">$file.csv";
	
	my %hash = ();
	my @OTU = ();
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\d+)\/(\d+)$/){ ## https://github.com/marbl/Mash/blob/master/doc/sphinx/tutorials.rst
			my $reference = $1;
			my $query = $2;
			my $mashdist= $3;
			my $pval = $4;
			my $similarity = "$5".'/'."$6";
			my $distance = 1-($5/$6);
			if (exists $hash{$reference}){
				$hash{$reference} .= "\,$distance";
			}
			else{
				$hash{$reference} .= "\,$distance";
				push (@OTU, $reference);
			}
		}
	}
	foreach (@OTU) {print OUT "\,$_";}
	print OUT "\n";
	while (my $taxa = shift@OTU){
		print OUT "$taxa"."$hash{$taxa}\n";
	}
}