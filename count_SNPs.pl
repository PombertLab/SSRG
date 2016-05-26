#!/usr/bin/perl
# v1.1

use strict;
use warnings;

my $usage = "perl count_SNPs.pl OUTPUT_PREFIX *.vcf";
die "$usage\n" unless@ARGV;

my $prefix=shift@ARGV;
open OUT, ">$prefix.SNPs.txt";

while (my $file = shift@ARGV){
	open IN, "<$file";
	my $snps = 0;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else {$snps++;}
	}
	print OUT "$file\t$snps\n";
}
