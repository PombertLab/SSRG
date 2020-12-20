#!/usr/bin/perl
## Pombert JF, Illinois tech - 2019
my $version = '1.2';
my $name = 'count_SNPs.pl';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Summarizes the number variants found in the VCF files
USAGE		count_SNPs.pl -o OUTPUT_PREFIX -v *.vcf

OPTIONS:
-o (--output)	Desired output file prefix ## File extension (.tsv) will be appended automatically
-v (--vcf)	VCF files to summarize
OPTIONS
die "$options\n" unless@ARGV;

my $output;
my @vcf;
GetOptions(
	'o|output=s' => \$output,
	'v|vcf=s@{1,}' => \@vcf
);

## Creating tab-delmited table
open OUT, ">$output.tsv";
while (my $file = shift@vcf){
	open IN, "<$file";
	my $snps = 0;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else {$snps++;}
	}
	print OUT "$file\t$snps\n";
}
