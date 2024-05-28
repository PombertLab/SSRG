#!/usr/bin/env perl
## Pombert JF, Illinois tech - 2019

my $version = '1.2a';
my $name = 'count_SNPs.pl';
my $updated = '2021-03-13';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $options = <<"OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Summarizes the number variants found in the VCF files

USAGE       ${name} -o OUTPUT_PREFIX -v *.vcf

OPTIONS:
-o (--output)    Desired output file prefix ## File extension (.tsv) will be appended automatically
-v (--vcf)       VCF files to summarize
--version        Show script version
OPTIONS

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my $output;
my @vcf;
my $sc_version;
GetOptions(
    'o|output=s' => \$output,
    'v|vcf=s@{1,}' => \@vcf,
    'version' => \$sc_version
);

#########################################################################
### Version
#########################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

#########################################################################
### Creating tab-delmited table
#########################################################################

open OUT, ">", "$output.tsv" or die "Can't create output file $output.tsv: $!\n";

while (my $file = shift@vcf){

    open IN, "<", "$file" or die "Can't read file $file: $!\n";

    my $snps = 0;
    while (my $line = <IN>){
        chomp $line;
        if ($line =~ /^#/){
            next;
        }
        else {
            $snps++;
        }
    }

    print OUT "$file\t$snps\n";

}
