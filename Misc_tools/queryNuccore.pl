#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2016
## Retrieve genomes/proteins from the NCBI Nucleotide database using a list of accession numbers (one per line)

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use LWP::Simple;

my $usage = "
USAGE = queryNuccore.pl [options]

EXAMPLE: queryNuccore.pl -fa -gbk -l accession.txt

OPTIONS:
-fa (--fasta)	Reference genome(s) in fasta format
-gb (--genbank)	Reference genome(s) in GenBank format
-sqn (--sequin)	Reference genome(s) in Sequin ASN format
-p (--proteins)	Protein sequences (amino acids)
-c (--cds)	Protein sequences (nucleotides)
-l (--list)	Accession numbers list, one accession per line
";
die "$usage\n" unless@ARGV;

## Options
my $fasta;
my $genbank;
my $proteins;
my $sequin;
my $cds;
my $list = '';

GetOptions(
	'fasta|fa' => \$fasta,
	'genbank|gb' => \$genbank,
	'proteins|p' => \$proteins,
	'sequin|sqn' => \$sequin,
	'cds|c' => \$cds,
	'list|l=s' => \$list,
);

my $start = localtime();
my $tstart = time;

open IN, "<$list";

my $efetch = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
my $rettype = undef;
my $retmode = undef;
my $db = 'nuccore';

while (my $accession = <IN>){
	chomp $accession;
	if ($accession =~ /^#/){next;} ## Discarding comments
	if ($fasta || $genbank || $proteins || $sequin || $cds){
		$retmode = 'text';
		if ($fasta){
			$rettype ='fasta';
			my $name = "$accession".'.fasta';
			my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
			system "echo Downloading $accession.fasta ...";
			getstore($URL, $name);
		}
		if ($genbank){
			$rettype ='gbwithparts';
			my $name = "$accession".'.gbk';
			my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
			system "echo Downloading $accession.gbk ...";
			getstore($URL, $name);
		}
		if ($sequin){
			$rettype ='';
			my $name = "$accession".'.sqn';
			my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
			system "echo Downloading $accession.sqn ...";
			getstore($URL, $name);
		}
		if ($proteins){
			$rettype ='fasta_cds_na';
			my $name = "$accession".'.fna';
			my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
			system "echo Downloading $accession.fna ...";
			getstore($URL, $name);
		}
		if ($cds){
			$rettype ='fasta_cds_aa';
			my $name = "$accession".'.faa';
			my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
			system "echo Downloading $accession.faa ...";
			getstore($URL, $name);
		}
	}
}

my $time_taken = time - $tstart;
my $end = localtime();

print "\nDownloading started on: $start\n";
print "Downloading ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";
