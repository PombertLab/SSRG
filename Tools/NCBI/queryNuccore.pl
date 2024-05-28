#!/usr/bin/perl
## Pombert JF,  Illinois Tech 2019

my $version = '0.5a';
my $name = 'queryNuccore.pl';
my $updated = '2024-05-28';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

#########################################################################
### Command line options
#########################################################################

my $usage = <<"OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Retrieve genomes/proteins from the NCBI Nucleotide database using a list of accession numbers
            (one per line). Based on NCBI's efetch tool - http://www.ncbi.nlm.nih.gov/books/NBK25499/

EXAMPLE     ${name} \\
              -db nuccore \\
              -o DATASETS \\
              -fa \\
              -gb \\
              -l accession.txt

OPTIONS:
-db (--database)  NCBI database to be queried [Default: nuccore]
-o (--outdir)     Output folder [Default: ./]
-fa (--fasta)     Reference genome(s) in fasta format
-gb (--genbank)   Reference genome(s) in GenBank format
-sqn (--sequin)   Reference genome(s) in Sequin ASN format
-p (--proteins)   Protein sequences (amino acids)
-c (--cds)        Protein sequences (nucleotides)
-l (--list)       Accession numbers list, one accession per line
-v (--version)    Show script version
OPTIONS

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

## Defining options
my $db = 'nuccore';
my $outdir = './';
my $fasta;
my $genbank;
my $proteins;
my $sequin;
my $cds;
my $list = '';
my $sc_version;
GetOptions(
    'db|database=s' => \$db,
    'o|outdir=s' => \$outdir,
    'fa|fasta' => \$fasta,
    'gb|genbank' => \$genbank,
    'p|proteins' => \$proteins,
    'sqn|sequin' => \$sequin,
    'c|cds' => \$cds,
    'l|list=s' => \$list,
    'v|version' => \$sc_version
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
### Output directory
#########################################################################

unless (-d $outdir){
    mkdir ($outdir,0755) or die "Can't create output folder $outdir: $!\n";
}

#########################################################################
### Downloading data from Nuccore
#########################################################################

my $start = localtime();
my $tstart = time;
my $efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi';
my $rettype;
my $retmode;
my $filename;
my $accession;

open IN, '<', $list or die "Can't read file $list: $!\n";

while ($accession = <IN>){

    chomp $accession;
    $accession =~ s/ //g;

    if ($accession =~ /^#/){ ## Discarding comments
        next;
    }

    if ($fasta || $genbank || $proteins || $sequin || $cds){

        $retmode = 'text';

        if ($fasta){
            $rettype = 'fasta';
            $filename = $accession . '.fasta';
            curl();
        }

        if ($genbank){
            $rettype = 'gbwithparts';
            $filename = $accession . '.gbk';
            curl();
        }

        if ($sequin){
            $rettype = '';
            $filename = $accession . '.sqn';
            curl();
        }

        if ($cds){
            $rettype = 'fasta_cds_na';
            $filename = $accession . '.fna';
            curl();
        }
        if ($proteins){
            $rettype = 'fasta_cds_aa';
            $filename = $accession . '.faa';
            curl();
        }
    }
}

my $time_taken = time - $tstart;
my $end = localtime();

print "\nDownloading started on: $start\n";
print "Downloading ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";

#########################################################################
### Subroutine(s)
#########################################################################

sub curl{
    print "Downloading $filename ...\n";
    my $URL = $efetch.'?db='.$db.'&id='.$accession.'&rettype='.$rettype.'&retmode='.$retmode;
    system "curl \\
      -o ${outdir}/$filename \\
      -L \'$URL\' ";
    print "\n\n";
}