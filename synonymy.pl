#!/usr/bin/perl
## Pombert Lab, 2017-2018 Illinois Tech
my $version = '0.4a';
my $name = 'synonymy.pl';
my $updated = '14/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

## Usage definition
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Sort SNPs per type (CDS, tRNA, rRNA or intergenic) and synonymous or non-synonymous, if applicable
NOTE		Does not work with introns yet
		
USAGE		${name} \\
		  -gcode 1 \\
		  -fa reference.fasta \\
		  -ref reference.gbk \\
		  -format gb \\
		  -vcf *.vcf \\
		  -o SYNONYMY \\
		  -p synonymy

OPTIONS:
-fa (--fasta)	Reference genome in fasta format
-r (--ref)	Reference genome annotation in GBK or GFF format
-f (--format)	Reference file format; gb or gff [default: gff]
-vcf (--vcf)	SNPs in Variant Calling Format (VCF)
-o (--outdir)	Output directory [Default: ./]
-p (--prefix)	Table prefix [Default: synonymy]
-gc (--gcode)	NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		# For complete list, see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-v (--verbose)	Adds verbosity
OPTIONS
die "\n$usage\n" unless@ARGV;

## GetOptions
my $fasta;
my $ref;
my $format = 'gff';
my @vcf;
my $outdir = './';
my $prefix = 'synonymy';
my $gc = '1';
my $verbose;
GetOptions(
	'fa|fasta=s' => \$fasta,
	'r|ref=s' => \$ref,
	'f|format=s' => \$format,
	'vcf=s@{1,}' => \@vcf,
	'o|outdir=s' => \$outdir,
	'p|prefix=s' => \$prefix,
	'gc|gcode=i' => \$gc,
	'v|verbose' => \$verbose
);

## Building hash of sequences
open FASTA, "<", "$fasta" or die "Can't read FASTA file $fasta: $!\n";
my %sequences; my $contig;
while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){ $contig = $1; }
	else{ $sequences{$contig} .= $line; }
}

## Creating output directory + output files
unless (-d $outdir){ mkdir ($outdir,0755) or die "Can't create folder $outdir: $!\n"; }
open FEAT, ">", "${outdir}/$prefix.features.tsv" or die "Can't create file ${outdir}/$prefix.features.tsv: $!\n";
open INTER, ">", "${outdir}/$prefix.intergenic.tsv" or die "Can't create file ${outdir}/$prefix.intergenic.tsv: $!\n";

my @feat_headers = (
	'Locus_tag',
	'Reference Contig',
	'VCF file',
	'Type',
	'Start',
	'End',
	'Strandedness',
	'Product',
	'Variant Position',
	'Reference',
	'Variant',
	'Synonymy',
	'Genetic Code',
	'RefCodon',
	'VarCodon',
	'Frequency',
	'E-value'
);
for (0..$#feat_headers-1){ print FEAT "$feat_headers[$_]\t"; }
print FEAT "$feat_headers[-1]\n";

my @inter_headers = (
	'Reference Contig',
	'VCF file',
	'Type',
	'Variant Position',
	'Reference',
	'Variant',
	'Frequency',
	'E-value'
);
for (0..$#inter_headers-1){ print INTER "$inter_headers[$_]\t"; }
print INTER "$inter_headers[-1]\n";

## Building hashes of features. Must simplify the data structure to one hash of features. Will be easier to adapt for input files. Working multiline would do it.
open REF, "<", "$ref" or die "Can't read reference file $ref: $!\n";
my %genes, my %features; my %hashes; ## hash of hashes for contigs and their feature locations
my $type; my $start; my $end; my $strand; my $id; my $locus_tag; my $parent; my $product; my $key;
if ($format eq 'gff'){
	while (my $line = <REF>){
		chomp $line;
		## Looking at GFF files
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^ID=/){next;}
		elsif ($line =~ /^(\S+)\t.*\t(gene|pseudogene)\t(\d+)\t(\d+)\t\.\t([+-])\t[.012]\tID=(\w+).*\;locus_tag=([A-Za-z0-9 _]+)/){ ## Working on locus_tags
			$contig = $1; $type = $2; $start = $3; $end = $4; $strand = $5; $id = $6; $locus_tag = $7;
			$genes{$id}[0] = $contig;
			$genes{$id}[1] = $start;
			$genes{$id}[2] = $end;
			$genes{$id}[3] = $strand;
			$genes{$id}[4] = $locus_tag;
		}
		elsif ($line =~ /^(\S+)\t.*\t(tRNA|rRNA|CDS|)\t(\d+)\t(\d+)\t\.\t([+-])\t[.012]\tID=(\w+).*Parent=([A-Za-z0-9_]+).*product=([A-Za-z0-9 -]+)/){	## Working on features
			$contig = $1; $type = $2; $start = $3; $end = $4; $strand = $5; $id = $6; $parent = $7; $product = $8;
			$features{$id}[0] = $type;
			$features{$id}[1] = $start;
			$features{$id}[2] = $end;
			$features{$id}[3] = $strand;
			$features{$id}[4] = $parent;
			$features{$id}[5] = $product;
			if ($verbose){ print "$contig\t$type\t$strand\t$locus_tag\t$product\n"; }
			for ($start..$end){
				push (@{$hashes{$contig}{$_}}, $id);
			}
		}
	}
}
elsif ($format eq 'gb'){
	my $gbk; while (my $line = <REF>){$gbk .= $line;}
	my @contigs = split ("\/\/\n", $gbk);
	while (my $cg = shift @contigs){
		my @data = split ("ORIGIN.*?\n", $cg); ## $data[0] => annotations, $data[1] => sequences
		my $sequence = $data[1]; $sequence =~ s/[0-9\s\n]//g;
		my $locus; ($locus) = ($data[0] =~ /LOCUS\s+(\S+)/);
		my $contig; ($contig) = ($data[0] =~ /VERSION\s+(\S+)/);
		if ($verbose){ print "Working on $locus version $contig\n"; }
		my $feat; ($feat) = ($data[0] =~ /FEATURES(.*)/s);
		my @features = split("     gene            ", $feat);
		while (my $line = shift@features){
			if ($line =~ /\s+(CDS|rRNA|tRNA)\s+.*?(\d+)\.\.(\d+)/m){
				$type = $1; $start = $2; $end = $3;
				if ($line =~ /complement/){$strand = '-';}
				else{$strand = '+';}
				($locus_tag) = ($line =~ /\/locus_tag="(.*?)"/s); $locus_tag =~ s/\s{2,}/ /g; $locus_tag =~ s/\n//g; ## removing new lines, if any
				($product) = ($line =~ /\/product="(.*?)"/s); $product =~ s/\s{2,}/ /g; $product =~ s/\n//g; ## removing new lines, if any
				if ($verbose){ print "$contig\t$type\t$strand\t$locus_tag\t$product\n"; }
				$genes{$locus_tag}[0] = $contig;
				$genes{$locus_tag}[1] = $start; $features{$locus_tag}[1] = $start;
				$genes{$locus_tag}[2] = $end; $features{$locus_tag}[2] = $end;
				$genes{$locus_tag}[3] = $strand; $features{$locus_tag}[3] = $strand;
				$genes{$locus_tag}[4] = $locus_tag; $features{$locus_tag}[4] = $locus_tag;
				$features{$locus_tag}[0] = $type;
				$features{$locus_tag}[5] = $product;
				for ($start..$end){
					push (@{$hashes{$contig}{$_}}, $locus_tag);
				}
			}
		}
	}
}
## Working on VCF files
my $codon; my $snp; my $revcodon; my $revsnp; my $locus;
my %gcodes; gcodes();
while (my $vcf = shift@vcf){
	open VCF, "<", "$vcf" or die "Can't read VCF file $vcf: $!\n";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		else{
			my @vcf = split("\t", $line);
			$contig = $vcf[0];
			my $position = $vcf[1];
			my $ref = $vcf[3];
			my $alt = $vcf[4];
			
			my @params = split(":", $vcf[9]);
			my $freq = $params[6];
			my $evalue = $params[7];
			
			if (exists $hashes{$contig}{$position}){
				for (0..$#{$hashes{$contig}{$position}}){
					$locus = $hashes{$contig}{$position}[$_];
					$type = $features{$hashes{$contig}{$position}[$_]}[0];
					$start = $features{$hashes{$contig}{$position}[$_]}[1];
					$end = $features{$hashes{$contig}{$position}[$_]}[2];
					$strand = $features{$hashes{$contig}{$position}[$_]}[3];
					$parent = $features{$hashes{$contig}{$position}[$_]}[4];
					$product = $features{$hashes{$contig}{$position}[$_]}[5];
					print FEAT "$genes{$parent}[4]\t$contig\t$vcf\t$type\t$start\t$end\t$strand\t$product\t$position\t$ref\t$alt\t";
					if ($type eq 'tRNA'){print FEAT "N\/A\tN\/A\tN\/A\tN\/A\t";}
					elsif ($type eq 'rRNA'){print FEAT "N\/A\tN\/A\tN\/A\tN\/A\t";}
					elsif ($type eq 'CDS'){
						if (($strand  eq '+') || ($strand  eq 'plus')){
							if (($position - $start) == 0){
								$codon = substr($sequences{$contig}, $position-1, 3);
								$snp = $codon; $snp =~ substr($snp, 0, 1, $alt);
								translation();
							}
							elsif (($position - $start) == 1){
								$codon = substr($sequences{$contig}, $position-2, 3);
								$snp = $codon; $snp =~ substr($snp, 1, 1, $alt);
								translation();
							}
							elsif (($position - $start) > 1){
								my $location = $position - $start + 1;
								if (($location % 3) == 0){ ## Third codon position
									$codon = substr($sequences{$contig}, $position-3, 3); ## this is OK
									$snp = $codon; $snp =~ substr($snp, 2, 1, $alt);
									translation();
								}
								elsif (($location % 3) == 1){ ## First codon position
									$codon = substr($sequences{$contig}, $position-1, 3);
									$snp = $codon; $snp =~ substr($snp, 0, 1, $alt);
									translation();
								}
								elsif (($location % 3) == 2){ ## Second codon position
									$codon = substr($sequences{$contig}, $position-2, 3);
									$snp = $codon; $snp =~ substr($snp, 1, 1, $alt);
									translation();
								}
							}
						}
						elsif (($strand  eq '-') || ($strand  eq 'minus')){
							if (($end - $position) == 0){
								$revcodon = substr($sequences{$contig}, $position-3, 3);
								$revsnp = $revcodon; $revsnp =~ substr($revsnp, 2, 1, $alt);
								rev(); translation();
							}
							elsif (($end - $position) == 1){
								$revcodon  = substr($sequences{$contig}, $position-2, 3);
								$revsnp = $revcodon; $revsnp =~ substr($revsnp, 1, 1, $alt);
								rev(); translation();
							}
							elsif (($end - $position) > 1){
								my $location = $end - $position + 1;
								if (($location % 3) == 0){ ## Third codon position
									$revcodon = substr($sequences{$contig}, $position-1, 3);
									$revsnp = $revcodon; $revsnp =~ substr($revsnp, 0, 1, $alt);
									rev(); translation();
								}
								elsif (($location % 3) == 1){ ## First codon position
									$revcodon = substr($sequences{$contig}, $position-3, 3);
									$revsnp = $revcodon; $revsnp =~ substr($revsnp, 2, 1, $alt);
									rev(); translation();
								}
								elsif (($location % 3) == 2){ ## Second codon position
									$revcodon = substr($sequences{$contig}, $position-2, 3);
									$revsnp = $revcodon; $revsnp =~ substr($revsnp, 1, 1, $alt);
									rev(); translation();
								}
							}
						}
					}
					print FEAT "$freq\t$evalue\n";
				}
			}
			else {print INTER "$contig\t$vcf\tintergenic\t$position\t$ref\t$alt\t$freq\t$evalue\n"}
		}
	}
}

## Subroutines
sub rev{
	$codon = reverse($revcodon); $codon=~ tr/ATGCNatgcn/TACGNtacgn/;
	$snp = reverse($revsnp); $snp=~ tr/ATGCNatgcn/TACGNtacgn/;	
}
sub translation{
	$codon = uc($codon); $snp = uc($snp);
	if ($codon =~ /[^ATGCatgc]/){print FEAT "ambiguous bases\t$gc\t$codon\t$snp\t";}
	elsif ($gcodes{$gc}{$codon} eq $gcodes{$gc}{$snp}){print FEAT "synonym\t$gc\t$codon\t$snp\t";}
	elsif ($gcodes{$gc}{$codon} ne $gcodes{$gc}{$snp}){print FEAT "non-syn\t$gc\t$codon\t$snp\t";}
}
sub gcodes{ ## NCBI Genetic codes
	%gcodes = (
		1 => { ## The Standard Code (transl_table=1)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',   
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		2 => { ## The Vertebrate Mitochondrial Code (transl_table=2)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => '*',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => '*',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		3 => { ## The Yeast Mitochondrial Code (transl_table=3)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'T', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'T', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'T', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'T', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		4 => { ## The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
		},
		5 => { ## The Invertebrate Mitochondrial Code (transl_table=5)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		6 => { ## The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C', 
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C', 
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => '*', 
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W', 
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R', 
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R', 
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R', 
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R', 
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S', 
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S', 
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R', 
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R', 
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G', 
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G', 
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G', 
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
		},
		9 => { ## The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		10 => { ## The Euplotid Nuclear Code (transl_table=10)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'C',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		11 => { ## The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		12 => { ## The Alternative Yeast Nuclear Code (transl_table=12)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'S', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		13 => { ## The Ascidian Mitochondrial Code (transl_table=13)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'G',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'G',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
		},
		14 => { ## The Alternative Flatworm Mitochondrial Code (transl_table=14)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
		},
		16 => { ## Chlorophycean Mitochondrial Code (transl_table=16)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		21 => { ## Trematode Mitochondrial Code (transl_table=21)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		22 => { ## Scenedesmus obliquus Mitochondrial Code (transl_table=22)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => '*', 'TAA' => '*', 'TGA' => '*',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		23 => { ## Thraustochytrium Mitochondrial Code (transl_table=23)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => '*', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		24 => { ## Pterobranchia Mitochondrial Code (transl_table=24)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'K',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		25 => { ## Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'G',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		26 => { ## Pachysolen tannophilus Nuclear Code (transl_table=26)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'A', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		27 => { ## Karyorelict Nuclear (transl_table=27)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		28 => { ## Condylostoma Nuclear (transl_table=28)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		29 => { ## Mesodinium Nuclear (transl_table=29)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Y', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		30 => { ## Peritrich Nuclear (transl_table=30)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		31 => { ## Blastocrithidia Nuclear (transl_table=31)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
	);
}