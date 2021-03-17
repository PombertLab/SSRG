#!/usr/bin/perl
## Pombert Lab, 2017-2018 Illinois Tech
my $version = '0.5a';
my $name = 'synonymy.pl';
my $updated = '16/03/2021';

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

## Creating output directory
unless (-d $outdir){ mkdir ($outdir,0755) or die "Can't create folder $outdir: $!\n"; }

## Building hashes of features
open REF, "<", "$ref" or die "Can't read reference file $ref: $!\n";
my %features; 
my %loci; ## loci -> contig -> genomic_locus (e.g. 79553) = corresponding feature
my $type; my $start; my $end; my $strand; my $gene_id; my $locus_tag; my $parent; my $product;

if ($format eq 'gff'){
	while (my $line = <REF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^ID=/){next;}

		my @columns = split ("\t", $line);

		$contig = $columns[0];
		$type = $columns[2];
		$start = $columns[3];
		$end = $columns[4];
		$strand = $columns[6];

		## Putting descriptions into an hash for easy parsing and access
		my $description = $columns[8];
		my %descriptors;
		my @descriptors = split (";", $description);
		foreach my $desc (@descriptors){
			my ($key, $value) = $desc =~ /(.*?)=(.*)/;
			$descriptors{$key} = $value;
		}

		$gene_id = $descriptors{'ID'};

		if ($type =~ /tRNA|rRNA|CDS/){
			$locus_tag = $descriptors{'locus_tag'};
			$product = $descriptors{'product'};

			## Populating features DB
			$features{$gene_id}{'type'} = $type;
			$features{$gene_id}{'start'} = $start;
			$features{$gene_id}{'end'} = $end;
			$features{$gene_id}{'strand'}= $strand;
			$features{$gene_id}{'locus_tag'} = $locus_tag;
			$features{$gene_id}{'product'} = $product;

			if ($verbose){ print "$contig\t$type\t$strand\t$locus_tag\t$product\n"; }

			## creating a hash of genomic loci and their associated gene_ID
			for ($start..$end){
				push (@{$loci{$contig}{$_}}, $gene_id);
			}
		}
	}
}

elsif ($format eq 'gb'){
	
	## Putting the gbk file in a single string, then split per conting
	my $gbk; while (my $line = <REF>){ $gbk .= $line; }
	my @contigs = split ("\/\/\n", $gbk);

	## Parsing contigs
	while (my $cg = shift @contigs){

		my @data = split ("ORIGIN.*?\n", $cg);

		## $data[0] => annotations
		my ($locus) = $data[0] =~ /LOCUS\s+(\S+)/;
		my ($contig) = $data[0] =~ /VERSION\s+(\S+)/;
		my ($feat) = $data[0] =~ /FEATURES(.*)/s;

		## $data[1] => sequences
		my $sequence = $data[1];
		$sequence =~ s/[0-9\s\n]//g;

		if ($verbose){ print "Working on $locus version $contig\n"; }

		## Working on features
		my @features = split("     gene            ", $feat);
		while (my $line = shift@features){
			if ($line =~ /\s+(CDS|rRNA|tRNA)\s+.*?(\d+)\.\.(\d+)/m){
				$type = $1;
				$start = $2;
				$end = $3;

				if ($line =~ /complement/){ $strand = '-'; }
				else{ $strand = '+'; }

				($locus_tag) = $line =~ /\/locus_tag="(.*?)"/s;
				$locus_tag =~ s/\s{2,}/ /g;
				$locus_tag =~ s/\n//g; ## removing new lines, if any

				($product) = $line =~ /\/product="(.*?)"/s;
				$product =~ s/\s{2,}/ /g;
				$product =~ s/\n//g; ## removing new lines, if any

				$features{$locus_tag}{'type'} = $type;
				$features{$locus_tag}{'start'} = $start;
				$features{$locus_tag}{'end'} = $end;
				$features{$locus_tag}{'strand'} = $strand;
				$features{$locus_tag}{'locus_tag'} = $locus_tag;
				$features{$locus_tag}{'product'} = $product;

				if ($verbose){ print "$contig\t$type\t$strand\t$locus_tag\t$product\n"; }

				## creating a hash of genomic loci and their associated gene_ID
				for ($start..$end){
					push (@{$loci{$contig}{$_}}, $locus_tag);
				}
			}
		}
	}
}

## Populating genetic code from subroutine
my %gcodes; gcodes();

## Working on VCF files
my $codon; my $snp; my $revcodon; my $revsnp; my $locus;
while (my $vcf = shift@vcf){
	
	## Reading VCF file
	open VCF, "<", "$vcf" or die "Can't read VCF file $vcf: $!\n";
	my $basename = fileparse($vcf);
	$basename =~ s/\.\w+$//;

	## Creating output tables per VCF
	my $feature_table = "${outdir}/$basename.$prefix.features.tsv";
	my $intergenic_table = "${outdir}/$basename.$prefix.intergenic.tsv";
	open FEAT, ">", "$feature_table" or die "Can't create $feature_table: $!\n";
	open INTER, ">", "$intergenic_table" or die "Can't create $intergenic_table: $!\n";

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
			
			if (exists $loci{$contig}{$position}){
				for (0..$#{$loci{$contig}{$position}}){

					$locus = $loci{$contig}{$position}[$_];

					$type = $features{$locus}{'type'};
					$start = $features{$locus}{'start'};
					$end = $features{$locus}{'end'};
					$strand = $features{$locus}{'strand'};
					$locus_tag = $features{$locus}{'locus_tag'};
					$product = $features{$locus}{'product'};

					print FEAT "$locus_tag\t$contig\t$vcf\t$type\t$start\t$end\t$strand\t$product\t$position\t$ref\t$alt\t";

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
	$codon = reverse($revcodon); $codon =~ tr/ATGCNatgcn/TACGNtacgn/;
	$snp = reverse($revsnp); $snp =~ tr/ATGCNatgcn/TACGNtacgn/;
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