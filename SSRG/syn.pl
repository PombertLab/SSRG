#!/usr/bin/perl
## Pombert Lab, 2017 Illinois Tech
## Sort SNPs per type (CDS, tRNA, rRNA or intergenic) and synonymous or non-synonymous, if applicable
## v0.1

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = "USAGE = syn.pl -fa reference.fasta -gff reference.gff -vcf *.vcf";
die "\n$usage\n\n" unless @ARGV;

my %aa = ('tca'=>'S','tcc'=>'S','tcg'=>'S','tct'=>'S','tcy'=>'S','tcr'=>'S','tcw'=>'S','tcs'=>'S','tck'=>'S','tcm'=>'S','tcb'=>'S','tcd'=>'S','tch'=>'S','tcv'=>'S','tcn'=>'S',
	'ttc'=>'F','ttt'=>'F','tty'=>'F','tta'=>'L','ttg'=>'L','ttr'=>'L','tac'=>'Y','tat'=>'Y','tay'=>'Y','taa'=>'X','tag'=>'X','tga'=>'X','tgc'=>'C','tgt'=>'C','tgy'=>'C',
	'tgg'=>'W','cta'=>'L','ctc'=>'L','ctg'=>'L','ctt'=>'L','cty'=>'L','ctr'=>'L','cts'=>'L','ctw'=>'L','ctk'=>'L','ctm'=>'L','ctb'=>'L','ctd'=>'L','cth'=>'L','ctv'=>'L','ctn'=>'L',
	'cca'=>'P','ccc'=>'P','ccg'=>'P','cct'=>'P','ccr'=>'P','ccy'=>'P','ccs'=>'P','ccw'=>'P','cck'=>'P','ccm'=>'P','ccb'=>'P','ccd'=>'P','cch'=>'P','ccv'=>'P','ccn'=>'P',
	'cac'=>'H','cat'=>'H','cay'=>'H','caa'=>'Q','cag'=>'Q','car'=>'Q','ata'=>'I','atc'=>'I','att'=>'I','atw'=>'I','aty'=>'I','ath'=>'I','atg'=>'M','aac'=>'N','aat'=>'N','aay'=>'N',
	'cga'=>'R','cgc'=>'R','cgg'=>'R','cgt'=>'R','cgr'=>'R','cgy'=>'R','cgs'=>'R','cgw'=>'R','cgk'=>'R','cgm'=>'R','cgb'=>'R','cgd'=>'R','cgh'=>'R','cgv'=>'R','cgn'=>'R',
	'aca'=>'T','acc'=>'T','acg'=>'T','act'=>'T','acr'=>'T','acy'=>'T','acs'=>'T','acw'=>'T','ack'=>'T','acm'=>'T','acb'=>'T','acd'=>'T','ach'=>'T','acv'=>'T','acn'=>'T',
	'aaa'=>'K','aag'=>'K','aar'=>'K','agc'=>'S','agt'=>'S','agy'=>'S','aga'=>'R','agg'=>'R','agr'=>'R','gac'=>'D','gat'=>'D','gay'=>'D','gaa'=>'E','gag'=>'E','gar'=>'E',
	'gta'=>'V','gtc'=>'V','gtg'=>'V','gtt'=>'V','gtr'=>'V','gty'=>'V','gts'=>'V','gtw'=>'V','gtk'=>'V','gtm'=>'V','gtb'=>'V','gtd'=>'V','gth'=>'V','gtv'=>'V','gtn'=>'V',
	'gca'=>'A','gcc'=>'A','gcg'=>'A','gct'=>'A','gcr'=>'A','gcy'=>'A','gcs'=>'A','gcw'=>'A','gck'=>'A','gcm'=>'A','gcb'=>'A','gcd'=>'A','gch'=>'A','gcv'=>'A','gcn'=>'A',
	'gga'=>'G','ggc'=>'G','ggg'=>'G','ggt'=>'G','ggr'=>'G','ggy'=>'G','ggs'=>'G','ggw'=>'G','ggk'=>'G','ggm'=>'G','ggb'=>'G','ggd'=>'G','ggh'=>'G','ggv'=>'G','ggn'=>'G',
);

## GetOptions
my $fasta;
my $gff;
my @vcf;
GetOptions(
	'fa=s' => \$fasta,
	'gff=s' => \$gff,
	'vcf=s@{1,}' => \@vcf
);

my $contig;
## Building hash of sequences
open FASTA, "<$fasta";
my %sequences;
while (my $line = <FASTA>){
	chomp $line;
	if ($line =~ /^>(\S+)/){$contig = $1;}
	else{$sequences{$contig} .= $line;}
}
## Building hashes of features
open GFF, "<$gff";	
my %features;
my %hashes;		## hash of hashes for contigs and their feature locations
while (my $line = <GFF>){
	chomp $line;
	if ($line =~ /^#/){next;}
	elsif ($line =~ /^ID=/){next;}
	elsif ($line =~ /^(\S+)\t\S+\t(tRNA|rRNA|CDS|)\t(\d+)\t(\d+)\t\.\t([+-])\t[.012]\tID=(\w+)/){
		$contig = $1; my $type = $2; my $start = $3; my $end = $4; my $strand = $5; my $locus = $6;
		$features{$locus}[0] = $type;
		$features{$locus}[1] = $start;
		$features{$locus}[2] = $end;
		$features{$locus}[3] = $strand;
		for ($start..$end){
			push (@{$hashes{$contig}{$_}}, $locus);
			#print "$contig\t$_\t$locus\n"; ## Debugging
		}
	}	
}
## Working on VCF files
my $codon; my $snp; my $revcodon; my $revsnp;
while (my $vcf = shift@vcf){
	open VCF, "<$vcf";
	open SORT, ">$vcf.sorted";
	while (my $line = <VCF>){
		chomp $line;
		if ($line =~ /^#/){next;}
		elsif ($line =~ /^(\S+)\t(\d+)\t\.\t(\w+)\t(\w+)/){
			$contig = $1; my $position = $2; my $ref = $3; my $alt = $4;
			if (exists $hashes{$contig}{$position}){
				for (0..$#{$hashes{$contig}{$position}}){
					my $locus = $hashes{$contig}{$position}[$_];
					my $type = $features{$hashes{$contig}{$position}[$_]}[0];
					my $start = $features{$hashes{$contig}{$position}[$_]}[1];
					my $end = $features{$hashes{$contig}{$position}[$_]}[2];
					my $strand  = $features{$hashes{$contig}{$position}[$_]}[3];
					print SORT "$type\t";
					print "$type\t$locus\t$start\t$end\t$strand\n";
					if ($type eq 'tRNA'){print SORT "Not applicable\t\t\t";}
					elsif ($type eq 'rRNA'){print SORT "Not applicable\t\t\t";}
					elsif ($type eq 'CDS'){
						if (($strand  eq '+') || ($strand  eq 'plus')){
							#~ my $seq = substr($sequences{$contig}, $start-1, $end-$start+1); ## Debugging
							#~ print "$seq\n";
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
							#~ my $seq = substr($sequences{$contig}, $start-1, $end-$start+1); ## Debugging
							#~ my $rev = reverse($seq);
							#~ $rev =~ tr/ATGCNatgcn/TACGNtacgn/;
							#~ print "$rev\n"; 
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
					#~ print "\ncodon = $codon\talt = $snp\n"; ## Debugging
					print SORT "$locus\t$line\n";
				}
			}
			else {print SORT "intergenic\tNot applicable\t\t\t\t$line\n"}
		}
	}
}

## Subroutines
sub rev{
	$codon = reverse($revcodon); $codon=~ tr/ATGCNatgcn/TACGNtacgn/;
	$snp = reverse($revsnp); $snp=~ tr/ATGCNatgcn/TACGNtacgn/;	
}
sub translation{
	$codon = lc($codon); $snp = lc($snp);
	if ($aa{$codon} eq $aa{$snp}){print SORT "synonym\t$codon\t$snp\t";}
	elsif ($aa{$codon} ne $aa{$snp}){print SORT "non-syn\t$codon\t$snp\t";}
}