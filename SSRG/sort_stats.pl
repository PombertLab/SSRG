#!/usr/bin/perl
## Pombert lab, 2017
## v1b

use strict;
use warnings;
use Math::Complex;

my $usage = "sort_stats.pl *.stats";
die "\n$usage\n" unless @ARGV;

## Initializing table
open OUT, ">StatsTable.tsv";
print OUT "\t";

## Creating headers for species/strains/isolates
my %species = ();
foreach(@ARGV){
	if ($_ =~ /fastq.(\S+)\.\w+\.\w+\.stats$/){
		if (exists $species{$1}){next;}
		else {$species{$1} = $1; print OUT "\t$1\t\t\t";}
	}
}
print OUT "\n";
my $num = scalar (keys %species);
my $col = "\t% cov\ttotal SNPs\tSNP per KB\t";
my $head = ($col)x($num);
print OUT "$head\n";

## Iterating through/parsing stats files
my $species = undef;
my $minc = 100;
my $maxsnp = 0;
my $maxsnpkb = 0;
my $sumsnp = '0';
my $count = '0';
my $snkbsum = '0';
while (my $file = shift@ARGV){
	open IN, "<$file";
	my $sp = undef;
	my $depth = undef;
	my $snps = undef;
	my $snkb = undef;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^FASTQ file\(s\) used: (\S+)\.\d+\.fastq/){$sp = $1;}
		elsif ($line =~ /^Sequencing breadth \(percentage of bases covered by at least one read\)\s+(\S+)\%/){$depth = $1; if ($depth < $minc){$minc = $depth;}}
		elsif ($line =~ /^Total number of SNPs found: (\d+)/){$snps = $1; $sumsnp += $1; $count++; if ($snps > $maxsnp){$maxsnp = $snps;}}
		elsif ($line =~ /^Average number of SNPs per Kb: (\S+)/){$snkb = $1; $snkbsum+= $1; if ($snkb > $maxsnpkb){$maxsnpkb = $snkb;}}
	}
	if (defined $species){
		if ($species eq $sp){
			print OUT "\t\t$depth\t$snps\t$snkb";
		}
		else{
			$species = $sp;
			print OUT "\n$species";
			print OUT "\t$depth\t$snps\t$snkb";
		}
	}
	else{
		$species = $sp;
		print OUT "$species";
		print OUT "\t$depth\t$snps\t$snkb";
	}
}
my $avsnp = sprintf ("%.2f", ($sumsnp/$count));
my $avsnb = sprintf ("%.2f", ($snkbsum/$count));
my $spp = sqrt($count);
print OUT "\n\n";
print OUT "Number of species found in node:\t$spp\n";
print OUT "Minimum genome coverage\t$minc \%\n";
print OUT "Max SNP total between species\t$maxsnp\n";
print OUT "Average SNP total between species\t$avsnp\n";
print OUT "Max SNP/kb total\t$maxsnpkb\n";
print OUT "Average SNP/kb between species\t$avsnb\n";
