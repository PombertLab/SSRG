#!/usr/bin/perl

## SSRG: Synthetic Short Read Generator
## Generates synthetic short reads in Fastq (Q33) format from multifasta files
## Pombert JF, Illinois Institute of Technology (2015)
## Iva Veseli, Illinois Institute of Technology (2016)
## Version 1.2

use strict;
use warnings;
use Getopt::Long qw(GetOptions);    


## Define command line options:
my $options = <<'END_OPTIONS';
OPTIONS: 
 -r (--readsize)	Synthetic reads size [default: 100].
                	Minimum size required for 50X coverage: 25 nt.
                	Minimum size required for 100X coverage: 50 nt.

 -c100 (--cov100)	Set sequencing depth to 100X [default: 50].

 -qs (--qscore)		Quality score associated with each base [default: 30].
 
 -q64          		Old Illumina Q64 FastQ format [default: Q33 (Sanger)].
END_OPTIONS

my $usage = "\nUSAGE = perl SSRG.pl [options] *.fasta\n\n$options";
die $usage unless(scalar@ARGV>=1);

my $start = localtime();
my $tstart = time;

## Declare options
my $winsize = 100;
my $cov100 = '';      ## Default is false so that 50X coverage is used 
my $slide = 4;        ## Determined by coverage and read size. For defaults of 100 bases and 50X, slide every 4 nucleotides.
my $qscore = 30;        
my $q64 = '';         ## Default is false so that Q33 format is used

GetOptions(
'readsize|r=i' => \$winsize,
'qscore|qs=i' => \$qscore,
'q64' => \$q64,
'cov100|c100' => \$cov100);

# Calculate the slide value - if read size is not divisible by 25, truncate decimal slide value to ensure at least 50X or 100X coverage
$slide = (2/50) * $winsize;       
if ($cov100){ $slide /= 2;}       # reduce slide to increase coverage to 100X if necessary
$slide = int($slide);             # truncate to integer value.
#print "slide = $slide\n";

my $warning1 = "Fatal Error: Read size must be an integer\n";
die $warning1 unless(int($winsize) == $winsize);
my $warning2 = "Fatal Error: Read size is too small\n";
die $warning2 unless($winsize > 0 && $slide > 0);


## If Q33 Q=0 is ASCII 33; if Q64 Q=0 is ASCII 64
my %ASCII = (
33 => '!', 34 => '"', 35 => '#', 36 => '$', 37 => '%', 38 => '&', 39 => "'", 40 => '(', 
41 => ')', 42 => '*', 43 => '+', 44 => ',', 45 => '-', 46 => '.', 47 => '/', 48 => '0',
49 => '1', 50 => '2', 51 => '3', 52 => '4', 53 => '5', 54 => '6', 55 => '7', 56 => '8',
57 => '9', 58 => ':', 59 => ';', 60 => '<', 61 => '=', 62 => '>', 63 => '?', 64 => '@',
65 => 'A', 66 => 'B', 67 => 'C', 68 => 'D', 69 => 'E', 70 => 'F', 71 => 'G', 72 => 'H', 
73 => 'I', 74 => 'J', 75 => 'K', 76 => 'L', 77 => 'M', 78 => 'N', 79 => 'O', 80 => 'P',
81 => 'Q', 82 => 'R', 83 => 'S', 84 => 'T', 85 => 'U', 86 => 'V', 87 => 'W', 88 => 'X',
89 => 'Y', 90 => 'Z', 91 => '[', 92 => "\\", 93 => ']', 94 => '^', 95 => '_', 96 => '`',
97 => 'a', 98 => 'b', 99 => 'c', 100 => 'd', 101 => 'e', 102 => 'f', 103 => 'g', 104 => 'h',
105 => 'i', 106 => 'j'
);


## take qscore, calculate appropriate offset, determine correct ASCII character for quality score
my $offset; ## ASCII offset for quality score
if ($q64)
{
    #print "Using Q64 format\n";
    $offset = $qscore + 64;
}
else
{
    #print "Using Q33 format\n";
    $offset = $qscore + 33;
}
my $score = $ASCII{$offset};  ## the artificial quality score to be assigned in the FastQ output


while (my $fasta = shift @ARGV) {
	open IN, "<$fasta" or die "cannot open $fasta";
	$fasta =~ s/\.fasta$//; $fasta =~ s/\.fsa$//; $fasta =~ s/\.fa$//; $fasta =~ s/\.fna$//; ## Removing file extensions
	open OUT, ">$fasta.$winsize.fastq";

        ##print run stats
        print "\nInput File:\t$fasta\n";
        print "Output File:\t$fasta.$winsize.fastq\n";
        print "Coverage:\t";
        if ($cov100){ print "100X\n";}
        else{ print "50X\n";}
        print "Quality score:\t$qscore\n";
        if ($q64){ print "Q64 format:\t$score\n";}
        else{ print "Q33 format:\t$score\n";}
        print "Read size:\t$winsize\n";
        print "Slide every:\t$slide nucleotides\n";
        print "\n\n";
        

	my @contigs = (); ## Initializing list of contigs per fasta file
	my %contigs = (); ## Initializing hash of contigs per fasta file: key = contig name; value = sequence
	my $name = undef; ## Contig name to be memorized
	my $count = 0; ## Read number counter to be auto-incremented
	
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(.*)$/){            ## if line is a header, save the new contig name and add it to hash
			$name = $1;
			push(@contigs, $name);
			$contigs{$name} = undef;
		}
		else{
			$contigs{$name} .= $line;   ## otherwise, append the line to the hash value of the current contig
		}
	}
	
        ## Parse each contig in the forward and reverse directions 
	while (my $todo = shift@contigs){
		chomp $todo;
		my $len = length($contigs{$todo});
		my $seq = $contigs{$todo};
                my $rev = reverse($seq);            ## get the reverse complement of the contig
                $rev =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
                my $i = 0;
		while($i <= ($len-$winsize)) { 
			my $read = substr($seq, $i, $winsize);
                        my $rev_read = substr($rev, $i, $winsize);
			$count++;
                        ## FastQ output:
			print OUT '@SYNTHREAD_'."$count\n";
			print OUT "$read\n";
			print OUT '+'."\n";
			print OUT "$score" x $winsize;   ## assign fake quality score
			print OUT "\n";
                        $count++;
                        print OUT '@SYNTHREAD_'."$count\n";
                        print OUT "$rev_read\n";
			print OUT '+'."\n";
			print OUT "$score" x $winsize;
			print OUT "\n";
                        $i += $slide;                       
		}
	}
}

my $end = localtime();
my $time_taken = time -$tstart;

print "\nSSRG started on: $start\n";
print "SSRG ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";

close IN;
close OUT;

exit;
