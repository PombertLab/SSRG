#!/usr/bin/perl
## SSRG: Synthetic Short Read Generator; generates synthetic short reads in Fastq (Q33) format from multifasta files
## Pombert Lab, Illinois Tech (2015-2018)
my $version = 'Version 1.5';

use strict; use warnings; use Getopt::Long qw(GetOptions);    

my $options = <<'END_OPTIONS';
FUNCTION = SSRG.pl generates synthetic short reads in Fastq (Q33) format from multifasta files

EXAMPLE (simple): SSRG.pl -f *.fasta
EXAMPLE (advanced): SSRG.pl -f *.fasta -r 250 -i 350 -s 20 -m rand -c 25 

OPTIONS: 
 -v (--version)		Display SSRG.pl version number
 -f (--fasta)		Fasta/multifasta input file
 -r (--readsize)	Synthetic reads size [default: 100]
 -m (--mode)		Sliding windows or random selection (slide or rand) [default: rand]
 -t (--type)		Single or paired ends (se or pe) [default: pe]
 -i (--insert)		Insert size for pe [default: 250]
 -s (--sdev)		Insert size standard deviation percentage [default: 10]
 -c (--coverage)	Set sequencing depth [default: 50]
 -qs (--qscore)		Quality score associated with each base [default: 30]
 -q64          		Used the old Illumina Q64 FastQ format instead of the default Q33 Sanger/Illumina 1.8+ encoding

END_OPTIONS
my $usage = "\nUSAGE = perl SSRG.pl [options] -f *.fasta\n\n$options"; die $usage unless(scalar@ARGV>=1);

my $start = localtime(); my $tstart = time;

## Declare options
my $ver;
my @fasta;
my $readsize = 100;
my $type = 'pe';
my $insert = 250;
my $sdev = 10;
my $mode = 'rand';
my $slide;
my $cov = '50';
my $qscore = 30;        
my $q64 = '';		## Default is false so that Q33 format is used
GetOptions(
	'v|version' => \$ver,
	'f|fasta=s@{1,}' => \@fasta,
	'readsize|r=s' => \$readsize,
	'm|mode=s' => \$mode,
	'i|insert=i' => \$insert,
	's|sdev=i' => \$sdev,
	't|type=s' => \$type,
	'qscore|qs=i' => \$qscore,
	'q64' => \$q64,
	'c|coverage=s' => \$cov
);
if ($ver){die "SSRG.pl $version\n";}
my $warning1 = "\nFatal Error: Read size must be an integer\n"; die $warning1 unless(int($readsize) == $readsize);
my $warning2 = "\nFatal Error: Insert size \($insert\) is smaller than read size \($readsize\). Please use a larger insert size.\n\n"; if ($type eq 'pe'){die $warning2 if ($insert < $readsize);}
if ($mode eq 'slide'){slidcov();}
my %ASCII = ASCII(); my $offset; my $score; qscore(); ## ASCII offset for quality score

## Working on files
my $fasta;		## Init fasta file
my $count = 0;	## Read number counter to be auto-incremented
my $read; my $rev_read; my $pe1; my $pe2; my $rev1; my $rev2; ## Init read variables
while ($fasta = shift @fasta) {
	open IN, "<$fasta" or die "cannot open $fasta";
	$fasta =~ s/\.fasta$//; $fasta =~ s/\.fsa$//; $fasta =~ s/\.fa$//; $fasta =~ s/\.fna$//; ## Removing file extensions
	if($type eq 'se'){open SE, ">$fasta.$readsize.SE.fastq";}
	elsif($type eq 'pe'){open R1, ">$fasta.$readsize.R1.fastq"; open R2, ">$fasta.$readsize.R2.fastq";}
	stats();

	my @contigs = ();	## Initializing list of contigs per fasta file
	my %contigs = ();	## Initializing hash of contigs per fasta file: key = contig name; value = sequence
	my $name = undef;	## Contig name to be memorized
	
	while (my $line = <IN>){	## Creating list of contigs and hash of sequences
		chomp $line;
		if ($line =~ /^>(.*)$/){
			$name = $1;
			push(@contigs, $name);
			$contigs{$name} = undef;
		}
		else{$contigs{$name} .= $line;}
	}
	
	while (my $todo = shift@contigs){ ## Parse each contig in the forward and reverse directions 
		chomp $todo;
		my $len = length($contigs{$todo});
		my $seq = $contigs{$todo};
		my $rev = reverse($seq);	## get the reverse complement of the contig
		$rev =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
		if (($mode eq 'slide') && (($type eq 'se'))){
			my $i = 0;
			while($i <= ($len-$readsize)) { 
				$read = substr($seq, $i, $readsize);
				$rev_read = substr($rev, $i, $readsize);
				$count++;
				se();
				$i += $slide;                       
			}
		}
		elsif (($mode eq 'slide') && (($type eq 'pe'))){
			my $i = 0;
			my $gap = $insert - ($readsize*2);
			my $dev = int($insert/$sdev);
			while($i <= ($len-$insert-($dev/2))){
				my $rand = int(rand($dev) - ($dev/2));
				$pe1 = substr($seq, $i, $readsize); $pe2 = substr($seq, $i+$gap+$rand+$readsize, $readsize); $pe2 = reverse($pe2); $pe2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
				$rev1 = substr($rev, $i, $readsize); $rev2 = substr($rev, $i+$gap+$rand+$readsize, $readsize); $rev2 = reverse($rev2); $rev2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
				pe();
				$i += $slide;                       
			}
		}
		elsif (($mode eq 'rand') && ($type eq 'se')){
			my $numreads = int(($len*$cov)/$readsize);	## Number of reads to achieve desired depth
			my $stdreads = int($numreads/2);			## Number of reads per strand +/-
			for (1..$stdreads){
				my $max = $len - $readsize -1;
				my $plus = int(rand($max));
				my $minus = int(rand($max));
				$read = substr($seq, $plus, $readsize);
				$rev_read = substr($rev, $minus, $readsize);
				se();
			}
		}
		elsif (($mode eq 'rand') && ($type eq 'pe')){
			my $numreads = int(($len*$cov)/$readsize);	## Number of reads to achieve desired depth
			my $stdreads = int($numreads/4);			## Number of reads per strand +/-
			my $gap = $insert - ($readsize*2);
			my $dev = int($insert/$sdev);
			for (1..$stdreads){
				my $max = $len - $insert -$dev -1;
				my $rand = int(rand($dev) - ($dev/2));
				my $plus = int(rand($max));
				my $minus = int(rand($max));
				$pe1 = substr($seq, $plus, $readsize); $pe2 = substr($seq, $plus+$gap+$rand+$readsize, $readsize); $pe2 = reverse($pe2); $pe2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
				$rev1 = substr($rev, $minus, $readsize); $rev2 = substr($rev, $minus+$gap+$rand+$readsize, $readsize); $rev2 = reverse($rev2); $rev2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
				pe();
			}
		}
	}
}

my $end = localtime(); my $time_taken = time - $tstart;
print "\nSSRG started on: $start\nSSRG ended on: $end\nTime elapsed: $time_taken seconds\n";
close IN; 
if($type eq 'se'){close SE;}
elsif($type eq 'pe'){close R1; close R2;}
exit;

## Subroutines
sub ASCII{	## If Q33 Q=0 is ASCII 33; if Q64 Q=0 is ASCII 64
	%ASCII = (33 => '!', 34 => '"', 35 => '#', 36 => '$', 37 => '%', 38 => '&', 39 => "'", 40 => '(', 
	41 => ')', 42 => '*', 43 => '+', 44 => ',', 45 => '-', 46 => '.', 47 => '/', 48 => '0',
	49 => '1', 50 => '2', 51 => '3', 52 => '4', 53 => '5', 54 => '6', 55 => '7', 56 => '8',
	57 => '9', 58 => ':', 59 => ';', 60 => '<', 61 => '=', 62 => '>', 63 => '?', 64 => '@',
	65 => 'A', 66 => 'B', 67 => 'C', 68 => 'D', 69 => 'E', 70 => 'F', 71 => 'G', 72 => 'H', 
	73 => 'I', 74 => 'J', 75 => 'K', 76 => 'L', 77 => 'M', 78 => 'N', 79 => 'O', 80 => 'P',
	81 => 'Q', 82 => 'R', 83 => 'S', 84 => 'T', 85 => 'U', 86 => 'V', 87 => 'W', 88 => 'X',
	89 => 'Y', 90 => 'Z', 91 => '[', 92 => "\\", 93 => ']', 94 => '^', 95 => '_', 96 => '`',
	97 => 'a', 98 => 'b', 99 => 'c', 100 => 'd', 101 => 'e', 102 => 'f', 103 => 'g', 104 => 'h',
	105 => 'i', 106 => 'j');
}
sub slidcov{
	$slide = int($readsize/$cov);	## Determined by coverage and read size. Truncate to integer value.
	if($type eq 'pe'){$slide = int($slide/2);}
	my $warning3 = "Fatal Error: Read size is too small\n"; die $warning2 unless($readsize >= 1 && $slide >= 1);
}
sub qscore{								## Take qscore, calculate appropriate offset, determine correct ASCII character for quality score
	if ($q64){$offset = $qscore + 64;}	## Print "Using Q64 format\n";
	else{$offset = $qscore + 33;}		## Print "Using Q33 format\n";
	$score = $ASCII{$offset};			## The artificial quality score to be assigned in the FastQ output
}
sub stats{
	print "\nInput File:\t$fasta\n";
	if ($type eq 'se'){print "Output File:\t$fasta.$readsize.SE.fastq\n";}
	else {print "Output Files:\t$fasta.$readsize.R1.fastq + $fasta.$readsize.R2.fastq\n";}
	print "Coverage:\t$cov"."X\n";
	print "Quality score:\t$qscore\n";
	if ($q64){print "Format:\t\tQ64 Illumina <= 1.7 ; Base quality set to $qscore \(ASCII character = $score\)\n";}
	else{print "Format:\t\tQ33 Sanger/Illumina 1.8+; Base quality set to $qscore \(ASCII character = $score\)\n";}
	print "Read size:\t$readsize nt\n";
	print "Type:\t\t";
	if ($type eq 'se'){print "Single ends\n";}
	elsif ($type eq 'pe'){print "Paired ends\n";}
	if ($mode eq 'slide'){print "Mode:\t\tSliding windows\nSlide every:\t$slide nucleotides\n";}
	elsif ($mode eq 'rand'){print "Mode:\t\tRandom selection\n";}
	print "\n\n";
}
sub se{ ## Print SE FastQ output:
	$count++;
	print SE '@SYNTHREAD_'."$count\n";
	print SE "$read\n";
	print SE '+'."\n";
	print SE "$score" x $readsize;   ## assign fake quality score
	print SE "\n";
	$count++;
	print SE '@SYNTHREAD_'."$count\n";
	print SE "$rev_read\n";
	print SE '+'."\n";
	print SE "$score" x $readsize;
	print SE "\n";
}
sub pe{ ## Print PE FastQ output:
	$count++;
	print R1 '@SYNTHREAD_'."$count".'_R1'."\n";
	print R1 "$pe1\n";
	print R1 '+'."\n";
	print R1 "$score" x $readsize;   ## assign fake quality score
	print R1 "\n";
	print R2 '@SYNTHREAD_'."$count".'_R2'."\n";
	print R2 "$pe2\n";
	print R2 '+'."\n";
	print R2 "$score" x $readsize;
	print R2 "\n";
	$count++;
	print R1 '@SYNTHREAD_'."$count".'_R1'."\n";
	print R1 "$rev1\n";
	print R1 '+'."\n";
	print R1 "$score" x $readsize;   ## assign fake quality score
	print R1 "\n";
	print R2 '@SYNTHREAD_'."$count".'_R2'."\n";
	print R2 "$rev2\n";
	print R2 '+'."\n";
	print R2 "$score" x $readsize;
	print R2 "\n";
}