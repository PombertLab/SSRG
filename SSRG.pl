#!/usr/bin/env perl
## Pombert Lab, Illinois Tech

my $version = '2.0a';
my $name = 'SSRG.pl';
my $updated = '2024-05-28';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use File::Path qw(make_path);
use PerlIO::gzip;

#########################################################################
### Command line options
#########################################################################

my $usage = <<"END_OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    SSRG (Synthetic Short Read Generator) generates short reads in FASTQ
            format (Q33) from multifasta files

COMMAND     ${name} \\
              -f *.fasta \\
              -o FASTQ \\
              -g \\
              -r 250 \\
              -i 350 \\
              -s 20 \\
              -c 25 \\
              -t PE

OPTIONS:
-f (--fasta)       Fasta/multifasta file(s)
-o (--outdir)      Output directory [Default: ./]
-g (--gzip)        Compress output FASTQ files with GZIP
-l (--list)        List of fasta file(s), one per line
-r (--readsize)    Desired read size [default: 150]
-t (--type)        Read type: single (SE) or paired ends (PE) [default: PE]
-i (--insert)      PE insert size [default: 350]
-s (--sdev)        PE insert size standard deviation (in percentage) [default: 10]
-c (--coverage)    Desired sequencing depth [default: 50]
-qs (--qscore)     Quality score associated with each base [default: 30]
-q64               Use the old Illumina Q64 FastQ format instead of the default Q33 Sanger/Illumina 1.8+ encoding
-v (--version)     Show script version
END_OPTIONS

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my @commands = @ARGV;

## Declare options
my @fasta;
my $outdir = './';
my $gzip;
my $list;
my $readsize = 150;
my $type = 'pe';
my $insert = 350;
my $sdev = 10;
my $slide;
my $cov = '50';
my $qscore = 30;
my $q64 = ''; ## Default is false so that Q33 format is used
my $sc_version;
GetOptions(
    'f|fasta=s@{1,}' => \@fasta,
    'o|outdir=s' => \$outdir,
    'g|gzip' => \$gzip,
    'l|list=s' => \$list,
    'r|readsize=i' => \$readsize,
    'i|insert=i' => \$insert,
    's|sdev=i' => \$sdev,
    't|type=s' => \$type,
    'qscore|qs=i' => \$qscore,
    'q64' => \$q64,
    'c|coverage=s' => \$cov,
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
### Command line checks
#########################################################################

my $warning1 = "\nFatal Error: Read size must be an integer\n";
my $warning2 = "\nFatal Error: Insert size \($insert\) is smaller than read size \($readsize\). Please use a larger insert size.\n\n";

die $warning1 unless (int($readsize) == $readsize);
if ($type eq 'pe'){
    die $warning2 if ($insert < $readsize);
}

#########################################################################
### Log
#########################################################################

my $start = localtime(); my $tstart = time;
my $date = `date`;
my $logfile = 'SSRG.log';

open LOG, '>>', $logfile or die "Can't create file $logfile: $!\n";
print LOG "\# $date";
print LOG "COMMAND LINE: $name @commands\n\n";

#########################################################################
### Checking ASCII offset for quality score
#########################################################################

my %ASCII = ASCII();
my $offset;
my $score;
qscore();

#########################################################################
### Output directory
#########################################################################

unless (-d $outdir){
    make_path( $outdir, { mode => 0755 } ) or die "Can't create folder $outdir: $!\n";
}

#########################################################################
### Working on files
#########################################################################

my $fasta;        ## Init fasta file
my $count = 0;    ## Read number counter to be auto-incremented
my $read;
my $rev_read;
my $pe1;
my $pe2;
my $rev1;
my $rev2;

if ($list){
    open LIST, '<', $list or die "Can't open $list: $!\n";
    print "\nSequences listed in $list:\n";
    while (my $line = <LIST>){
        chomp $line;
        my $fasta = `printf $line`; ## Converting shortcuts (e.g. ~, $HOME) to absolute paths
        print "$fasta\n";
        push (@fasta, $fasta);
    }
}

while ($fasta = shift @fasta) {

    open IN, '<', $fasta or die "Can't open $fasta: $!\n";
    $fasta =~ s/\.\w+$//; ## Removing file extensions
    $type = lc($type);

    my $basename = fileparse($fasta);
    my $filename = "${outdir}/$basename.$readsize";

    if($type eq 'se'){
        if ($gzip){
            my $se_gz = "$filename.SE.fastq.gz";
            open SE, ">:gzip", $se_gz or die "Can't create file $se_gz: $!\n";
        }
        else{
            my $se_fq = "$filename.SE.fastq";
            open SE, ">", $se_fq or die "Can't create file $se_fq: $!\n";
        }
    }
    elsif($type eq 'pe'){
        if ($gzip){
            my $r1_gz = "$filename.R1.fastq.gz";
            my $r2_gz = "$filename.R2.fastq.gz";
            open R1, ">:gzip", $r1_gz or die "Can't create file $r1_gz: $!\n";
            open R2, ">:gzip", $r2_gz or die "Can't create file $r2_gz: $!\n";
        }
        else{
            my $r1_fq = "$filename.R1.fastq";
            my $r2_fq = "$filename.R2.fastq";
            open R1, ">", $r1_fq or die "Can't create file $r1_fq: $!\n";
            open R2, ">", $r2_fq or die "Can't create file $r2_fq: $!\n";
        }
    }
    stats();

    my @contigs = ();    ## Initializing list of contigs per fasta file
    my %contigs = ();    ## Initializing hash of contigs per fasta file: key = contig name; value = sequence
    my $name = undef;    ## Contig name to be memorized

    while (my $line = <IN>){    ## Creating list of contigs and hash of sequences
        chomp $line;
        if ($line =~ /^>(.*)$/){
            $name = $1;
            push(@contigs, $name);
            $contigs{$name} = undef;
        }
        else{
			$contigs{$name} .= $line;
		}
    }

    while (my $todo = shift@contigs){ ## Parse each contig in the forward and reverse directions 

        chomp $todo;
        my $len = length($contigs{$todo});
        my $seq = $contigs{$todo};
        my $rev = reverse($seq);    ## get the reverse complement of the contig
        $rev =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;

        if ($type eq 'se'){

            my $numreads = int(($len*$cov)/$readsize);    ## Number of reads to achieve desired depth
            my $stdreads = int($numreads/2);            ## Number of reads per strand +/-

            for (1..$stdreads){
                my $max = $len - $readsize -1;
                my $plus = int(rand($max));
                my $minus = int(rand($max));
                $read = substr($seq, $plus, $readsize);
                $rev_read = substr($rev, $minus, $readsize);
                se();
            }
        }

        elsif ($type eq 'pe'){

            my $numreads = int(($len*$cov)/$readsize);    ## Number of reads to achieve desired depth
            my $stdreads = int($numreads/4);            ## Number of reads per strand +/-
            my $gap = $insert - ($readsize*2);
            my $dev = int($insert/$sdev);

            for (1..$stdreads){

                my $max = $len - $insert -$dev -1;
                my $rand = int(rand($dev) - ($dev/2));
                my $plus = int(rand($max));
                my $minus = int(rand($max));

                $pe1 = substr($seq, $plus, $readsize);
                $pe2 = substr($seq, $plus+$gap+$rand+$readsize, $readsize);
                $pe2 = reverse($pe2);
                $pe2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;
                $rev1 = substr($rev, $minus, $readsize);
                $rev2 = substr($rev, $minus+$gap+$rand+$readsize, $readsize);
                $rev2 = reverse($rev2);
                $rev2 =~ tr/ATGCNRYSWKMBVDHatgcnryswkmbvdh/TACGNYRWSMKVBHDtacgnyrwsmkvbhd/;

                pe();

            }
        }
    }
}

my $end = localtime();
my $time_taken = time - $tstart;
print "\nSSRG started on: $start\nSSRG ended on: $end\n";
print "Time elapsed: $time_taken seconds\n";

close IN; 

#########################################################################
### Closing out
#########################################################################

if ($type eq 'se'){
    if ($gzip eq ':gzip'){ binmode SE, ":gzip(none)"; }
    close SE;
}
elsif ($type eq 'pe'){
    if ($gzip eq ':gzip'){ 
        binmode R1, ":gzip(none)";
        binmode R2, ":gzip(none)";
    }
    close R1;
    close R2;
}
exit(0);

#########################################################################
### Subroutines
#########################################################################

sub ASCII{    ## If Q33 Q=0 is ASCII 33; if Q64 Q=0 is ASCII 64
    %ASCII = (
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
}

sub qscore{                              ## Take qscore, calculate appropriate offset, determine correct ASCII character for quality score
    if ($q64){                           ## Print "Using Q64 format\n";
		$offset = $qscore + 64;
	}
    else{                                ## Print "Using Q33 format\n";
		$offset = $qscore + 33;
	}
    $score = $ASCII{$offset};            ## The artificial quality score to be assigned in the FastQ output
}

sub stats{

    print "\nInput File:\t$fasta\n";

    if ($type eq 'se'){
        if ($gzip){
            print "Output File:\t$fasta.$readsize.SE.fastq.gz\n";
        }
        else {
            print "Output File:\t$fasta.$readsize.SE.fastq\n";
        }
    }
    else {
        if ($gzip){
            print "Output Files:\t$fasta.$readsize.R1.fastq.gz + $fasta.$readsize.R2.fastq.gz\n";
        }
        else{
            print "Output Files:\t$fasta.$readsize.R1.fastq + $fasta.$readsize.R2.fastq\n";
        }
    }

    print "Coverage:\t$cov"."X\n";
    print "Quality score:\t$qscore\n";

    if ($q64){
        print "Format:\t\tQ64 Illumina <= 1.7 ; Base quality set to $qscore \(ASCII character = $score\)\n";
    }
    else{
        print "Format:\t\tQ33 Sanger/Illumina 1.8+; Base quality set to $qscore \(ASCII character = $score\)\n";
    }

    print "Read size:\t$readsize nt\n";
    print "Type:\t\t";

    if ($type eq 'se'){
		print "Single ends\n";
	}
    elsif ($type eq 'pe') {
		print "Paired ends\n";
	}
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
    print R1 '@SYNTHREAD_'."$count".' 1'."\n";
    print R1 "$pe1\n";
    print R1 '+'."\n";
    print R1 "$score" x $readsize;   ## assign fake quality score
    print R1 "\n";
    print R2 '@SYNTHREAD_'."$count".' 2'."\n";
    print R2 "$pe2\n";
    print R2 '+'."\n";
    print R2 "$score" x $readsize;
    print R2 "\n";

    $count++;
    print R1 '@SYNTHREAD_'."$count".' 1'."\n";
    print R1 "$rev1\n";
    print R1 '+'."\n";
    print R1 "$score" x $readsize;   ## assign fake quality score
    print R1 "\n";
    print R2 '@SYNTHREAD_'."$count".' 2'."\n";
    print R2 "$rev2\n";
    print R2 '+'."\n";
    print R2 "$score" x $readsize;
    print R2 "\n";

}