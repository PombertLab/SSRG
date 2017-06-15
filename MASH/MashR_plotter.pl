#!/usr/bin/perl
# Pombert lab 2017
# Plot Mash distances from CSV files generated by MashToDistanceCSV.pl
# Requires R, R-devel
# Requires the Rtsne, gplots, ggplot2, ggfortify, cluster, plotly (for MDS), RColorBrewer, and ape R packages
# E.g. to install the Rtsne package in R: type install.packages("Rtsne")
# Version 0.4d

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

## Usage definition
my $usage = "\nUSAGE = MashR_plotter.pl [OPTIONS]\n
EXAMPLE (simple): MashR_plotter.pl -i Mash.mashdist.csv -o image_01
EXAMPLE (advanced): MashR_plotter.pl -type cluster -m tsne -i Mash.mashdist.csv -R Mash.R -t pdf -o image_01 -pe 30 -it 500";
my $hint = "Type MashR_plotter.pl -h (--help) for list of options\n";
die "$usage\n\n$hint\n" unless@ARGV;

## Defining options
my $options = <<'END_OPTIONS';

OPTIONS:
--help (-h)		Display this list of options
--type (-t)		Plot type: cluster, tree, heatmap [Default: cluster]
--method (-m)		Dimensionality reduction method for clusters (mds, tsne) [Default: mds]
--input (-i)		Input file [Default: Mash.mashdist.csv]
--rscript (-R)		R script output name generated [Default: Mash.R]
--format (-f)		Image format (pdf, svg, jpeg, png) [Default: pdf]
--output (-o)		Output plot name [Default: plot]
--resolution (-res)	Resolution (in DPI) [Default: 300]
--labels (-lb)		Displays labels

## Clustering options
--cluster (-cl)		Clustering method (pam, fanny, kmeans, clara) [Default: pam]
--nclust (-nc)		Number of clusters desired [Default: 10]

## t-SNE options
--perplexity (-pe)	Perplexity [Default: 30]
--iterations (-it)	Maximum number of iterations [Default: 500]
--dimensions (-di)	Number of dimensions [Default: 2]
--cmode (cm)		Color mode: rainbow, heat, terrain, topo, cm or none [Default: rainbow]

## Phylogenetic tree options
--treetype (-tt)	Tree type: phylogram, cladogram, fan, unrooted or radial [Default: phylogram]
--distmeth (-dm)	Distance method: nj (neighbor-joining) or upgma [Default: nj]
--outgroup (-og)	Desired outgroup from the distance matrix (e.g. -og outgroup_name)
--newick (-nw)		Phylogenetic tree ouput in Newick format [Default: tree.tre]

## Heatmap options
--colors		Desired library(RColorBrewer) colors in order of decreasing similarity [Default: white yellow red]
--shades		Number of color shades desired [Default: 300]
--separator (-sep)	Switch off row/column separators
--clcol			Switch off column clustering
--clrow			Switch off row clustering

## R plotter options
--width (-wd)		Plot width [Default: 16]
--height (-he)		Plot height [Default: 10]
--fonts			Fonts [Default: Times]
--fontsize (-fs)	Plot font size [Default: 16]
--symbol (-pch)		R Plot PCH Symbols (for cluster graphs) [Default: 1]
--edges (-ed)		Draw edges (for cluster graphs)
--xrange (-x)		X-axis width (for cluster graphs) [Default: 0.005]

END_OPTIONS

my $help ='';
my $type = 'cluster';
my $tt = 'phylogram'; my $dm = 'nj'; my $nw = 'tree.tre';
my $wd = '16'; my $he = '10'; my $fs = '16'; ## Width, height and font size
my $res = 300;
my $method = 'mds';
my $cluster = 'pam';
my $labels;
my $nclust = 10;
my $cmode = 'rainbow';
my $og;
my $input = 'Mash.mashdist.csv';
my $rscript = 'Mash.R';
my $format = 'pdf';
my $perplexity = '30';
my $iterations = '500';
my $dimensions = '2';
my $output = "plot";
my $symbol = '1';
my $edge = '';
my $fonts = "Times";
my $xrange = '0.005';
my @colors = ();
my $shades = '300';
my $sep;
my $clcol; my $clrow;

GetOptions(
	'h|help' => \$help,
	't|type=s' => \$type,
	'tt|treetype=s' => \$tt,
	'dm|distmeth=s' => \$dm,
	'og|outgroup=s' => \$og,
	'nw|newick=s' => \$nw,
	'wd=i' => \$wd,
	'he=i' => \$he,
	'resolution|res=i' => \$res,
	'fs|fontsize=i' => \$fs,
	'method|m=s' => \$method,
	'cluster|cl=s' => \$cluster,
	'nclust|nc=i' => \$nclust,
	'cmode|cm=s' => \$cmode,
	'labels|lb'	=> \$labels,
	'input|i=s' => \$input,
	'rscript|R=s' => \$rscript,
	'format|f=s' => \$format,
	'perplexity|pe=i' => \$perplexity,
	'iterations|it=i' => \$iterations,
	'dimensions|di=i' => \$dimensions,
	'output|o=s' => \$output,
	'symbol|pch=i' => \$symbol,
	'edges|ed' => \$edge,
	'fonts=s' => \$fonts,
	'xrange|x=s' => \$xrange,
	'colors=s@{1,}' => \@colors,
	'shades=i' => \$shades,
	'sep|separator' => \$sep,
	'clcol' => \$clcol,
	'clrow' => \$clrow,
);

if (scalar @colors == 0){@colors = ('white', 'yellow', 'red');}
if ($help){die "$usage\n$options";}

if ($type eq 'cluster'){
	Rhm();
	Format();
	if ($method eq 'mds'){
		Rnames();
		print OUT 'library(ggplot2)'."\n".'library(ggfortify)'."\n".'library(cluster)'."\n".'library(plotly)'."\n";
		if(($cluster eq 'pam')||($cluster eq 'fanny')||($cluster eq 'clara')){
			print OUT 'p <- autoplot('."$cluster".'(mat_data, '."$nclust".'), frame = T, frame.type = \'t\', main = \'Genetic Similarity\'';
			if ($labels){print OUT ", label = as.integer(1)";}
			print OUT ")\n";
		}
		elsif($cluster eq 'kmeans'){
			print OUT 'p <- autoplot('."$cluster".'(mat_data, '."$nclust".'), data = mat_data, frame = T, frame.type = \'t\', main = \'Genetic Similarity\'';
			if ($labels){print OUT ", label = as.integer(1)";}
			print OUT ")\n";
		}
		print OUT "print(p)\n";
	}
	elsif ($method eq 'tsne'){
		print OUT 'library(Rtsne)'."\n".'library(cluster)'."\n";
		Rnames();
		if ($cmode eq 'rainbow'){print OUT 'palette('."$cmode".'('."$nclust".'))'."\n";}
		elsif ($cmode eq 'none'){}
		else{print OUT 'palette('."$cmode".'.colors('."$nclust".'))'."\n";}
		print OUT 'tsne <- Rtsne(distance_matrix[,-1], dims = '."$dimensions".', perplexity='."$perplexity".', check_duplicates=FALSE, is_distance=TRUE, verbose=TRUE, max_iter = '."$iterations".')$Y'."\n";
		if ($edge){$edge = ' t="o",';}
		my $id = '$clustering'; if ($cluster eq 'kmeans'){$id = '$cluster'}
		if ($cmode eq 'none'){
			print OUT 'plot(tsne, main="tsne",'."$edge".' pch='."$symbol".')'."\n";
			if ($labels){print OUT 'text(tsne, cex=0.4, labels=distance_matrix$OTU)'."\n";}
		}
		else{
			print OUT 'plot(tsne, main="tsne",'."$edge".' col='."$cluster".'(mat_data, '."$nclust".')'."$id".', pch='."$symbol".')'."\n";
			if ($labels){print OUT 'text(tsne, cex=0.4, labels=distance_matrix$OTU, col='."$cluster".'(mat_data, '."$nclust".')'."$id".')'."\n";}
		}
	}
	close IN;
	close OUT;
}
elsif ($type eq 'tree'){
	Rhm();
	Format();
	print OUT 'row.names(distance_matrix) <- distance_matrix[, 1]'."\n";
	print OUT 'distance_matrix <- distance_matrix[, -1]'."\n";
	print OUT 'library(ape)'."\n";
	if ($og){
		Dist();
		print OUT 'rooted <- root(tree, "'."$og".'", node = NULL, resolve.root = TRUE)'."\n";
		print OUT 'plot(rooted, "'."$tt".'")'."\n";
		print OUT 'write.tree(rooted, file="'."$nw".'")'."\n";
	}
	else {
		Dist();
		print OUT 'plot(tree, "'."$tt".'")'."\n";
		print OUT 'write.tree(tree, file="'."$nw".'")'."\n";
	}
	close IN;
	close OUT;
}
elsif ($type eq 'heatmap'){
	Rhm();
	Format();
	Rnames();
	print OUT 'library(gplots)'."\n".'library(RColorBrewer)'."\n";
	print OUT 'colors <- colorRampPalette(c(';
	for (0..$#colors-1){print OUT '"'."$colors[$_]\", ";}
	print OUT "\"$colors[$#colors]\"".'))(n = '."$shades".')'."\n";
	print OUT 'heatmap.2(mat_data, density.info="none", trace="none", col=colors, margins =c('."$he".','."$wd".')';
	unless ($sep) {print OUT ', rowsep = 1:nrow(mat_data), colsep = 1:ncol(mat_data)';}
	if (($clcol) || ($clrow)){
		if (($clcol) && ($clrow)){print OUT ', Colv="NA", Rowv="NA", dendrogram="none"';}
		elsif ($clcol){print OUT ', Colv="NA", dendrogram="row"';}
		elsif ($clrow){print OUT ', Rowv="NA", dendrogram="col"';}
	}
	else {print OUT ', dendrogram="both"';}
	print OUT ')'."\n";
	print OUT 'dev.off()'."\n";
	close IN;
	close OUT;
}

## Running R script
system "chmod a+x $rscript\; ./$rscript";

## Subroutines
sub Format {
	my $wide = 3300; my $tall = 3300;
	my $scaling = ($res/300);
	$wide = int($wide*$scaling); $tall = int($tall*$scaling*($he/$wd)); print "saling = $scaling,  wide = $wide, tall = $tall\n";
	if ($format eq 'pdf'){print OUT "$format".'(file="'."$output\.$format".'", useDingbats=FALSE, family="'."$fonts".'", pointsize='."$fs".', width='."$wd".',height='."$he".')'."\n";}
	elsif ($format eq 'svg'){print OUT "$format".'(file="'."$output\.$format".'", width='."$wd".', height='."$he".')'."\n";}
	elsif (($format eq 'jpeg')||($format eq 'png')){print OUT "$format".'(file="'."$output\.$format".'", width='."$wide".', height='."$tall".', res = '."$res".')'."\n";}
}

sub Rhm { ## Generating R script, headers and matrix
	open IN, "<$input" or die "Cannot open $input\n";
	open OUT, ">$rscript";
	print OUT '#!/usr/bin/Rscript'."\n";
	print OUT 'message ("'."\nPlotting $input\:\n".'")'."\n";
	print OUT 'distance_matrix <- read.csv("'."$input".'")'."\n";
}

sub Rnames {
	print OUT 'rnames <- distance_matrix[,1]'."\n";
	print OUT 'mat_data <- data.matrix(distance_matrix[,2:ncol(distance_matrix)])'."\n";
	print OUT 'rownames(mat_data) <- rnames'."\n";
}

sub Dist{
	if ($dm eq 'nj') {print OUT 'tree <- nj(as.dist(distance_matrix))'."\n";}
	elsif ($dm eq 'upgma') {print OUT 'tree <- upgma(as.dist(distance_matrix))'."\n";}
}