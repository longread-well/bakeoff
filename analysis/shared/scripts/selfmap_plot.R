library( argparse )
library( RSQLite )

# Imported from shared_functions.R

fix_chromosome <- function( chromosome ) {
    chromosome = as.character( chromosome )
	chromosome = gsub( "chr", "", chromosome, fixed = T )
    n = nchar( chromosome )
    w = which( n == 1 ) ;
    chromosome[w] = sprintf( "0%s", chromosome[w] )
    return( chromosome )
}


parse_ranges <- function( range_spec, db = NULL ) {
    elts = strsplit( range_spec, split = ",", fixed = T )[[1]]
    ranges = data.frame()
    for( i in 1:length( elts )) {
        bits = strsplit( elts[i], split = ":", fixed = T )[[1]] ;
		well.formed = FALSE ;
		if( length( bits ) == 1 && !is.null( db )) {
			# Assume it's a gene name
			genePos = dbGetQuery( db, sprintf( "SELECT chrom, min( txStart )+1 AS start, max( txEnd )+1 AS end FROM refGene_20160609 WHERE name2 == '%s'", bits[1] ))
			if( nrow( genePos ) == 1 ) {
				genePos$chromosome = fix_chromosome( genePos$chrom )
				ranges = rbind(
					ranges,
					data.frame(
						chromosome = genePos$chromosome,
						start = genePos$start,
						end = genePos$end
					)
				)
			}
			well.formed = TRUE ;
		} else if( length( bits) == 2 ) {
			ranges = rbind(
				ranges,
				data.frame(
					chromosome = bits[1],
					start = as.integer( strsplit( bits[2], split = "-", fixed = T )[[1]][1] ),
					end = as.integer( strsplit( bits[2], split = "-", fixed = T )[[1]][2] )
				)
			)
			well.formed = TRUE ;
		}
		if( !well.formed ) {
            cat( sprintf( "Malformed range spec \"%s\".  Quitting.\n", elts[i] ))
            quit("no") ;
        }
    }
	ranges
}

plot.genes <- function(
    chromosome,
    region,
    local.genes,
    height_in_inches,
    exons = get.exons( local.genes ),
    vertical = FALSE,
	colours = list(
		gene = "blue",
		arrows = rgb( 0.5, 0.5, 0.5, 0.7 ),
		exon = "blue",
		text = "black"
	),
	verbose = FALSE,
	label.cex = NULL,
	spacer = NULL,
	...
)
{
	w = which( local.genes$chromosome == chromosome & local.genes$txEnd >= region[1] & local.genes$txStart <= region[2] )
	if( verbose ) {
		print(w)
	}
    if( length(w) > 0 ) {
		local.genes = local.genes[w,]
		local.genes$lty = 1
		stopifnot( nrow( local.genes ) > 0 )
        local.genes = local.genes[ order( local.genes$txStart ),, drop = FALSE ]
        local.genes$y = NA ;
        local.genes$y[1] = 1 ;
        if( nrow( local.genes ) > 1 ) {
			if( is.null( spacer )) {
				spacer = ( region[2] - region[1] ) / 10 ;
			}
            maxes = ( local.genes[1,]$txEnd + spacer )
            for( i in 2:nrow( local.genes )) {
                for( l in 1:length( maxes )) {
                    if( local.genes$txStart[i] >= maxes[l] ) {
                        local.genes$y[i] = l ;
                        maxes[l] = local.genes$txEnd[i] + spacer ;
                        break ;
                    }
                }
                if( is.na( local.genes$y[i] )) {
                    maxes = c( maxes, local.genes$txEnd[i] + spacer )
                    local.genes$y[i] = length( maxes ) ;
                }
            }
        }

		if( is.null( label.cex )) {
			# We fit about 7 gene names per inch at cex = 1.  After that we need to start scaling.
			label.cex = min( 1, ( 6 * height_in_inches ) / max( local.genes$y, na.rm = T ) )
		}
		local.genes$label.cex = label.cex
		wPseudoGene = which( local.genes$cdsStart == local.genes$cdsEnd )
		if( length( wPseudoGene ) > 0 ) {
			local.genes$lty[ wPseudoGene ] = 2
			local.genes$label.cex[ wPseudoGene ] = label.cex * 0.8
		}

        exons$y = local.genes$y[ match( exons$name, local.genes$name )]
		if( verbose ) {
			print( exons )
		}
        height_of_genes = height_in_inches / max( local.genes$y )

        strands = c("+","-" ) ;
        level = c( 80, 70 ) ;
		if( vertical ) {
			plot( region[1], 0, pch = '', ylim = region, xlim = c( 0, max( 2, max(local.genes$y)+1 ) ), xlab = '', ylab = '', xaxt = 'n', ... ) ;
		} else {
			plot( 0, region[1], pch = '', xlim = region, ylim = c( 0, max( 2, max(local.genes$y)+0.5 ) ), xlab = '', ylab = '', yaxt = 'n', ... ) ;
		}

        arrow.sep = ( region[2] - region[1] ) / 200 ;
		arrow.voffset = 0

        relative.lengths = ( local.genes$txEnd - local.genes$txStart ) / ( region[2] - region[1] )

        local.genes$mark1 = ( 0.25 * local.genes$txStart + 0.75 * local.genes$txEnd ) ;
        local.genes$mark2 = ( 0.5 * local.genes$txStart + 0.5 * local.genes$txEnd ) ;
        local.genes$mark3 = ( 0.75 * local.genes$txStart + 0.25 * local.genes$txEnd ) ;
        local.genes$sign = 1 ;
        local.genes$sign[ which( local.genes$strand == '-' ) ] = -1 ;

	gene.height = 0.4
	exon.height = 0.25

	if( vertical ) {
	segments(
	    y0 = local.genes$txStart, y1 = local.genes$txEnd,
	    x0 = local.genes$y, x1 = local.genes$y,
	    col = colours[['gene']],
			lty = local.genes$lty
	)
	} else {
		segments(
		    x0 = local.genes$txStart, x1 = local.genes$txEnd,
		    y0 = local.genes$y, y1 = local.genes$y,
		    col = colours[['gene']],
				lty = local.genes$lty
		)
	}
	w = which( local.genes$sign == 1 )
	arrow.top = exon.height + 0.1
	if( length(w) > 0 ) {
		segment.length = (region[2] - region[1])/120 ;
		if( vertical ) {
			segments( y0 = local.genes$txStart[w], y1 = local.genes$txStart[w], x0 = local.genes$y[w] - exon.height, x1 = local.genes$y[w] + arrow.top, col = colours[['gene']] )
			segments(
			    y0 = c( local.genes$txStart[w], local.genes$txStart[w], local.genes$txStart[w] + segment.length/2, local.genes$txStart[w] + segment.length/2 ),
			    y1 = c( local.genes$txStart[w], local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length ),
			    x0 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top + 0.05, local.genes$y[w] + arrow.top - 0.05 ),
			    x1 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top ),
			    col = colours[['gene']],
				xpd = NA
			)
		} else {
			segments( x0 = local.genes$txStart[w], x1 = local.genes$txStart[w], y0 = local.genes$y[w] - exon.height, y1 = local.genes$y[w] + arrow.top, col = colours[['gene']] )
			segments(
			    x0 = c( local.genes$txStart[w], local.genes$txStart[w], local.genes$txStart[w] + segment.length/2, local.genes$txStart[w] + segment.length/2 ),
			    x1 = c( local.genes$txStart[w], local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length, local.genes$txStart[w] + segment.length ),
			    y0 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top + 0.05, local.genes$y[w] + arrow.top - 0.05 ),
			    y1 = c( local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top, local.genes$y[w] + arrow.top ),
			    col = colours[['gene']],
						xpd = NA
			)
		}
	}

	wBigEnough = which( ( local.genes$txEnd - local.genes$txStart ) > arrow.sep * 2 ) ;
	if( length(wBigEnough) > 0 ) {
		arrows = data.frame()
		for( i in wBigEnough ) {
			arrows = rbind(
				arrows,
				data.frame(
					name2 = local.genes$name2[i],
					y = local.genes$y[i],
					x = seq( from = local.genes$txStart[i] + arrow.sep, to = local.genes$txEnd[i], by = arrow.sep ),
					sign = local.genes$sign[i]
				)
			)
		}

		if( vertical ) {
			segments(
				y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
				x0 = arrows$y + 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
				col = colours[['arrow']]
			)

			segments(
				y0 = arrows$x, y1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
				x0 = arrows$y - 0.2 + arrow.voffset, x1 = arrows$y + arrow.voffset,
				col = colours[['arrow']]
			)

			wExon = which( exons$name2 %in% local.genes$name2[wBigEnough] )
			if( length(wExon) > 0 ) {
				print( exons )
				exons$exonMid = (exons$exonStart + exons$exonEnd)/2
				segmentLength = (region[2] - region[1])/1000 ;
				exons$fixedExonStart = pmin( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
				exons$fixedExonEnd = pmax( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
				rect(
					ybottom = exons$fixedExonStart[wExon],
					ytop = exons$fixedExonEnd[wExon],
					xleft = exons$y[wExon] - 0.25,
					xright = exons$y[wExon] + 0.25,
					border = NA,
					col = colours[['exon']]
				)
			}
		} else {
			segments(
				x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
				y0 = arrows$y + 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
				col = colours[['arrow']]
			)

			segments(
				x0 = arrows$x, x1 = arrows$x + arrows$sign * 0.5 * arrow.sep,
				y0 = arrows$y - 0.2 + arrow.voffset, y1 = arrows$y + arrow.voffset,
				col = colours[['arrow']]
			)

			wExon = which( exons$name2 %in% local.genes$name2[wBigEnough] )
			if( length(wExon) > 0 ) {
				exons$exonMid = (exons$exonStart + exons$exonEnd)/2
				segmentLength = (region[2] - region[1])/1000 ;
				exons$fixedExonStart = pmin( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
				exons$fixedExonEnd = pmax( exons$exonStart, exons$exonMid - segmentLength, exons$exonMid + segmentLength ) ;
				rect(
					xleft = exons$fixedExonStart[wExon],
					xright = exons$fixedExonEnd[wExon],
					ybottom = exons$y[wExon] - 0.25,
					ytop = exons$y[wExon] + 0.25,
					border = NA,
					col = colours[['exon']]
				)
			}
		}
	}
		if( vertical ) {
			text( local.genes$y, local.genes$txEnd, label = local.genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, srt=90, col = colours[[ 'text' ]], font = 3 )
		} else {
			text( local.genes$txEnd, local.genes$y, label = local.genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, col = colours[[ 'text' ]], font = 3, xpd = NA )
		}
    } else {
#        plot.new()
        plot( region[1], 0, pch = '', xlim = region, ylim = c( 0, 1 ), xlab = '', ylab = 'genes' ) ;
    }
}

get.exons <- function( genes ) {
	result = data.frame()
	for( i in 1:nrow( genes )) {
		gene = genes[i,]
		if( gene$exonCount > 0 ) {
			result = rbind(
				result,
				data.frame(
					name = gene$name,
					name2 = gene$name2,
					txStart = gene$txStart,
					txEnd = gene$txEnd,
					exonStart = as.integer( strsplit( gene$exonStarts, split = "," )[[1]] ),
					exonEnd = as.integer( strsplit( gene$exonEnds, split = "," )[[1]] )
				)
			)
		}
	}
	return( result )
}

#############################
# Declare options
{
	parser <- ArgumentParser( description = 'Compare two genomes' )
	parser$add_argument(
		"--s1",
		type = "character",
		help = "First sequence file, in format name=<filename>",
		required = TRUE
	)
	parser$add_argument(
		"--chromosome1",
		type = "character",
		help = "Chromosome identifier of first sequence file (to match genes)",
		required = TRUE
	)
	parser$add_argument(
		"--genes1",
		type = "character",
		help = "First genes file",
		default = "/well/longread/projects/reference/GRCh37/UCSC/2019-08-27-GRCh37-wgEncodeGencodeBasicV19.tsv"
	)
	parser$add_argument(
		"--s2",
		type = "character",
		help = "Second sequence file, in format name=<filename>",
		required = TRUE
	)
	parser$add_argument(
		"--chromosome2",
		type = "character",
		help = "Chromosome identifier of second sequence file (to match genes)",
		required = TRUE
	)
	parser$add_argument(
		"--genes2",
		type = "character",
		help = "Second genes file",
		default = "/well/longread/projects/reference/GRCh38/UCSC/2019-08-27-GRCh38-wgEncodeGencodeBasicV31.tsv"
	)
	parser$add_argument(
		"--k",
		type = "integer",
		help = "k-mer size for upper panel",
		default = 100
	)
	parser$add_argument(
		"--range1",
		type = "character",
		help = "Range in sequence 1 to plot"
	)
	parser$add_argument(
		"--range2",
		type = "character",
		help = "Range in sequence 2 to plot"
	)
	parser$add_argument(
		"--output",
		type = "character",
		help = "Path to output file"
	)
}


opts = parser$parse_args()

cat( "Options are:\n" )
print( opts )



#############################
# Load data

load.genes = function( filename, condense = TRUE, protein.coding.only = TRUE ) {
	if( is.null( filename )) {
		return(
			data.frame(
				chrom = c(),
				txStart = c(),
				txEnd = c(),
				name2 = c(),
				chromosome = c()
			)
		)
	}
	# Load genes
	gene <- read.delim(
	filename,
		header=TRUE,
		as.is=TRUE
	);

	if( condense ) {
		gene <- gene[order(gene$txEnd - gene$txStart,decreasing=TRUE),];  #Get just longest transcript
		gene <- gene[ !duplicated( gene$name2 ), ];
	}
	gene <- gene[ !is.na(gene$txStart), ];

	if( protein.coding.only ) {
		gene <- gene[ gene$cdsStart != gene$cdsEnd, ]
	}

	chromosome =  gsub( "^chr", "", gene$chrom )
	w1 = which( nchar( chromosome ) == 1 )
	chromosome[ w1 ] = sprintf( "0%s", chromosome[w1] )
	gene$chromosome = chromosome
	return( gene ) ;
}

geneLists = list(
	'GRCh37' = "/well/longread/projects/reference/GRCh37/UCSC/2019-08-27-GRCh37-wgEncodeGencodeBasicV19.tsv",
	'GRCh38' = "/well/longread/projects/reference/GRCh38/UCSC/2019-08-27-GRCh38-wgEncodeGencodeBasicV31.tsv"
)


sequence1 = strsplit( opts$s1, split = '=' )[[1]]; names( sequence1 ) = c( "name", "filename" )
sequence2 = strsplit( opts$s2, split = '=' )[[1]]; names( sequence2 ) = c( "name", "filename" )

if( is.null( opts$genes1 )) {
	opts$genes1 = geneLists[[sequence1['name'] ]]
}
if( is.null( opts$gene2 )) {
	if( sequence2['name'] %in% names( geneLists ) ) {
		opts$genes2 = geneLists[[sequence2['name'] ]]
	}
}

cat( "Loading genes...\n" )
genes = list(
	first = load.genes( opts$genes1 ),
	second = load.genes( opts$genes2 )
)

run.selfmap <- function( sequence1, sequence2, k, chromosome1, chromosome2 ) {
	output = tempfile()
	cmd = sprintf( 'selfmap_v2.1-dev -sequence %s=%s %s=%s -kmer-size %d -o %s',
		sequence1['name'], sequence1['filename'],
		sequence2['name'], sequence2['filename'], k, output
	)
	cat( "Running \"", cmd, "\"...\n", sep = '' )
	system( cmd )
	cat( "Loading \"", output, "\"...\n", sep = '' )
	X = read.table( output, hea=T, as.is=T )
	X = X[ grep( sprintf( "^%s:", sequence1['name'] ), X$chromosome ), ]
	X = X[ grep( sprintf( "^%s:", sequence2['name'] ), X$other_chromosome ), ]
	X$chromosome = chromosome1
	X$other_chromosome = chromosome2
	cat( "...ok, read:\n" )
	print( head( X ))
	return(X)
}
cat( "Computing shared k-mers...\n" )
X1 = run.selfmap( sequence1, sequence2, opts$k, fix_chromosome( opts$chromosome1 ), fix_chromosome( opts$chromosome2 ))
X2 = run.selfmap( sequence1, sequence2, opts$k / 2, fix_chromosome( opts$chromosome1 ), fix_chromosome( opts$chromosome2 ))

if( !is.null( opts$range1 )) {
	range1 = parse_ranges( opts$range1 )
} else {
	range1 = data.frame(
		chromosome = opts$chromosome1,
		start = min( c( X1$position )),
		end = max( c( X1$position ) + opts$k )
	)
}

if( !is.null( opts$range2 )) {
	range2 = parse_ranges( opts$range2 )
} else {
	range2 = data.frame(
		chromosome = opts$chromosome1,
		start = min( c( X1$other_position )),
		end = max( c( X1$other_position ) + opts$k )
	)
}

myGrid <- function( main.divisions = 1, sub.divisions = 4,	main.col = "grey40", sub.col = "grey80", lty = 3, horiz = TRUE, vert = TRUE ) {
	xaxp = par( "xaxp" )
	yaxp = par( "yaxp" )
	if( horiz ) {
		abline( h = seq( from = yaxp[1], to = yaxp[2], length = (sub.divisions*yaxp[3]+1)), col = sub.col, lty = lty )
		abline( h = seq( from = yaxp[1], to = yaxp[2], length = (main.divisions*yaxp[3]+1)), col = main.col, lty = lty )
	}
	if( vert ) {
		abline( v = seq( from = xaxp[1], to = xaxp[2], length = (sub.divisions*xaxp[3]+1)), col = sub.col, lty = lty )
		abline( v = seq( from = xaxp[1], to = xaxp[2], length = (main.divisions*xaxp[3]+1)), col = main.col, lty = lty )
	}
}

#############################
# Plot


X1$col = 'black'
X1$col[ which( X1$other_orientation == '-' )] = 'red'
X1$pch = '.'

X2$col = 'grey'
X2$col[ which( X2$other_orientation == '-' )] = 'indianred2'


for( i in 1:1 ) {
	chromosomes = c(
		fix_chromosome( as.character( range1$chromosome[i] ) ),
		fix_chromosome( as.character( range2$chromosome[i] ) )
	)
	region_names = c(
		sprintf( "%s: %s:%d-%d", sequence1['name'], chromosomes[1], range1$start[i], range1$end[i] ),
		sprintf( "%s: %s:%d-%d", sequence2['name'], chromosomes[2], range2$start[i], range2$end[i] )
	)

	data = list(
		X2[
			which(
				X2$chromosome == chromosomes[1] & X2$position >= range1$start[1] & X2$position <= range1$end[i]
				& X2$other_chromosome == chromosomes[2] & X2$other_position >= range2$start[i] & X2$other_position <= range2$end[i]
			),
		],
		X1[
			which(
				X1$chromosome == chromosomes[1] & X1$position >= range1$start[1] & X1$position <= range1$end[i]
				& X1$other_chromosome == chromosomes[2] & X1$other_position >= range2$start[i] & X1$other_position <= range2$end[i]
			),
		]
	)
	names(data) = sprintf( "k=%d", c( as.integer(opts$k/2), as.integer(opts$k )))

	pdf( file = opts$output, width = 8, height = 8 )
	layout( matrix( c( 1, 2, 3, 4 ), nrow = 2, ncol = 2 ), heights = c( 5, 2 ), widths = c( 5.3, 2 ))
	par( mar = c( 2, 2, 2, 1 ))

	# Dot plot

	# Upper dot plot
	plot(
		0, 0, col = 'white',
		xlab = region_names[1],
		ylab = region_names[2],
		xaxt="n",
		xlim = c( range1$start, range1$end ),
		ylim = c( range2$start, range2$end )
	)
	for( j in 1:length( data )) {
		points(
			data[[j]]$position, data[[j]]$other_position,
			pch = '.', col = data[[j]]$col
		)
	}
	legend( "topleft", bty = 'n', legend = names( data ), pch = 19, col = sapply( 1:length(data), function(k) { data[[k]]$col[1] } ))
	axis( side = 3 )
	axis( side = 2 )
	axis( side = 1 )

	# Lower dot plot
	#points( XL$other_position, XL$position, pch = '.', col = XL$col )

	# grid
	#abline( a = 0, b = 1, col = "grey" )
	myGrid( main.col = "grey10", sub.col = "grey50" )


	# genes
	{
		par( mar = c( 2, 2, 0, 1 ), xaxt = 's', bty = 'o' )
		plot.genes( chromosomes[1], c( range1$start, range1$end ), genes[[1]], height_in_inches = 2 )
		myGrid( horiz = FALSE, main.col = "grey10", sub.col = "grey80" )

		par( mar = c( 2, 0, 2, 2 ), xaxt = 'n', yaxt = 'n', bty = 'o' )
		plot.genes( chromosomes[2], c( range2$start, range2$end ), genes[[2]], height_in_inches = 2, vertical = T )
		axis( side = 4 )
		myGrid( vert = FALSE, main.col = "grey10", sub.col = "grey80" )
	}

	dev.off()
}
