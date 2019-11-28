plot.ucsc.genes <- function(
    genes,
    region,
    height_in_inches = 1,
    exons = get.exons( genes ),
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
) {
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
	chromosome = region$chromosome
	region = c( start = region$start, end = region$end )
	w = which( genes$chromosome == chromosome & genes$txEnd >= region[1] & genes$txStart <= region[2] )
	if( verbose ) {
		print(w)
	}
    if( length(w) > 0 ) {
		genes = genes[w,]
		genes$lty = 1
		stopifnot( nrow( genes ) > 0 )
        genes = genes[ order( genes$txStart ),, drop = FALSE ]
        genes$y = NA ;
        genes$y[1] = 1 ;
        if( nrow( genes ) > 1 ) {
			if( is.null( spacer )) {
				spacer = ( region[2] - region[1] ) / 10 ;
			}
            maxes = ( genes[1,]$txEnd + spacer )
            for( i in 2:nrow( genes )) {
                for( l in 1:length( maxes )) {
                    if( genes$txStart[i] >= maxes[l] ) {
                        genes$y[i] = l ;
                        maxes[l] = genes$txEnd[i] + spacer ;
                        break ;
                    }
                }
                if( is.na( genes$y[i] )) {
                    maxes = c( maxes, genes$txEnd[i] + spacer )
                    genes$y[i] = length( maxes ) ;
                }
            }
        }
		
		if( is.null( label.cex )) {
			# We fit about 7 gene names per inch at cex = 1.  After that we need to start scaling.
			label.cex = min( 1, ( 6 * height_in_inches ) / max( genes$y, na.rm = T ) )
		}
		genes$label.cex = label.cex
		wPseudoGene = which( genes$cdsStart == genes$cdsEnd )
		if( length( wPseudoGene ) > 0 ) {
			genes$lty[ wPseudoGene ] = 2
			genes$label.cex[ wPseudoGene ] = label.cex * 0.8
		}
        
        exons$y = genes$y[ match( exons$name, genes$name )]
		if( verbose ) {
			print( exons )
		}
        height_of_genes = height_in_inches / max( genes$y )
        
        strands = c("+","-" ) ;
        level = c( 80, 70 ) ;
		if( vertical ) {
			plot( region[1], 0, pch = '', ylim = region, xlim = c( 0, max( 2, max(genes$y)+1 ) ), xlab = '', ylab = '', xaxt = 'n', ... ) ;
		} else {
			plot( 0, region[1], pch = '', xlim = region, ylim = c( 0, max( 2, max(genes$y)+0.5 ) ), xlab = '', ylab = '', yaxt = 'n', ... ) ;
		}

        arrow.sep = ( region[2] - region[1] ) / 200 ;
		arrow.voffset = 0
		
        relative.lengths = ( genes$txEnd - genes$txStart ) / ( region[2] - region[1] )

        genes$mark1 = ( 0.25 * genes$txStart + 0.75 * genes$txEnd ) ;
        genes$mark2 = ( 0.5 * genes$txStart + 0.5 * genes$txEnd ) ;
        genes$mark3 = ( 0.75 * genes$txStart + 0.25 * genes$txEnd ) ;
        genes$sign = 1 ;
        genes$sign[ which( genes$strand == '-' ) ] = -1 ;

	gene.height = 0.4
	exon.height = 0.25

	if( vertical ) {
	segments(
	    y0 = genes$txStart, y1 = genes$txEnd,
	    x0 = genes$y, x1 = genes$y,
	    col = colours[['gene']],
			lty = genes$lty
	)
	} else {
		segments(
		    x0 = genes$txStart, x1 = genes$txEnd,
		    y0 = genes$y, y1 = genes$y,
		    col = colours[['gene']],
				lty = genes$lty
		)
	}
	w = which( genes$sign == 1 ) 
	arrow.top = exon.height + 0.1
	if( length(w) > 0 ) {
		segment.length = (region[2] - region[1])/120 ;
		if( vertical ) {
			segments( y0 = genes$txStart[w], y1 = genes$txStart[w], x0 = genes$y[w] - exon.height, x1 = genes$y[w] + arrow.top, col = colours[['gene']] )
			segments(
			    y0 = c( genes$txStart[w], genes$txStart[w], genes$txStart[w] + segment.length/2, genes$txStart[w] + segment.length/2 ),
			    y1 = c( genes$txStart[w], genes$txStart[w] + segment.length, genes$txStart[w] + segment.length, genes$txStart[w] + segment.length ),
			    x0 = c( genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top + 0.05, genes$y[w] + arrow.top - 0.05 ),
			    x1 = c( genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top ),
			    col = colours[['gene']],
				xpd = NA
			)
		} else {
			segments( x0 = genes$txStart[w], x1 = genes$txStart[w], y0 = genes$y[w] - exon.height, y1 = genes$y[w] + arrow.top, col = colours[['gene']] )
			segments(
			    x0 = c( genes$txStart[w], genes$txStart[w], genes$txStart[w] + segment.length/2, genes$txStart[w] + segment.length/2 ),
			    x1 = c( genes$txStart[w], genes$txStart[w] + segment.length, genes$txStart[w] + segment.length, genes$txStart[w] + segment.length ),
			    y0 = c( genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top + 0.05, genes$y[w] + arrow.top - 0.05 ),
			    y1 = c( genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top, genes$y[w] + arrow.top ),
			    col = colours[['gene']],
						xpd = NA
			)
		}
	}
	
	wBigEnough = which( ( genes$txEnd - genes$txStart ) > arrow.sep * 2 ) ;
	if( length(wBigEnough) > 0 ) {
		arrows = data.frame()
		for( i in wBigEnough ) {
			arrows = rbind(
				arrows,
				data.frame(
					name2 = genes$name2[i],
					y = genes$y[i],
					x = seq( from = genes$txStart[i] + arrow.sep, to = genes$txEnd[i], by = arrow.sep ),
					sign = genes$sign[i]
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

			wExon = which( exons$name2 %in% genes$name2[wBigEnough] )
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

			wExon = which( exons$name2 %in% genes$name2[wBigEnough] )
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
			text( genes$y, genes$txEnd, label = genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, srt=90, col = colours[[ 'text' ]], font = 3 )
		} else {
			text( genes$txEnd, genes$y, label = genes$name2, adj = c( -0.1, 0.5 ), cex = label.cex, col = colours[[ 'text' ]], font = 3, xpd = NA )
		}
    } else {
#        plot.new()
        plot( region[1], 0, pch = '', xlim = region, ylim = c( 0, 1 ), xlab = '', ylab = 'genes' ) ;
    }
}
