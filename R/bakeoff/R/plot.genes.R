plot.genes = function(
	genes,
	region,
	aesthetic = list(
		heights = c(
			gene = 0.4,
			exon = 0.15,
			arrow = 0.35,
			label = 0.7
		),
		colour = c(
			gene = 'black',
			exon = 'grey',
			arrow = 'black'
		)
	),
	verbose = FALSE
) {
	genes$layout_level = NA
	genes = genes[ order( genes$chromosome, genes$start ), ]
	wGene = which( genes$feature == 'gene' & genes$chromosome == region$chromosome )
	wExon = which( genes$feature == 'exon' & genes$chromosome == region$chromosome )
	wCDS = which( genes$feature == 'CDS' & genes$chromosome == region$chromosome )
	genes$layout_level[ wGene ] =  layout.intervals( genes[ wGene, ], spacer = (region$end - region$start) / 10 )
	genes$layout_level[ wExon ] = genes[wGene,]$layout_level[ match( genes$ID[wExon], genes$ID[wGene] )]
	if( verbose ) {
		cat( "GENES:\n" )
		print( genes[wGene,] )
		cat( "EXONS:\n" )
		print( genes[wExon,] )
	}
	blank.plot(
		xlim = c( region$start, region$end ), ylim = c( -0.1, max( max( genes$layout_level[wGene] ) + 0.1, 2.5 ) ),
		xlab = spr( "Position on chromosome %s", region$chromosome )
	)
	# plot a line for the gene
	segments( x0 = genes$start[wGene], x1 = genes$end[wGene], y0 = genes$layout_level[wGene], y1 = genes$layout_level[wGene] )

	min.size = (region$end - region$start) / 100

	plot.arrows <- function( genes, arrow.length ) {
		if( nrow(genes) == 0 ) {
			return ;
		}
		pos = genes$start
		pos[ which(genes$strand == "-") ] = genes$end[ which(genes$strand == "-") ]
		sign = rep( 1, nrow( genes ))
		sign[ which(genes$strand == "-") ] = -1
		segments(
			x0 = c( pos, pos, pos + sign * ( arrow.length - min.size/4), pos + sign * ( arrow.length - min.size/4)),
			x1 = c( pos, pos + sign * arrow.length, pos + sign * arrow.length, pos + sign * arrow.length ),
			y0 = c(
				genes$layout_level - aesthetic$height['arrow'],
				genes$layout_level + aesthetic$height['arrow'],
				genes$layout_level + aesthetic$height['arrow'] - 0.1,
				genes$layout_level + aesthetic$height['arrow'] + 0.1
			),
			y1 = rep( genes$layout_level + aesthetic$height['arrow'], 4 ),
			col = aesthetic$colour['arrow']
		)
	}
	plot.exons <- function( genes, exons, min.size ) {
		w = which( genes$end - genes$start > min.size )
		if( length(w) > 0 ) {
			ids = genes$ID[w]
			wPlotExons = which( exons$ID %in% ids )
			rect(
				xleft = exons$start[wPlotExons],
				xright = exons$end[wPlotExons],
				ybottom = exons$layout_level[wPlotExons] - aesthetic$height['exon'],
				ytop = exons$layout_level[wPlotExons] + aesthetic$height['exon'],
				col = aesthetic$colour['exon'],
				border = NA
			)
		}
	}

	plot.arrows( genes[wGene,], min.size )
	plot.exons( genes[wGene,], genes[wExon,], min.size )
	
	text(
		genes$end[wGene] + min.size,
		genes$layout_level[wGene],
		genes$symbol[wGene],
		font = 3,
		adj = c( 0, 0.5 ),
		cex = aesthetic$height['label']
	)
}

