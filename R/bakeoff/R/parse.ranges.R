parse.ranges <- function( range_specs, db = NULL ) {
	ranges = data.frame()
	elts = range_specs
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
