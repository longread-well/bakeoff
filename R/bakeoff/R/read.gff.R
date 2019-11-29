read.gff <- function( filename ) {
	result = read.table(
		filename, head = F, comment = '#', sep = "\t", quote = "", stringsAsFactors = FALSE
	)
	# This is GFF https://en.wikipedia.org/wiki/General_feature_format
	colnames( result ) = c(
		"sequence",
		"source",
		"feature",
		"start",
		"end",
		"score",
		"strand",
		"phase",
		"attributes"
	)
	return( result )
}

