library( argparse )
library( hspf )

parse_arguments <- function() {
        parser = ArgumentParser( description = 'Plot PC loadings' )
        parser$add_argument(
                "--name",
                type = "character",
                help = "Name of region",
                required = TRUE
        )
        parser$add_argument(
                "--k",
                type = "integer",
                help = "Kmer size",
                required = TRUE
        )
        parser$add_argument(
                "--expected_haploid",
                type = "double",
                help = "Expected (haploid) kmer count",
                required = TRUE
        )
        parser$add_argument(
                "--regional",
                type = "character",
                help = "kmer counts from regional assembly",
                required = TRUE
        )
        parser$add_argument(
                "--whole_genome",
                type = "character",
                help = "kmer counts from whole genome assembly",
                required = TRUE
        )
        parser$add_argument(
                "--validation",
                type = "character",
                help = "kmer counts from validation data",
                required = TRUE
        )
	parser$add_argument(
		"--repeats",
		type = "character",
		help = "Repeatmasker file",
		required = TRUE
	)
        parser$add_argument(
                "--output",
                type = "character",
                help = "Name of output png file",
                required = TRUE
        )
        return( parser$parse_args())
}

args = parse_arguments()

load.repeatmasker <- function( filename ) {
	column_names = c(
		'sw_score', 'perc_div', 'perc_del', 'perc_ins',
		'query', 'query_begin', 'query_end', 'query_left', 'strand',
		'matching_repeat',
		'repeat_class', 'repeat_begin', 'repeat_end', 'repeat_left', 'ID', 'extra'
	)

	X = read.table(
		pipe(
			sprintf(
				"cat %s | sed -e 's/  */ /g' | sed -e 's/^ //'",
				filename
			)
		),
		header = F,
		col.names = column_names,
		skip = 3,
		fill = T,
		sep = " "
	)
	X$repeat_class = as.character( X$repeat_class )
	X$repeat_class2 = sapply( X$repeat_class, function(s) { strsplit(s, split = "/", fixed = T )[[1]][1] } )
	return(X)
}

data = list()
for( what in c( 'regional', 'validation', 'whole_genome' )) {
	data[[what]] = scan( pipe( sprintf( "cut -d' ' -f 2 %s", args[[what]] )), what = integer() )
}

stopifnot( length(data[[1]]) == length(data[[2]]) )
stopifnot( length(data[[1]]) == length(data[[3]]) )

repeatmasker = load.repeatmasker( args$repeats )

plot.data = data.frame(
	position = 1:length( data[[1]] ),
	k = args$k,
	regional = data[['regional']],
	whole_genome = data[['whole_genome']],
	validation = data[['validation']]
)	

plot.data$colour = 'black'
plot.data$colour[ plot.data$whole_genome > plot.data$regional ] = 'grey'

# Exclude highly repetetive kmers
# i.e. those with apparently more than 10 fold coverage in the validation data
w = which( plot.data$validation < args$expected_haploid * 20 )
plot.data = plot.data[w,]

png( file = args$output, width = 2400, height = 1000)

layout( matrix( 1:4, ncol = 1 ), heights = c( 0.1, 1, 1, 0.3 ))
par( cex = 2 )
par( mar = c( 0.1, 4, 0.1, 1 ) )

cat( "Plotting title...\n" )
plot( 0, 0, col = 'white', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = c( 0, 1 ), ylim = c( 0, 1 ) )
text( 0.5, 0.25, sprintf( "%s assembly vs. short-read validation data; k=%d", args$name, args$k ), font = 2, xpd = NA )

ylim = c( 0, max( plot.data$regional ))

cat( "Plotting raw values...\n" )

plot(
	plot.data$position,
	pmin( plot.data$validation / args$expected_haploid, ylim[2] ),
	pch = '.',
	col = 'grey',
	xlab = '', ylab = 'kmer multiplicity',
	ylim = ylim, xaxt = 'n', bty = 'n'
)
points( plot.data$position, pmin( plot.data$regional * 2, ylim[2] ), pch = '.', col = 'red' )
grid()
legend( "topright", pch = 19, col = c( "grey", "red" ), legend = c( "Validation kmers", "HV31 assembly (assuming homozygous)" ))

cat( "Plotting ratios...\n" )
par( mar = c( 0.1, 4, 0.1, 1 ))
plot(
	plot.data$position,
	pmin( plot.data$validation / plot.data$regional / args$expected_haploid / 2, 5 ),
	pch = '.',
	col = plot.data$colour,
	xlab = '', ylab = 'Ratio of multiplicity', ylim = c( 0, 5 ),
	bty = 'n',
	xaxt= 'n'
)
grid()
abline( h = 0:20, col = 'grey' )

legend(
	"topright",
	pch = 19,
	col = c( 'black', 'grey' ),
	legend = c( "Only in region", 'Also outside region' )
)

if( 0 ) {
	cat( "Plotting repeats...\n" )
	par( mar = c( 0.1, 4, 0.1, 1 ))

	repeat_classes = sort( unique( repeatmasker$repeat_class2 ) )
	print( repeat_classes )
	repeatmasker$repeat_class2 = factor( as.character( repeatmasker$repeat_class2 ), levels = repeat_classes )
	plot( 0, 0, col = 'white', xlab = '', ylab = '', bty = 'n', xlim = range( plot.data$position ), ylim = c( 1, length( repeat_classes) ), xaxt = 'n', yaxt = 'n' )
	rect(
		xleft = repeatmasker$query_begin,
		xright = repeatmasker$query_end,
		ybottom = as.integer( repeatmasker$repeat_class2 ) - 0.4,
		ytop = as.integer( repeatmasker$repeat_class2 ) + 0.4,
		col = 'grey'
	)
	text(
		min( plot.data$position ) - 10000,
		1:length( repeat_classes ),
		repeat_classes,
		adj = 1,
		cex = 0.5
	)
}

cat( "Plotting axis...\n" )
par( mar = c( 4, 4, 0.1, 1 ))
plot( 0, 0, col = 'white', xlab = 'Position in region', ylab = '', bty = 'n', xlim = range( plot.data$position ), ylim = c( 0, 1 ), yaxt = 'n' )

dev.off()

