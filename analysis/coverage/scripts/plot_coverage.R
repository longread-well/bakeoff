library( argparse )
library( wesanderson )
library( bakeoff )

ROOT = Sys.getenv( "BAKEOFF_ROOT" )

parser <- ArgumentParser( description = 'Plot coverage for one of the regions' )
parser$add_argument(
        '--region',
        type = "character",
        help = 'Region name to plot',
        required = TRUE
)
parser$add_argument(
        '--platforms',
        type = "character",
        nargs = '+',
	default = c( "nanopore", "pacbio-ccs", "pacbio-clr", "10x" ),
        help = 'Platforms to plot'
)
parser$add_argument(
        '--binsizes',
        type = "integer",
        nargs = '+',
	default = c( 100, 25, 10 ),
        help = 'Bin sizes to plot'
)
parser$add_argument(
	"--sample",
	type = "character",
	default = "JK_HV31",
	help = "sample to plot"
)
parser$add_argument(
	"--genes",
	type = "character",
	default = sprintf( "%s/projects/reference/GRCh38/UCSC/2019-08-27-GRCh38-refGene.tsv", ROOT ),
	help = "Path to genes"
)
parser$add_argument(
	"--range",
	type = "character",
	help = "Range to plot"
)
parser$add_argument(
        '--output',
        type = "character",
        nargs = 1,
        help = 'Where to put plot',
        required = TRUE
)

args = parser$parse_args()

replace <- function( what, fromto ) {
	result = what
	for( name in names( fromto )) {
		result = gsub( name, fromto[name], result, fixed = T )
	}
	return( result )
}

load.coverage.data = function( sample, region, platforms, binsizes ) {
	data = data.frame()
	for( platform in platforms ) {
		for( size in binsizes ) {
			filename = replace(
				"results/{sample}/{region}/{sample}_{region}_{platform}_coverage_{size}kb.bed",
				c(
					"{sample}" = sample,
					"{region}" = region,
					"{platform}" = platform,
					"{size}" = size
				)
			)
			cat( sprintf( "Loading coverage from \"%s\"...\n", filename ))
			X = read.delim( filename, hea=F, as.is=T, sep = "\t" )
			colnames(X) = c( "chromosome", "start", "end", "coverage" )
			data = rbind(
				data,
				cbind(
					data.frame(
						sample = sample,
						region = region,
						platform = platform,
						size = size
					),
					X
				)
			)
		}
	}
	return( data )
}

data = load.coverage.data( args$sample, args$region, args$platforms, args$binsizes )
data$coverage_per_base = data$coverage / ( data$size * 1000 )

cat( "Data is:\n" )
print( head( data ))

blank.plot <- function( xlim = c(0,1), ylim = c(0,1) ) {
	plot( 0, 0, col = 'white', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n', xlim = xlim, ylim = ylim )
}

plot.coverage <- function( data, region, args ) {
	N = length( args$platform )
	par( mar = c( 1, 2, 0.1, 1 ))
	region = c( region$start, region$end )
	for( platform in args$platform ) {
		w = which( data$platform == platform & data$size == 25 )
		cat( "Plotting:\n" )
		print( head(data[w,]))
		max.y = min( max( data$coverage_per_base[w], na.rm = T ), mean( data$coverage_per_base[w], na.rm = T ) * 3 )
		ylim = c( 0, max.y )
		blank.plot( xlim = region, ylim = ylim )
		rect(
			xleft = data$start[w],
			xright = data$end[w],
			ybottom = rep( 0, length(w)),
			ytop = pmin( data$coverage_per_base[w], max.y ),
			col = 'grey',
			border = NA
		)
		axis(2)
		grid()
		text( region[1], ylim[2] * 0.75, platform, cex = 2, font = 2, xpd = NA, adj = c( 0, 0 ) )	
		print( "HERE" )
	}
}

cat( "Loading genes...\n" )
genes = read.delim( sprintf( "%s/projects/reference/GRCh38/UCSC/2019-08-27-GRCh38-refGene.tsv", ROOT ), head = T, as.is=T, sep = "\t" )
genes$chromosome = genes$chrom
genes$name = genes$name2
genes$size = genes$txEnd - genes$txStart
genes = genes[ genes$cdsEnd > genes$cdsStart, ]
genes = genes[ order( genes$size, decreasing = T ), ]
genes = genes[ -which(duplicated( genes$name )), ]
print( head( genes ))

if( !is.null( args$range )) {
	plot.region = parse.ranges( args$range )
} else {
	plot.region = data.frame( chromosome = unique( data$chromosome ), start = min( data$start ), end = max( data$end ))
}

genes = genes[ which( genes$chromosome == plot.region$chromosome & genes$txEnd >= plot.region$start & genes$txStart <= plot.region$end ), ]

pdf( file = args$output, width = 6, height = 1 * length(args$platform ))
NP = length(args$platform)
NG = nrow( genes )
layout(
	matrix( 1:(NP+2), ncol = 1 ),
	heights = c( .5, rep( 1, NP ), 2 )
)
par( mar = c( 1, 2, 1, 1 ))
blank.plot()
text( 0.5, 0.5, sprintf( "%s: %s (%s:%d-%d)", args$sample, args$region, plot.region$chromosome, plot.region$start, plot.region$end ), font = 2, cex = 1, xpd = NA )
plot.coverage( data, plot.region, args )
par( mar = c( 4, 2, 0.1, 1 ))
plot.ucsc.genes( genes, region = plot.region, bty = 'n', height_in_inches = 1 )
grid()
dev.off()

