args = commandArgs( trailingOnly = TRUE )
input = args[1] 		# histogram file
output = args[2]		# output pdf file
sample = args[3]		# sample identifier
platform = args[4]		# Name of platform
k = args[5]			# kmer size
bq = args[6]			# bq filter

X = read.table( input, hea=F, as.is=T );

height.at.1 = X[1,2]

dip = min( which( X[-1,2] > X[-nrow(X),2] ))
height.of.dip = X[dip,2]
peak = which.max( X[dip:nrow(X),2] ) + dip - 1
height.of.peak = X[peak,2]

errors = sum( as.numeric(X[1:dip,2]) * as.numeric(X[1:dip,1]) )
repeats = sum( X[X[,1] > peak * 1.5,2] )
solid = sum( as.numeric( X[(dip+1):nrow(X),2] ) * as.numeric( X[(dip+1):nrow(X),1] ) )

xlim = c( 0, max( 50, peak * 2 ))

pdf( file = output, width = 6, height = 3 );
par( mfrow = c( 1, 2 ))

par( mar = c( 4.1, 1.1, 1.1, 0 ))
plot( 0, 0, col = 'white', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', bty = 'n' )
platformText = paste( strwrap( platform, width = 25 ), collapse = "\n" )

legend(
	"right",
	c(
		sprintf( "Sample: %s", sample ),
		sprintf( "Platform: %s", platformText ),
		sprintf( "k: %s; bq: %s", k, bq ),
		sprintf( "Errors (E): %.2g", errors ),
		sprintf( "Solid (S): %.2g", solid ),
		sprintf( "possible repeats (R): %.2g", repeats ),
		sprintf( "Repeat propn: %.2g%%", 100 * repeats / solid ),
		sprintf( "propn error kmers: %.2g%%", 100 * errors / ( errors + solid )),
		sprintf( "propn error bases %.2g%%", (100/as.integer(k)) * errors / ( errors + solid ))
	),
	bty = "n",
	cex = 1
);

par( mar = c( 4.1, 1.1, 1.1, 0.5 ))
plot( X[,1], X[,2], type = "l", xlab = "Coverage", ylab = "Count", xlim = xlim, ylim = c( 0, height.of.peak * 1.25 ) ); 
markers = data.frame(
	x = c( dip, peak / 2, peak, round(1.5 * peak), 2 * peak ),
	label = c( "dip", "1", "2", "3", "4" )
)
markers$y = X[markers$x,2] + 0.05 * height.of.peak
markers$transformed_y = X[markers$x,2]*X[markers$x,1] + 0.05 * height.of.peak

segments(
	x0 = markers$x, x1 = markers$x,
	y0 = rep(0, nrow(markers)), y1 = markers$y,
	lty = 2
)
text(
	markers$x[c(1,3,5)],
	markers$y[c(1,3,5)],
	sprintf( "%d", markers$x[c(1,3,5)] ),
	cex = 0.75,
	pos = 3
)

grid();

dev.off()

