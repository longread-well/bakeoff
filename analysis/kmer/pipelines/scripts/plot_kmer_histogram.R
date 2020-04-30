args = commandArgs( trailingOnly = TRUE )
input = args[1]
output = args[2]
platform = args[3]
k = args[4]
bq = args[5]

X = read.table( input, hea=F, as.is=T );

height.at.1 = X[1,2]
dip = which.min( X[1:20,2] )
height.of.dip = X[dip,2]
peak = which.max( X[dip:nrow(X),2] ) + dip - 1
height.of.peak = X[peak,2]

if( peak > 30 ) {
	xlim = c( 0, 200 )
} else {
	xlim = c( 0, 100 )
}

pdf( file = output, width = 8, height = 4 );
par( mfrow = c( 1, 2 ))
plot( pmin( X[,1], xlim[2] ), X[,2], type = "l", xlab = "Coverage", ylab = "Count" ); 
markers = data.frame(
	x = c( dip, peak / 2, peak, round(1.5 * peak), 2 * peak ),
	label = c( "dip", "1", "2", "3", "4" )
)
markers$y = X[markers$x,2] + 0.05 * height.of.peak

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
legend(
	"topright",
	c(
		sprintf( "Platform: %s", platform ),
		sprintf( "k: %s; bq: %s", k, bq ),
		sprintf( "count at 1 (c1): %.2g", height.at.1 ),
		sprintf( "...at dip: %.2g", height.of.dip),
		sprintf( "...at peak (cp): %.2g", height.of.peak ),
		sprintf( "cp/c1: %.2g", height.of.peak / height.at.1 )
	),
	bty = "n"
);

m = max( X[ X[,1] > 10, 2 ] );
plot( pmin( X[,1], xlim[2]), X[,2], type = "l", xlab = "Coverage", ylab = "Count (zoomed)", ylim = c( 0, m * 1.25 ), xlim = c( 0, peak * 1.5  ), main = "...zoom") ;
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

