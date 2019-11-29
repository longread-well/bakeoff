.onLoad <- function(libname, pkgname) {
	cat( "\n" )
	cat( "++ Welcome to the bakeoff R package.\n" )
	cat( "\n" )
	options <- options()
	bakeoff.opts <- list(
		bakeoff.root = Sys.getenv( "BAKEOFF_ROOT" )
	)
	if( bakeoff.opts$bakeoff.root == "" && !'bakeoff.root' %in% names( options ) ) {
		cat( "!! The BAKEOFF_ROOT environment variable is not set, and neither is `options( 'bakeoff.root' )`.\n" )
		cat( "!! Please use `options( bakeoff.root = <value> )` or set the BAKEOFF_ROOT environment variable, before package loading.\n" )
		cat( "\n" )
	}
	
	toset <- !(names( bakeoff.opts ) %in% names( options ))
	if( any( toset ) ) {
		options( bakeoff.opts[toset] )
	}
	invisible()
}

.onUnload <- function() {
	options( 'bakeoff.root' = NULL )
}
