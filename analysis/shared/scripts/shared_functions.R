load.sample.file <- function( filename, na.strings = "NA", sep = " " ) {
	# Read header and types
	samples = read.table( filename, as.is = T, hea=T, sep = sep )
	sample.header = colnames( samples )
	types = as.character( samples[1,] )
	# Now read the data rows
	samples = read.delim( filename, as.is = T, hea=F, skip = 2, sep = sep, na.strings = na.strings )
	names( samples ) = sample.header
	
	# In case IDs happened to be numerical, make the first two columns character vectors
	samples[,1] = as.character( samples[,1] )
	samples[,2] = as.character( samples[,2] )

	# Convert discrete levels to factors.
	for( i in 1:ncol( samples )) {
		if( types[i] == 'D' | types[i] == 'B' ) {
			samples[,i] = as.factor( samples[,i] )
		}
	}
	
	attr( samples, 'types' ) = types ;
	
	return( samples ) ;
}

save.sample.file <- function( samples, filename, sep = " " ) {
	stopifnot( length( attr( samples, "types" )) == ncol( samples ))
	write( colnames(samples), filename, ncol = 1000 )
	write( attr(samples, "types"), filename, ncol = 1000, append = T )
	write.table( samples, filename, col.names = F, row.names = F, quote = F, append = T )
}

fix_chromosome <- function( chromosome ) {
    chromosome = as.character( chromosome )
	chromosome = gsub( "chr", "", chromosome, fixed = T )
    n = nchar( chromosome )
    w = which( n == 1 ) ;
    chromosome[w] = sprintf( "0%s", chromosome[w] )
    return( chromosome )
}

load.sample.file <- function( filename, na.strings = "NA", sep = " " ) {
	# Read header and types
	samples = read.table( filename, as.is = T, hea=T, sep = sep )
	sample.header = colnames( samples )
	types = as.character( samples[1,] )
	# Now read the data rows
	samples = read.delim( filename, as.is = T, hea=F, skip = 2, sep = sep, na.strings = na.strings )
	names( samples ) = sample.header
	
	# In case IDs happened to be numerical, make the first two columns character vectors
	samples[,1] = as.character( samples[,1] )
	samples[,2] = as.character( samples[,2] )

	# Convert discrete levels to factors.
	for( i in 1:ncol( samples )) {
		if( types[i] == 'D' | types[i] == 'B' ) {
			samples[,i] = as.factor( samples[,i] )
		}
	}
	
	return( samples ) ;
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

load.genetic.map <- function( chromosomes = c( sprintf( "%d", 1:22 ), "X_PAR1", "X_PAR2", "X_nonPAR" ) ) {
    map <- data.frame()
    for( chromosome in chromosomes ) {
        padded_chromosome = chromosome ;
        if( chromosome %in% c( sprintf( "%d", 1:9 ) )) {
            padded_chromosome = paste( "0", chromosome, sep = "" )
        }
        path = "/kwiat/1/galton/malariagen/human/reference/mathgen.stats.ox.ac.uk/impute/2013-02-14_ALL_1000G_phase1integrated_v3"
        A <- cbind(
            chromosome = padded_chromosome,
            read.table(
                sprintf( "%s/genetic_map_chr%s_combined_b37.txt", path, chromosome ),
                hea = T, as.is = T
            )
        ) ;
        map = rbind( map, A )
    }
    map$chromosome = as.character( map$chromosome )
    map$chromosome[ which( map$chromosome == "X_nonPAR" ) ] = "0X"
    
    return(map) ;
}

load.genes <- function( filename, condense = TRUE ) {
    # Load genes
    gene <- read.delim(
	filename,
        header=TRUE,
        as.is=TRUE
    );
	gene <- gene[ -grep( "hap", gene$chrom ), ]
	
	if( condense ) {
	    gene <- gene[order(gene$txEnd - gene$txStart,decreasing=TRUE),];  #Get just longest transcript
	    gene <- gene[ !duplicated( gene$name2 ), ];
	}
    gene <- gene[ !is.na(gene$txStart), ];
	
    chromosome =  gsub( "^chr", "", gene$chrom )
    w1 = which( nchar( chromosome ) == 1 )
    chromosome[ w1 ] = sprintf( "0%s", chromosome[w1] )
    gene$chromosome = chromosome
    return( gene ) ;
}

#
# Given
# - a list of hapdb databases
# - a list of the samples we want from each db
# this function loads all the data (haplotypes, genotypes, or intensities) from each database
# returning data with missing values for variants / db combinations that are not present.
load.all <- function( ld.dbs, ld.samples = NULL, cohorts = NULL, what = "haplotypes", ... ) {
	if( is.null( cohorts )) {
		cohorts = names( ld.dbs )
	}
	if( is.null( ld.samples )) {
		ld.samples = list()
		for( name in cohorts ) {
			ld.samples[[ name ]] = ld.dbs[[ name ]]$samples$identifier
		}
	}

	if( what == 'haplotypes' ) {
		load.fn = load.haplotypes
	} else if( what == 'genotypes' ) {
		load.fn = load.genotypes
	} else if( what == 'intensities' ) {
		load.fn = load.intensities
	} else {
		stop( "You must specify what=\"haplotypes\", \"genotypes\", or \"intensities\"." )
	}

    haps = list( variant = data.frame(), data = matrix(), samples = data.frame() ) ;
    for( cohort in cohorts ) {
        cohort.data = load.fn( ld.dbs[[ cohort ]], samples = ld.samples[[ cohort ]], ... )
		cohort.data$variant$key = paste( cohort.data$variant$chromosome, cohort.data$variant$position, cohort.data$variant$alleleA, cohort.data$variant$alleleB, cohort.data$variant$rsid, sep = ":" )
		cohort.data$variant$key = paste( cohort.data$variant$chromosome, cohort.data$variant$position, cohort.data$variant$alleleA, cohort.data$variant$alleleB, cohort.data$variant$rsid, sep = ":" )
        #cohort.data$variant$key = paste( cohort.data$variant$chromosome, cohort.data$variant$position, cohort.data$variant$alleleA, cohort.data$variant$alleleB, sep = ":" )
        if( nrow( haps$variant ) == 0 ) {
            haps$variant = cohort.data$variant
            haps$data = cohort.data$data
			if( 'dosage' %in% names( cohort.data )) {
				haps$dosage = cohort.data$dosage
			}
        } else {
            w = which( !cohort.data$variant$key %in% haps$variant$key )
            if( length(w) > 0 ) {
                haps$variant = rbind( haps$variant, cohort.data$variant[w,] )
                haps$data = rbind( haps$data, matrix( NA, nrow = length(w), ncol = ncol( haps$data )) )
				if( 'dosage' %in% names( cohort.data )) {
					haps$dosage = rbind( haps$dosage, matrix( NA, nrow = length(w), ncol = ncol( haps$dosage )) )
				}
            }
            M = match( cohort.data$variant$key, haps$variant$key ) ;
            stopifnot( length( which( is.na( M ))) == 0 ) ;
            new.m = matrix( NA, nrow(haps$variant), ncol( cohort.data$data ), dimnames = list( haps$variant$rsid, colnames( cohort.data$data ) ) )
            new.m[M,] = cohort.data$data ;
            haps$data = cbind( haps$data, new.m )
			if( 'dosage' %in% names( cohort.data )) {
	            dosage.m = matrix( NA, nrow(haps$variant), ncol( cohort.data$dosage ), dimnames = list( haps$variant$rsid, colnames( cohort.data$dosage ) )) ;
	            dosage.m[M,] = cohort.data$dosage ;
				haps$dosage = cbind( haps$dosage, dosage.m )
			}
        }
        haps$samples = rbind( haps$samples, cohort.data$samples )
    }
    gc()
    return( haps )
}

get.available.annotations <- function( chromosome, region ) {
	ANALYSISDIR = Sys.getenv( "MG_ANALYSISDIR" )
	list(
		vep = list(
			name = "VEP",
			filename = sprintf( "%s/annotation/sqlite/meta_analysis.sqlite.symlink", ANALYSISDIR ),
			query = sprintf(
				paste( sep = ' ',
					"SELECT * FROM Variant V",
					"INNER JOIN VEP V2 ON (V2.variant_id == V.id) AND",
					"( Consequence IN",
					"( 'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification' )",
					"OR Consequence IN ( 'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant', 'regulatory_region_ablation' ))",
					"WHERE chromosome == '%s' AND position BETWEEN %d AND %d"
				),
				chromosome, region[1], region[2]
			),
			chromosome.column = 'chromosome',
			position.column = 'position',
			name.column = 'Consequence'
		),
		gwas2 = list(
			name = "NHGRI-EBI GWAS Catalog",
			filename = sprintf( "%s/annotation/sqlite/2015-06-30_NHGRI-EBI_GWAS_Catalog.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM Catalog WHERE CHR_ID == '%s' AND B37_CHR_POS BETWEEN %d AND %d AND `DISEASE_TRAIT` != 'Malaria'", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "CHR_ID",
			position.column = "B37_CHR_POS",
			name.column = "DISEASE_TRAIT"
		),
		eqtl = list(
			name = "GTEx v6p",
			filename = sprintf( "%s/annotation/sqlite/GTEx_Analysis_v6p.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM eqtl WHERE chromosome == '%s' AND position BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "position",
			name.column = "gene_id"
		),
		eqtl1 = list(
			name = "cis eqtl",
			filename = sprintf( "%s/annotation/sqlite/Westra_et_al_eqtls.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT SNPChr AS chromosome, snp142_position AS position, SNPName AS rsid, HUGO, pvalue, FDR FROM cis WHERE SNPChr == '%s' AND snp142_position BETWEEN %d AND %d AND FDR < 0.01", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "position",
			name.column = "HUGO"
		),
		eqtl2 = list(
			name = "trans eqtl",
			filename = sprintf( "%s/annotation/sqlite/Westra_et_al_eqtls.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT SNPChr AS chromosome, snp142_position AS position, SNPName AS rsid, HGNCName, pvalue FROM trans WHERE chromosome == '%s' AND snp142_position BETWEEN %d AND %d AND FDR < 0.01", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "position",
			name.column = "HGNCName"
		),
		ABS1 = list(
			name = "Ancient balancing selection loci (*)",
			filename = "/kwiat/1/galton/malariagen/human/reference/articles/Leffler_et_al_Multiple_Instances_of_Ancient_Balancing_Selection_Shared_Between_Humans_and_Chimpanzees/TableS4_SharedHaplotypeRegions.csv",
			chromosome.column = "Chr (hg19)",
			position.column = "Position (hg19)",
			name.column = "rs#",
			reference = "Leffler et al Science (2013)"
		),
		ABS2 = list(
			name = "Ancient balancing selection loci (*)",
			filename = "/kwiat/1/galton/malariagen/human/reference/articles/Leffler_et_al_Multiple_Instances_of_Ancient_Balancing_Selection_Shared_Between_Humans_and_Chimpanzees/TableS5_SharedCodingSNPs.csv",
			chromosome.column = "Chr (hg19)",
			position.column = "Position (hg19)",
			name.column = "rs#",
			reference = "Leffler et al Science (2013)"
		),
		rnaseq.CD34 = list(
			name = "RNA expression in CD34+ cells (‡)",
			filename = sprintf( "%s/annotation/sqlite/erythrocyte_precursor_rnaseq_coverage.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM CD34 WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			height.column = "coverage",
			transform = log10,
			reference = "Li et al, Blood 124)"
		),
		rnaseq.erythroblast = list(
			name = "RNA expression in erythroblasts (†)",
			filename = sprintf( "%s/annotation/sqlite/erythrocyte_precursor_rnaseq_coverage.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM erythroblast WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			height.column = "coverage",
			transform = log10
		),
		chromhmm.proerythroblast = list(
			name = "Chromatin state (proerythroblast)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM adult_15 A INNER JOIN states S ON S.state == A.state WHERE chromosome == '%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.liver = list(
			name = "Chromatin state (liver)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E066_liver A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.spleen = list(
			name = "Chromatin state (spleen)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E113_spleen A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.T = list(
			name = "Chromatin state (primary T)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E034_primary_T A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.Thelper = list(
			name = "Chromatin state (primary T helper memory)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E037_CD4_Memory_Primary_Cells A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.B = list(
			name = "Chromatin state (primary B)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E032_primary_B A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.monocytes = list(
			name = "Chromatin state (primary monocytes)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E029_primary_monocytes A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.cd34 = list(
			name = "Chromatin state (hsc)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E035_cd34_stem_cells A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chromhmm.K562 = list(
			name = "Chromatin state (K562)",
			filename = sprintf( "%s/annotation/sqlite/roadmap_chromhmm_segments.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT chromosome, start, end, A.state AS state, color FROM E123_K562_Leukemia A INNER JOIN states S ON S.state == A.state WHERE chromosome == 'chr%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = "color"
		),
		chipseq.GATA1.erythroblast = list(
			name = "GATA1; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM970258_GATA1-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0
		),
		chipseq.NFE2.erythroblast = list(
			name = "NFE2; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM908059_NFE2-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0
		),
		chipseq.TAL1.erythroblast = list(
			name = "TAL1; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM908055_TAL1-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0,
			reference = "Xu et al, Dev. Cell 2012"
		),
		chipseq.CTCF.erythroblast = list(
			name = "CTCF; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM908065_CTCF-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0,
			reference = "Xu et al, Dev. Cell 2012"
		),
		chipseq.RAD21.erythroblast = list(
			name = "RAD21; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM908067_RAD21-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0,
			reference = "Xu et al, Dev. Cell 2012"
		),
		chipseq.PolII.erythroblast = list(
			name = "PolII; adult ProEs (†)",
			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM `GSM908069_PolII-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "END",
			height.column = "score",
			transform = log10,
			min = 0,
			reference = "Xu et al, Dev. Cell 2012"
		),
#
#		chipseq.H3K4me1.erythroblast = list(
#			name = "H3K4me1 (enhancer,tss downstream); ProEs (†)",
#			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
#			query = sprintf( "SELECT * FROM `GSM908035_H3K4me1-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
#			chromosome.column = "chromosome",
#			position.column = "start",
#			end.position.column = "END",
#			height.column = "score",
#			transform = log10,
#			min = 0,
#			reference = "Xu et al, Dev. Cell 2012"
#		),
#		chipseq.H3K4me2.erythroblast = list(
#			name = "H3K4me2 (promoter,enhancer); ProEs (†)",
#			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
#			query = sprintf( "SELECT * FROM `GSM908037_H3K4me2-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
#			chromosome.column = "chromosome",
#			position.column = "start",
#			end.position.column = "END",
#			height.column = "score",
#			transform = log10,
#			min = 0,
#			reference = "Xu et al, Dev. Cell 2012"
#		),
#		chipseq.H3K4me3.erythroblast = list(
#			name = "H3K4me3 (promoter,tss); ProEs (†)",
#			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
#			query = sprintf( "SELECT * FROM `GSM908039_H3K4me3-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
#			chromosome.column = "chromosome",
#			position.column = "start",
#			end.position.column = "END",
#			height.column = "score",
##			transform = log10,
#			min = 0,
#			reference = "Xu et al, Dev. Cell 2012"
#		),
#		chipseq.H3K27me3.erythroblast = list(
#			name = "H3K27me3 (repression); ProEs (†)",
#			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
#			query = sprintf( "SELECT * FROM `GSM908043_H3K27me3-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
#			chromosome.column = "chromosome",
#			position.column = "start",
#			end.position.column = "END",
#			height.column = "score",
#			transform = log10,
#			min = 0,
#			reference = "Xu et al, Dev. Cell 2012"
#		),
#		chipseq.H3K36me3.erythroblast = list(
#			name = "H3K36me3 (elongation); ProEs (†)",
#			filename = sprintf( "%s/annotation/sqlite/proerythroblast_chipseq.sqlite.symlink", ANALYSISDIR ),
#			query = sprintf( "SELECT * FROM `GSM908047_H3K36me3-A_peaks` WHERE chromosome == 'chr%s' AND start BETWEEN %d AND %d", gsub( '^0', '', chromosome ), region[1], region[2] ),
#			chromosome.column = "chromosome",
#			position.column = "start",
#			end.position.column = "END",
#			height.column = "score",
#			transform = log10,
#			min = 0,
#			reference = "Xu et al, Dev. Cell 2012"
#		),
		astle.rbc = list(
			name = "Haemotopetic traits",
			filename = sprintf( "%s/annotation/sqlite/astle_et_al_haematopoetic_trait_summary_statistics.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT * FROM Association A INNER JOIN Variant V ON V.id == A.variant_id WHERE A.chromosome == '%s' AND A.position BETWEEN %d AND %d AND P < 8.31E-9", gsub( '^0', '', chromosome ), region[1], region[2] ),
			chromosome.column = 'chromosome',
			position.column = 'position',
			name.column = 'trait'
		),
		TADs = list(
			name = "TADs",
			filename = sprintf( "%s/annotation/sqlite/Javierre_et_al_TADs.sqlite.symlink", ANALYSISDIR ),
			query = sprintf( "SELECT *, '128,128,128' AS colour FROM `TADs` WHERE chromosome == '%s' AND start < %d AND end > %d", gsub( '^0', '', chromosome ), region[2], region[1] ),
			name.column = 'cell_type',
			chromosome.column = "chromosome",
			position.column = "start",
			end.position.column = "end",
			colour.column = 'colour',
			reference = "Javierre et al"
		)
	)
}

load.annotation.data <- function(
    annotations,
	range = NULL,
	verbose = FALSE
) {
    annotation.data = data.frame()
    for( name in names( annotations )) {
        annotation = annotations[[ name ]]
		cat( sprintf( "Loading annotations for \"%s\"...\n", name ) )
		if( 'query' %in% names( annotation )) {
			if( verbose ) {
				cat( sprintf( "Running query: \"%s\"...\n", annotation$query ) )
			}
			db = dbConnect( dbDriver( "SQLite" ), annotation$filename, flags = SQLITE_RO )
			X = dbGetQuery( db, annotation$query )
			dbDisconnect( db )
			rm(db)
		} else {
	        if( substring( annotation$filename, nchar( annotation$filename ) - 3, nchar( annotation$filename ) ) == '.csv' ) {
	            X = read.csv( annotation$filename, as.is = TRUE, header = T, check.names = F ) ;
	        } else {
	            X = read.delim( annotation$filename, as.is=TRUE, header=TRUE, check.names = F ) ;
	        }
		}
        if( "position.offset" %in% names( annotation ) ) {
            X[, annotation$position.column ] = X[, annotation$position.column ] + annotation$position.offset ;
        }

		if( verbose ) {
			print( head(X) )
		}
		score = rep( NA, nrow(X) )
		if( "score.column" %in% names( annotation )) {
			if( "min.score" %in% names( annotation ) ) {
				X = X[ X[, annotation$score.column ] > annotation$min.score, , drop = F ]
			}
			score = X[, annotation$score.column ] / max( X[, annotation$score.column ], na.rm = T )
		}
		if( verbose ) {
			print( head(X) )
		}
		if( !'name.column' %in% names( annotation )) {
			annotation$name.column = 'name'
			if( !'name' %in% colnames(X)) {
				X[, 'name'] = annotation$name
			}
		}
		if( nrow(X) == 0 ) {
			cat( "ZERO ROWS.\n" )
		}


		if( nrow(X) > 0 ) {
		if( "end.position.column" %in% names( annotation )) {
	        if( "position.offset" %in% names( annotation ) ) {
	            X[, annotation$end.position.column ] = X[, annotation$end.position.column ] + annotation$position.offset ;
	        }
			if( "height.column" %in% names( annotation )) {
				if( "transform" %in% names( annotation )) {
					X[, annotation$height.column ] = annotation$transform( X[, annotation$height.column ] )
				}
				if( "min" %in% names( annotation )) {
					X = X[ X[, annotation$height.column ] >= annotation$min, ]
				}
				if( "max" %in% names( annotation )) {
					X = X[ X[, annotation$height.column ] <= annotation$max, ]
				}
				if( "upper.threshhold" %in% names( annotation )) {
					X[, annotation$height.column ] = pmin( X[, annotation$height.column], annotation$upper.threshhold )
				}
				max.height = max( X[, annotation$height.column ], na.rm = T )
				Z = data.frame(
					type = "continuous",
					label = annotation$name,
					name = X[, annotation$name.column ],
					chromosome = fix_chromosome( X[, annotation$chromosome.column ] ),
					position = X[, annotation$position.column ],
					end_position = X[, annotation$end.position.column ],
					height = X[, annotation$height.column ] / max.height,
					colour = NA,
					score = score 
				)
				if( verbose ) {
					print( head( Z ))
				}
				annotation.data = rbind( annotation.data, Z )
			} else if( 'colour.column' %in% names(annotation)) {
				colours = sapply(
					X[,annotation$colour.column],
					function(s) {
						bits = as.integer( strsplit( s, split = ",", fixed = T )[[1]] )
						rgb( bits[1], bits[2], bits[3], maxColorValue = 255 )
					}
				)
				annotation.data = rbind(
		            annotation.data,
		            data.frame(
		                type = "segment",
		                label = annotation$name,
		                name = X[, annotation$name.column ],
		                chromosome = fix_chromosome( X[, annotation$chromosome.column ] ),
		                position = X[, annotation$position.column ],
						end_position = X[, annotation$end.position.column ],
						height = NA,
						colour = colours,
						score = score 
					)
				)
			} else {
				annotation.data = rbind(
		            annotation.data,
		            data.frame(
		                type = "segment",
		                label = annotation$name,
		                name = X[, annotation$name.column ],
		                chromosome = fix_chromosome( X[, annotation$chromosome.column ] ),
		                position = X[, annotation$position.column ],
						end_position = X[, annotation$end.position.column ],
						height = NA,
						colour = rgb( 0.6, 0.1, 1, score ),
						score = score 
		            )
		        )
			}
		} else {
	        annotation.data = rbind(
	            annotation.data,
	            data.frame(
	                type = "point",
	                label = annotation$name,
	                name = X[, annotation$name.column ],
	                chromosome = fix_chromosome( X[, annotation$chromosome.column ] ),
	                position = X[, annotation$position.column ],
					end_position = NA,
					height = NA,
					colour = NA,
					score = score 
	            )
	        )
		}
		}
		if( verbose ) {
			cat( sprintf( "Loaded annotations for \"%s\", they are now:\n", name ) ) ;
			print( tail( annotation.data ) )
		}
    }
	if( !is.null( range )) {
		annotation.data = annotation.data[
			which(
				annotation.data$chromosome == as.character( range$chromosome )
				& annotation.data$position >= range$start
				& annotation.data$position <= range$end
			),
			, drop = FALSE
		]
	}
	cat( "Annotation data loaded.\n" )
    return( annotation.data ) ;
}

test.d_prime <- function() {
    # All possible quadruples of haplotypes.
    # Frequency is either 0, 0.25, 0.5, 0.75, or 1
    # Up to permutation the only pair of SNPs where |D'| could be less than 1 is
    # 1 1 0 0
    # 1 0 1 0
    # and these two rows are orthogonal, so D' = 0 there.
    # Thus every entry should be NA, 0, 1 or -1.
    test1 = as.matrix( expand.grid( 0:1, 0:1, 0:1, 0:1 ) )
    compute.d_prime_matrix( test1 )$D
    compute.d_prime_matrix( test1 )$Dprime
    compute.d_prime_matrix( test1 )$frequency

    # test2: all possible quintuples of haplotypes, now we can start to see nonzero |D'| != 1.
    # Up to permutation or flip of alleles the only haplotype combinations with D' != 0 or 1 shoul dbe
    # 1 1 1 0 0
    # 1 0 0 1 0
    # with D = 0.2 - ( 0.6 * 0.4 ) = -0.04
    # and D' = D / ( 6/25) = -0.16666

    # or
    # 1 1 1 0 0
    # 1 1 0 1 0
    # with D = 0.4 - ( 0.6 * 0.6 ) = 0.04
    # and D' = D / ( 6/25) = 0.166666
    # So the only values in our matrix should be +/- 1 and +/1 0.16666
    test2 = as.matrix( expand.grid( 0:1, 0:1, 0:1, 0:1, 0:1 ) )
    compute.d_prime_matrix( test2 )$D
    compute.d_prime_matrix( test2 )$Dprime
    compute.d_prime_matrix( test2 )$frequency
}

plot.region <- function(
	name, local, hit, tags, region, ld.breaks, ld.colours, Dprime.breaks, Dprime.colours,
	pt.cex = 1, annotate.tags = TRUE, annotate.rsids = c(), rsid.cex = 1, main = "",
	what = "mean_bf",
	scale.fn = log10,
	ylab = "log10 model-averaged BF",
	ymax = NULL,
	legends = TRUE,
	axis = "top",
	...
) {
    #Colour for LD
    col = get.ld.colours( local$r_squared, ld.breaks, ld.colours )
    Dprime.col = get.ld.colours( local$absolute_Dprime, Dprime.breaks, Dprime.colours )

	if( what == 'pvalue' ) {
		ymin = 1
	} else {
		ymin = 0
	}

	print( ymax )
	print( head( local[,what] ))
	if( is.null( ymax )) {
		if( what == "maller_et_al_posterior" ) {
			ymax = 1.05 * max( scale.fn( local[, what] ) );
		} else {
			ymax = 1.05 * max(5,max( scale.fn( local[, what] ), na.rm = T ));
		}
	} else {
		raw.scale.fn = scale.fn
		scale.fn = function(x) {
			pmin( ymax, raw.scale.fn( x ))
		}
		ymax = max( ymax, 1.05 * max(5,max( scale.fn( local[, what] ), na.rm = T )));
	}
	print( what )
	print( ymax )
    plot(
        local$position, scale.fn( local[, what] ),
        xlim = region, ylim=c(ymin,ymax),
        ylab = ylab, xlab="",
        main = main,
        pch=19, col=col,
        xaxt = 'n',
		cex = pt.cex,
		...
    ) ;
    ticks = axTicks(3)
	if( axis == "top" ) {
		axis( side = 3, at = ticks, labels = sprintf( "%.2fMb", as.numeric( ticks ) / 1000000 ), tcl = 0.5, mgp =c(3,0,0))
	} else if( axis == "bottom" ) {
		axis( side = 3, at = ticks, labels = sprintf( "%.2fMb", as.numeric( ticks ) / 1000000 ), tcl = 0.5, mgp =c(3,0,0))
	}

    typed = local$typed
    wDprime = which( !is.na( Dprime.col ) & Dprime.col != Dprime.colours[1] )
    if( length( wDprime ) > 0 ) {
        points(
            local$position[wDprime], scale.fn( local[, what][wDprime] ), pch=19, col = Dprime.col[ wDprime ], cex = pt.cex / 2.5
        ) ;
    }

    points( local$position[typed], scale.fn( local[, what][typed] ),pch=3);
    text( hit$position, scale.fn( hit[, what] ),hit$rsid, pos=3, cex = rsid.cex );
	if( length( annotate.rsids ) > 0 ) {
		annotate.rsids.M = match( annotate.rsids, local$rsid )
		text( local$position[ annotate.rsids.M ], scale.fn( local[, what][ annotate.rsids.M ] ), annotate.rsids, pos = 3, cex = rsid.cex )
	}

    if( nrow( tags ) > 0 ) {
		if( annotate.tags ) {
			text( tags$position, scale.fn( tags[, what] ), tags$rsid, pos=1, col="darkgrey", cex = rsid.cex ) ;
		}
		if( !'cex' %in% colnames( tags )) {
			tags$cex = pt.cex
		}
        points( tags$position, scale.fn( tags[, what] ), pch = tags$pch, cex = tags$cex ) ;
    }

	if( legends ) {
		legend(
			x = 'topright',
			legend = c( 'Imputed variant', 'Omni 2.5M variant', 'Sequenom-typed variant', 'Imputed SV' ),
			pch = c( 19, 3, 2, 1 ),
			cex = 1.1,
			bty = 'n'
		)
		legend(
			x = 'topleft',
			legend = c(
				sprintf( "r² <= %.1f", ld.breaks[2] ),  sprintf( "r² > %.1f", ld.breaks[ -c( 1, length( ld.breaks ) ) ] ),
				sprintf( "|D'| <= %.1f", Dprime.breaks[2] ),  sprintf( "|D'| > %.1f", Dprime.breaks[ -c( 1, length( Dprime.breaks ) ) ] )
			),
			col =c( ld.colours[-length(ld.colours)], Dprime.colours[-length(Dprime.colours)] ),
			pch = 19,
			pt.cex = c(
				rep( pt.cex, length( ld.breaks ) -1 ),
				rep( pt.cex / 2.5, length( Dprime.breaks ) - 1 )
			),
			cex = 1.1,
			ncol = 2,
			bty = 'n'
		)
	}
}

get.ld.colours <- function( values, breaks = seq( 0, 1, length=6 ), colours = c( "lightgrey", rev( heat.colors(6))[-1] ) ) {
    stopifnot( length( breaks ) == length( colours ) )
    result = rep( colours[1], length( values ))
    for(b in 2:length( breaks )) {
        result[ (values > breaks[b] & values <= breaks[b+1]) ] <- colours[b]
    }
    return( result )
}

plot.annotations <- function(
	region,
	annotation.data,
	ld.breaks = NULL,
	ld.colours = NULL,
	Dprime.breaks = NULL,
	Dprime.colours = NULL,
	pt.cex = 1,
	levels = NULL,
	focus.position = NULL,
	label.position = 'above',
	...
) {
	if( !is.null( ld.breaks )) {
		col = get.ld.colours( annotation.data$r_squared, ld.breaks, ld.colours )
	} else {
		col = rep( "black", nrow( annotation.data ))
	}

	if( !is.null( Dprime.breaks )) {
		Dprime.col = get.ld.colours( annotation.data$absolute_Dprime, Dprime.breaks, Dprime.colours )
	} else {
		Dprime.col = rep( "black", nrow( annotation.data ))
	}

	if( is.null( levels )) {
		annotation.data$label = as.factor( annotation.data$label ) ;
		levels = levels( annotation.data$label )
	} else {
		annotation.data$label = factor( as.character( annotation.data$label ), levels = levels )
	}
	annotation.data$level = as.integer( annotation.data$label )
#	annotation.data = annotation.data[order( annotation.data$level, annotation.data$r_squared, decreasing = F ),]
    print( table( annotation.data$label, annotation.data$level ) )
    plot( 0, 0, xlim = region, ylim = c( 0.5, length(levels) + 1 ), xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xpd = NA, xaxs = 'i' )
    segments( x0 = region[1], x1 = region[2], y0 = 1:length(levels), y1 = 1:length( levels ), col = rgb( 0, 0, 0, 0.2 ))
	if( !is.null( focus.position )) {
		segments( x0 = focus.position, x1 = focus.position, y0 = 0.5, y1 = length(levels) + 0.5, col = rgb( 0, 0, 0, 0.2 ))
	}
    if( nrow( annotation.data ) > 0 ) {
	wPoint = which( annotation.data$type == 'point' )# & ( col != ld.colours[1] | Dprime.col != Dprime.colours[1] ) )
	if( length( wPoint ) > 0 ) {
		r2_cex = rep( pt.cex, length( wPoint ))
		r2_cex[ which( col[wPoint]== ld.colours[1]  & Dprime.col[wPoint] == Dprime.colours[1] ) ] = (pt.cex / 3)
	    points(
			annotation.data$position[ wPoint ], annotation.data$level[ wPoint ],
			pch = 16, col = col[ wPoint ], cex = r2_cex
		) ;
	}
	wDprime = which( annotation.data$type == 'point' & !is.na( Dprime.col ) & Dprime.col != Dprime.colours[1] )
	if( length( wDprime ) > 0 ) {
		points(
			annotation.data$position[wDprime], annotation.data$level[wDprime],
			pch=19, col = Dprime.col[ wDprime ], cex = pt.cex/3
		) ;
	}
	wLabel = which( annotation.data$type == 'point' & annotation.data$r_squared > 0.1 )
	if( length( wLabel > 0 )) {
		for( level in unique( annotation.data$level[ wLabel ] )) {	
			wLevel = intersect( wLabel, which( annotation.data$level == level ) )
			if( length( wLevel ) < 10 ) {
				if( level == 1 ) {
					text.pos = 1
				} else {
					text.pos = 3
				}
				text(
					annotation.data$position[ wLevel ], annotation.data$level[ wLevel ],
					annotation.data$name[ wLevel ],
					pos = text.pos, cex = 0.75,
					xpd = NA
				)
			}
		}
	}
	wSegment = which( annotation.data$type == 'segment' )
	if( length( wSegment ) > 0 ) {
	    rect(
			xleft = annotation.data$position[wSegment], xright = annotation.data$end_position[wSegment],
			ybottom = annotation.data$level[wSegment] - 0.15, ytop = annotation.data$level[wSegment] + 0.15,
			col = annotation.data$colour[wSegment], cex = pt.cex,
			border = NA
		) ;
	}
	wBlock = which( annotation.data$type == 'block' )
	if( length( wBlock ) > 0 ) {
	    rect(
			xleft = annotation.data$position[wBlock], xright = annotation.data$end_position[wBlock],
			ybottom = annotation.data$level[wBlock] - 0.25, ytop = annotation.data$level[wBlock] + 0.25,
			col = annotation.data$col,
			border = NA
		) ;
	}
	wBar = which( annotation.data$type == 'continuous' )
	if( length( wBar ) > 0 ) {
	    rect(
			xleft = annotation.data$position[wBar], xright = annotation.data$end_position[wBar],
			ybottom = annotation.data$level[wBar], ytop = annotation.data$level[wBar] + 0.5 * annotation.data$height[wBar],
			col = rgb( 0.5, 0.5, 0.5, 1 ),
			border = rgb( 0.5, 0.5, 0.5, 1 )
		) ;
	}
	}
	if( label.position == 'above' ) {
		text( region[1], (1:length( levels ))+0.5, label = levels, adj = c( 0, 0.5 ), cex = 1 )
	} else if( label.position == 'left' ) {
		text( region[1] - (region[2] - region[1]) / 50, (1:length( levels )), label = levels, adj = c( 1, 0.5 ), cex = 1, xpd = NA )
	} else if( label.position == 'right' ) {
		text( region[2] + (region[2] - region[1]) / 50, (1:length( levels )), label = levels, adj = c( 0, 0.5 ), cex = 1, xpd = NA )
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
) {
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

make.forest.plot <- function( model.fit, left = -1, right = 1, point.scale = 2, margin = c( 0.5, 0.5 ), populations = c( "Gambia", "Malawi", "Kenya" ) ) {
    if( model.fit$chromosome == '0X' ) {
        models = expand.grid( c( "males", "females" ), populations, c( "add", "rec", "dom", "het" ) )
        models = sprintf( "%s:%s:%s", models[,2], models[,1], models[,3] )
        betas = as.numeric( model.fit[, sprintf( "%s:beta_1", models ) ] )
        ses = as.numeric( model.fit[, sprintf( "%s:se_1", models ) ] )
        names( betas ) = names(ses) = models ;
        stopifnot( length( betas ) == length( ses ))
        stopifnot( length( betas ) == 4*3*2 )
    } else {
        models = expand.grid( populations, c( "add", "rec", "dom", "het" ) )
        models = sprintf( "%s:%s", models[,1], models[,2] )
        betas = as.numeric( model.fit[, sprintf( "%s:beta_1", models ) ] )
        ses = as.numeric( model.fit[, sprintf( "%s:se_1", models ) ] )
        names( betas ) = names(ses) = models ;
        stopifnot( length( betas ) == length( ses ))
        stopifnot( length( betas ) == 4*3 )
    }

    frequency.names = rep( sprintf( "%s:B_allele_frequency", populations ), 4 )
    frequency = data.frame(
        name = frequency.names,
        value = as.numeric( model.fit[, frequency.names ] )
    )

    data = data.frame(
        name = models,
        beta = betas,
        beta_lower = betas - 1.96 * ses,
        beta_upper = betas + 1.96 * ses
    )
    data$beta_clipped = pmax( pmin( data$beta, right ), left )
    data$beta_lower_clipped = pmax( data$beta_lower, left )
    data$beta_upper_clipped = pmin( data$beta_upper, right )

    frequency$colour = NA
    frequency$colour[ grep( "Gambia", frequency$name )] = rgb( 1, 0.2, 0.2 )
    frequency$colour[ grep( "Malawi", frequency$name )] = rgb( 0.5, 1, 0.7 )
    frequency$colour[ grep( "Kenya", frequency$name )] = rgb( 0.3, 0.3, 1 )
    
    data$colour = NA
    data$colour[ grep( "Gambia", data$name )] = rgb( 1, 0.2, 0.2 )
    data$colour[ grep( "Malawi", data$name )] = rgb( 0.5, 1, 0.7 )
    data$colour[ grep( "Kenya", data$name )] = rgb( 0.3, 0.3, 1 )
    
    L = 18 ;
    vpos = c( (L-0):(L-2), (L-4):(L-6), (L-8):(L-10), (L-12):(L-14) )
    
    if( ( model.fit$chromosome == '0X' ) ) {
        beta.vpos = vpos[ sort( rep( 1:length( vpos ), 2 ))]
        beta.vpos[ seq( from = 1, by = 2, length = length( vpos )) ] = beta.vpos[ seq( from = 1, by = 2, length = length( vpos )) ] + (0.2)
        beta.vpos[ seq( from = 2, by = 2, length = length( vpos )) ] = beta.vpos[ seq( from = 2, by = 2, length = length( vpos )) ] - (0.2)
    } else {
        beta.vpos = vpos
    }
    
    divider = c( L-3, L-7, L-11, L-15 )

    {
        plot(
            x = data$beta, y = beta.vpos,
            xlim = c( left - margin[1], right + margin[2] ),
            ylim = c( 0, L ),
            col = data$colour,
            pch = 16,
            cex = 0.8,
            yaxt = 'n',
            ylab = '',
            xaxt = 'n',
            xlab = "OR",
            bty = 'n',
            main = sprintf( "%s %s/%s", model.fit$rsid, model.fit$alleleA, model.fit$alleleB )
        )
        
        points(
            x = rep( right + 0.1, nrow( frequency ) ),
            y = vpos,
            cex = point.scale * sqrt( frequency$value ),
            col = frequency$colour,
            pch = 19
        )
        
        text(
            x = rep( right + 0.17, nrow( frequency ) ),
            y = vpos,
            sprintf( "%.1f%%", frequency$value * 100 ),
            adj = 0
        )
        
    #    text( 0, L+1, sprintf( "%s %s/%s", model.fit$rsid, model.fit$alleleA, model.fit$alleleB ), font = 2 )
        text( x = left - margin[1], y = vpos, label = gsub( ":B_allele_frequency", "", frequency.names ), adj = 0 )
        segments( x0 = left, x1 = right, y0 = beta.vpos, y1 = beta.vpos, col = rgb( 0, 0, 0, 0.2 ) )

        axis( side = 1, at = c( log( 1/2 ), log( 1 ), log(2) ), labels = c( "½", "1", "2" ) )
        segments( x0 = 0, x1 = 0, y0 = min( vpos ) - 0.5, y1 = max(vpos)+0.5, lty = 2, col = "red" )
        segments( x0 = data$beta_lower_clipped, x1 = data$beta_upper_clipped, y0 = beta.vpos, y1 = beta.vpos )
        w = which( data$beta_lower < data$beta_lower_clipped )
        if( length(w) > 0 ) {
            segments( x0 = left-0.2, y0 = beta.vpos[w] , x1 = left, y1 = beta.vpos[w], lty = 2 )
        }
        w = which( data$beta_upper > data$beta_upper_clipped )
        if( length(w) > 0 ) {
            segments( x0 = right, y0 = beta.vpos[w] , x1 = right + 0.2, y1 = beta.vpos[w], lty = 2 )
        }
        text( left, divider, c( "add", "rec", "dom", "het" ), cex = 0.8 )
    
        points(
            seq( from = left + 0.2, to = right - 0.2, length = 5 ),
            rep( 1, 5 ),
            pch = 16,
            cex = sqrt( seq( from = 0.1, to = 0.9, by = 0.2 )) * point.scale
        )

        text(
            seq( from = left + 0.2, to = right - 0.2, length = 5 ),
            rep( 0.5, 5 ),
            label = sprintf( "%.1f", seq( from = 0.1, to = 0.9, by = 0.2 ) ),
            adj = c( 0.5, 1 )
        )
    
        text( left + 0.1, 0.75, "Frequency:", adj = 1 )
    }
}

annotate.list.with.nearest.gene <- function( top, genes ) {
    for( name in c(
        'nearest_gene',
        'nearest_gene_distance'
        )
    ) {
        top[, name ] = NA
    }
    for( i in 1:nrow( top )) {
        chromosome = as.character( top$chromosome[i] )
        region = c( top$min[i], top$max[i] )
        proximal.genes <- genes[ genes$chromosome == chromosome, ]
        
        # Compute distance to each gene
        proximal.genes$distance = pmax( pmax( proximal.genes$txStart - top$position[i], top$position[i] - proximal.genes$txEnd ), 0 )
        proximal.genes <- proximal.genes[ order( proximal.genes$distance ), ]
        wMin = which.min( proximal.genes$distance )
        if( length( wMin ) > 0 ) {
            top[i,"nearest_gene"] = sprintf( "%s", proximal.genes$name2[wMin] )
            top[i,"nearest_gene_distance"] = sprintf( "%s(%.0fkb)", proximal.genes$name2[wMin], proximal.genes$distance[wMin]/1000 )
        }
	}
	return( top )
}

qqplot <- function( ps, plot.new = TRUE, limit = 50, CI = 0.95, thin = TRUE, scale = "chisq", df = 1, xlab = NULL, ylab = NULL, labels = NULL, ... ) {
    cat( "Sorting p-values...\n" )
	wNotNA = which( !is.na( ps ))
	ps = ps[wNotNA]
	p.value.order = order( ps )
    ps = ps[p.value.order]
	if( !is.null( labels )) {
		labels = labels[wNotNA,]
		labels = labels[p.value.order,]
	}
	if( scale == "chisq" ) {
	    lambda = qchisq( median( ps ), df = df, lower.tail = F )
	    lambda = lambda / qchisq( 0.5, df = df, lower.tail = F )

	    qs = qchisq( ps, df = df, lower.tail = F )
	    expected = (1:length(ps))/(length(ps)+1)
	    expected = qchisq( expected, df = df, lower.tail = F )
		
		what.to.plot = "χ² statistic"
		
	} else if( scale == "log10" ) {
	    lambda = -log10( median( ps ) )
	    lambda = lambda / -log10( 0.5 )

	    qs = -log10( ps )
	    expected = (1:length(ps))/(length(ps)+1)
	    expected = -log10( expected )

		what.to.plot = "-log10 P-value"
	} else {
		stop( "Only chisq and log10 scales are supported." )
	}

    cat( "Computing quantiles...\n" )
    lower.quantile = qbeta( 1-(1-CI)/2, (1:length(ps)), length(ps):1 )
    upper.quantile = qbeta( (1-CI)/2, (1:length(ps)), length(ps):1 )
	if( scale == "chisq" ) {
	    lower.quantile = qchisq( lower.quantile, df = df, lower.tail = F )
	    upper.quantile = qchisq( upper.quantile, df = df, lower.tail = F )
	} else if( scale == "log10" ) {
	    lower.quantile = -log10( lower.quantile )
	    upper.quantile = -log10( upper.quantile )
	} else {
		stop( sprintf( "Unrecognised scale \"%s\".", scale ))
	}
    # We're not much interested in anything that lies near the diagonal line.  So don't plot them.
    cat( "Thinning...\n" )
    if( thin ) {
		wDiag = which( qs <= upper.quantile & qs >= lower.quantile )
        wNotDiag = which( qs > upper.quantile | qs < lower.quantile )
        wDiag = sample( wDiag, min( 10000, length( wDiag ) ))
		qs = qs[ c( wDiag, wNotDiag )]
        expected = expected[ c( wDiag, wNotDiag ) ]
        lower.quantile = lower.quantile[ c( wDiag, wNotDiag ) ]
        upper.quantile = upper.quantile[ c( wDiag, wNotDiag ) ]
        o = order( expected )
        qs = qs[ o ]
        expected = expected[ o ]
        lower.quantile = lower.quantile[o]
        upper.quantile = upper.quantile[o]
    }
	if( is.null( xlab )) {
		xlab = sprintf( "Expected %s", what.to.plot )
	}
	if( is.null( ylab )) {
		ylab = sprintf( "Observed %s", what.to.plot )
	}

    cat( "Plotting...\n" )
    pch = rep( 20, length( qs )) # = '.'
	cex = rep( 0.6, length( qs ))
    pch[ which( qs > upper.quantile | qs < lower.quantile )] = 20
	cex[ which( qs > upper.quantile | qs < lower.quantile )] = 0.9
    if( plot.new ) {
		limits = c( max( pmin( expected, limit ) ), max( pmin( qs, limit ) ))
		plot( 0, 0, col = 'white', xlim = c( 0, limits[1] ), ylim = c( 0, limits[2] ), pch = 19, xlab = xlab, ylab = ylab, cex = cex, ... )
		grid()
	}
	points( pmin( expected, limit ), pmin( qs, limit ), pch = pch, xlab = xlab, ylab = ylab, cex = cex, ... )
	polygon(
		x = c( pmin( expected, limit ), rev( pmin( expected, limit ) )),
		y = c( pmin( lower.quantile, limit ), rev( pmin( upper.quantile, limit ) )),
		col = rgb( 0, 0, 0, 0.2 ),
		border = NA
	)
	#points( pmin( expected, limit ), pmin( lower.quantile, limit ), type = 'l', lty = 2 )
	#points( pmin( expected, limit ), pmin( upper.quantile, limit ), type = 'l', lty = 2 )
	if( !is.null( labels )) {
		wLabelled = which( !is.na( labels$label ))
		text( pmin( expected, limit )[wLabelled], pmin( qs, limit )[wLabelled], labels = labels$label[wLabelled], pos = labels$pos[wLabelled], cex = labels$cex[wLabelled] )
	}
    abline( a = 0, b = 1, col = "red", lty = 2 )
    return( lambda )
}

compute.ld <- function( haps, focus = 1:nrow(haps), variant.names = rownames( haps ) ) {
    # haps is a 0-1 matrix with L SNPs (in rows) and N haplotypes (in columns).
    # Since the values are 0, 1, we rely on the fact that a*b=1 iff a and b are 1.
    # assume focus.i contains l values in range 1..L
    focus.i = focus
    L = nrow( haps )
    l = length( focus.i )
    # p11 = lxL matrix.  i,jth entry is probability of 11 haplotype for ith focal SNP against jth SNP.
    focus.hap = haps[ focus.i, , drop = FALSE ]
    p11 <- ( focus.hap %*% t( haps )) / ncol( haps )
    # p1. = lxL matrix.  ith row is filled with the frequency of ith focal SNP.
    p1. <- matrix( rep( rowSums( focus.hap ) / ncol( haps ), L ), length( focus.i ), L, byrow = FALSE )
    # p.1 = lxL matrix.  jth column is filled with the frequency of jth SNP.
    frequency = rowSums( haps ) / ncol( haps )
    p.1 <- matrix( rep( frequency, length( focus.i ) ), length( focus.i ), L, byrow = TRUE )

    # Compute D
    D <- p11 - p1. * p.1

    # Compute D'
    denominator = pmin( p1.*(1-p.1), (1-p1.)*p.1 )
    wNeg = which(D < 0)
	if( length( wNeg ) > 0 ) {
		denominator[ wNeg ] = pmin( p1.*p.1, (1-p1.)*(1-p.1) )[wNeg]
	}
    denominator[ denominator == 0 ] = NA
    Dprime = D / denominator 

    # Compute correlation, this result should agree with cor( t(haps ))
    R = D / sqrt( p1. * ( 1 - p1. ) * p.1 * ( 1 - p.1 ))
    R[ which( ( p1. * ( 1 - p1. ) * p.1 * ( 1 - p.1 )) == 0 ) ] = NA
    if( !is.null( variant.names )) {
        rownames( D ) = rownames( Dprime ) = rownames( R ) = variant.names[ focus.i ]
        colnames( D ) = colnames( Dprime ) = colnames( R ) = variant.names
        names( frequency ) = variant.names
    }

    return( list( D = D, Dprime = Dprime, frequency = frequency, R = R ) ) ;
}

# Modified version of grid() which gives major divisions as thicker lines
# and minor divisions as thinner lines.
myGrid <- function(
    main.divisions = 1,
    sub.divisions = 4,
    main.col = "grey40",
    sub.col = "grey80",
    lty = 3,
    horiz = TRUE,
    vert = TRUE
) {
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

inthinnerate <- function(
	data,
	distance = "0.125cM+25kb",
	rank.column = NULL,
	genetic.map = sprintf( "%s/reference/mathgen.stats.ox.ac.uk/impute/2014-10-07-1000G_Phase3/genetic_map_chr#_combined_b37.txt", Sys.getenv( 'MG_PROJECTDIR' ) ),
	genes = sprintf( "%s/reference/genome-mysql.cse.ucsc.edu/2016-06-09/UCSC_hg19_2016-06-09_refGene.tsv", Sys.getenv( 'MG_PROJECTDIR' ) ),
	excluded.regions = c(),
	extra.args = "",
	picked.only = TRUE
) {
    a = tempfile( pattern = c( "in", "out" ), fileext = ".csv" )
	excluded.region.opts = ""
	if( length( excluded.regions ) > 0 ) {
		excluded.region.opts = sprintf( "-excl-range %s", paste( excluded.regions, collapse = " " ) )
	}
	if( picked.only ) {
		extra.args = paste( extra.args, '-suppress-excluded' )
	}
	data$chromosome = gsub( '^0', '', data$chromosome )
	if( is.null( rank.column )) {
		write.csv( data[,c( "rsid", "rsid", "chromosome", "position", "alleleA", "alleleB" ), ], file = a[1], row.names = F, quote = F )
	    system( sprintf( 'inthinnerator_v2.0.1 -min-distance %s -map %s -genes %s -g %s -o %s %s %s', distance, genetic.map, genes, a[1], a[2], excluded.region.opts, extra.args ) )
	} else {
		write.csv( data[,c( "rsid", "rsid", "chromosome", "position", "alleleA", "alleleB", gsub( '^-', '', rank.column )), ], file = a[1], row.names = F, quote = F )
		cmd = 
	    sprintf( 'inthinnerator_v2.0-dev -min-distance %s -map %s -genes %s -rank %s -rank-column %s -g %s -o %s %s %s', distance, genetic.map, genes, a[1], rank.column, a[1], a[2], excluded.region.opts, extra.args )
		print(cmd)
		system(cmd)
	}
	thinned = read.table( sprintf( '%s.000', a[2] ), hea=T, comment = '#', as.is = T )
	thinned = thinned[ which( thinned$result == "picked" ), ]
	print( head( thinned ) )
	result = data[ match( paste( thinned$rsid, fix_chromosome( thinned$chromosome ), thinned$position ), paste( data$rsid, fix_chromosome( data$chromosome ), data$position ) ), ]
	result$inthinnerator_result = thinned$result
	result$inthinnerator_pick_index = thinned$pick_index
	result$cM_from_start_of_chromosome = thinned$cM_from_start_of_chromosome
	result$region_lower_bp = thinned$region_lower_bp
	result$region_upper_bp = thinned$region_upper_bp
	result$nearest_gene = thinned$nearest_gene_in_region
	result$all_genes_in_region = thinned$all_genes_in_region
	return( result )
}

load.qctool.relatedness.matrix <- function( filename, sample_ids = NULL ) {
	X = read.csv( filename, comment = '#', head = T, as.is= T )
	rows = nrow(X)
	# This is lower diagonal (including diagonal)
	# so rows = N(N+1)/2 where N is number of samples
	# so N (N+1) = ( 2 * rows )
	# which is solved using the quadratic formula
	if( is.null( sample_ids )) {
		sample_ids = sort( unique( X$sample_1 ))
	}
	N = length( sample_ids )
	M = matrix( NA, N, N, dimnames = list( sample_ids, sample_ids ) )
	X$i = match( X$sample_1, sample_ids )
	X$j = match( X$sample_2, sample_ids )
	X = X[ which( !is.na( X$i ) & !is.na( X$j )), ]
	indices = as.matrix( X[, c( 'i', 'j' )] )
	M[ indices ] = X$value
	M[ upper.tri( M, diag = F )] = t(M)[ upper.tri( M, diag = F ) ]
	return( M )
}

################################################
# ggplot2 Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL, widths = NULL ) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
	if( is.null( widths )) {
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout) )))
	} else {
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), widths = widths )))
	}

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

load.project.callsets <- function( catalogue.filename, project, callsets = NULL ) {
    catalogue.db = dbConnect( dbDriver( "SQLite" ), args$catalogue )
    sql = sprintf( "SELECT * FROM ProjectCallsetView P INNER JOIN Callset C ON C.id == P.callset_id WHERE project_acronym == '%s'", args$project )
    if( !is.null( callsets ) ) {
        sql = paste(
            sql,
            sprintf(
                " AND callset_acronym IN ( '%s' )",
                paste( args$callsets, collapse = "','" )
            ),
            sep = ' '
        )
    }

    return( dbGetQuery( catalogue.db, sql ) )
}

compute.manhattan_plot_positions <- function(
        chromosome,
        position,
        seperation = 10000000
) {
        O = order( chromosome, position, decreasing = F )
        chromosome = chromosome[O]
        position = position[O]
        plot_pos = position
        plot_pos[-1] = plot_pos[-1] - plot_pos[1:( length( plot_pos ) - 1 )]
        plot_pos[ which( plot_pos < -100 ) ] = seperation
        plot_pos = cumsum( plot_pos )
        return( plot_pos )
}

SNPHWE <- function( obs_hets, obs_hom1, obs_hom2 ) {
   stopifnot(obs_hom1 >= 0 && obs_hom2 >= 0 && obs_hets >= 0)
   if( obs_hom1 + obs_hom2 + obs_hets == 0 ) {
	return( NA );
	}

   # total number of genotypes
   N <- obs_hom1 + obs_hom2 + obs_hets
   
   # rare homozygotes, common homozygotes
   obs_homr <- min(obs_hom1, obs_hom2)
   obs_homc <- max(obs_hom1, obs_hom2)

   # number of rare allele copies
   rare  <- obs_homr * 2 + obs_hets

   # Initialize probability array
   probs <- rep(0, 1 + rare)

   # Find midpoint of the distribution
   mid <- floor(rare * ( 2 * N - rare) / (2 * N))
   if ( (mid %% 2) != (rare %% 2) ) mid <- mid + 1

   probs[mid + 1] <- 1.0
   mysum <- 1.0

   # Calculate probablities from midpoint down 
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr

   while ( curr_hets >=  2)
      {
      probs[curr_hets - 1]  <- probs[curr_hets + 1] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0)  * (curr_homc + 1.0))
      mysum <- mysum + probs[curr_hets - 1]

      # 2 fewer heterozygotes -> add 1 rare homozygote, 1 common homozygote
      curr_hets <- curr_hets - 2
      curr_homr <- curr_homr + 1
      curr_homc <- curr_homc + 1
      }    

   # Calculate probabilities from midpoint up
   curr_hets <- mid
   curr_homr <- (rare - mid) / 2
   curr_homc <- N - curr_hets - curr_homr
   
   while ( curr_hets <= rare - 2)
      {
      probs[curr_hets + 3] <- probs[curr_hets + 1] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
      mysum <- mysum + probs[curr_hets + 3]
         
      # add 2 heterozygotes -> subtract 1 rare homozygtote, 1 common homozygote
      curr_hets <- curr_hets + 2
      curr_homr <- curr_homr - 1
      curr_homc <- curr_homc - 1
      }    
 
    # P-value calculation
    target <- probs[obs_hets + 1]

    #plo <- min(1.0, sum(probs[1:obs_hets + 1]) / mysum)

    #phi <- min(1.0, sum(probs[obs_hets + 1: rare + 1]) / mysum)

    # This assignment is the last statement in the fuction to ensure 
    # that it is used as the return value
    p <- min(1.0, sum(probs[probs <= target])/ mysum)
}

read.pslx <- function( filename ) {
	X = read.table( filename, skip = 5, header = F, as.is = T, sep = "\t" )
	colnames(X) = c(
		'match', 'mismatch', 'rep.match', 'Ns', 'Q_gap_count', 'Q_gap_bases', 'T_gap_count', 'T_gap_bases', 'strand',
		'Q_name', 'Q_size', 'Q_start', 'Q_end', 'T_name', 'T_size', 'T_start', 'T_end', 'blockCount', 'blockSizes',
		'qStarts', 'tStarts'
	)
	return(X)
}

draw_bracket <- function( x0, x1, x2 = NULL, y0, y1, y2 = NULL, vertical = FALSE, ... ) {
	if( vertical ) {
		stopifnot( !is.null( x2 ))
		segments(
			x0 = c( x0, x0, x1, x1 ),
			x1 = c( x1, x1, x1, x2 ),
			y0 = c( y0, y1, y0, (y0+y1)/2 ),
			y1 = c( y0, y1, y1, (y0+y1)/2 ),
			...
		)
	} else {
		stopifnot( !is.null( y2 ))
		segments(
			x0 = c( x0, x1, x0, (x0+x1)/2 ),
			x1 = c( x0, x1, x1, (x0+x1)/2 ),
			y0 = c( y0, y0, y1, y1 ),
			y1 = c( y1, y1, y1, y2 ),
			...
		)
	}
}

# betas is a dxN matrix of effect sizes, where N is a number of cohorts and d a number of effects per cohort
# Vs is a Nxdxd array of parameter variances + covariances.
compute.frequentist_meta_analysis <- function( betas, Vs, side = NULL ) {
    meta.V = matrix( 0, nrow = dim( betas )[1], ncol = dim( betas )[1], dimnames = list( dimnames( betas )[[1]], dimnames( betas )[[1]] ) )
    pre.meta.beta = matrix( 0, nrow = nrow( betas ), ncol = 1, dimnames = list( dimnames( betas )[[1]], "meta" ) )
    for( i in 1:ncol(betas) ) {
        W = solve( as.matrix( Vs[i,,] ))
        pre.meta.beta = pre.meta.beta + W %*% betas[,i,drop=F]
        meta.V = meta.V + W
    }
    meta.V = solve( meta.V )
    meta.se = sqrt( diag( meta.V ))
    meta.beta = meta.V %*% pre.meta.beta
	chi_squared = sum( meta.beta * pre.meta.beta )
    pvalue = pchisq( chi_squared, df = nrow( betas ), lower.tail = F )

    if( is.null( side ) ) {
        effect.pvalues = 2*pnorm( abs( meta.beta ), mean = 0, sd = meta.se, lower.tail = F )
    } else {
        side = sign(side)
        effect.pvalues = pnorm( (side * meta.beta), mean = 0, sd = meta.se, lower.tail = F )
    }
    return(
        list(
            meta_beta = meta.beta,
            meta_V = meta.V,
            meta_se = meta.se,
            meta_beta_lower = meta.beta - 1.96 * meta.se,
            meta_beta_upper = meta.beta + 1.96 * meta.se,
            meta_pvalue = pvalue,
            wald_pvalues = effect.pvalues,
            side = side
        )
    ) ;
}


load.hitplot.data <- function(
	db, range, analyses,
	verbose = FALSE,
	columns = list(
		bf = 'mean_bf',
		typed = 'typed',
		info = 'info'
	),
	tables = list(
		bf = 'CombinedMetaFiltered',
		typed = 'ImputeType'
	)
) {
	if( class( db ) == "character" ) {
		db = dbConnect( dbDriver( 'SQLite' ), db )
	}
	
	sql = paste(
		'SELECT analysis, BF.variant_id AS variant_id, BF.chromosome AS chromosome, BF.position AS position, BF.rsid AS rsid, CAST( "%s" AS FLOAT) AS mean_bf, "%s" AS typed, "%s" AS info',
		'FROM %s BF',
		'INNER JOIN %s Count',
		'ON Count.variant_id == BF.variant_id',
		'WHERE "prior_id" IN ("prior_id", 1) AND BF.chromosome == \'%s\' AND BF.position BETWEEN %d AND %d',
		sep = ' '
	)

	if( 'analysis' %in% names( opts )) {
		sql = sprintf(
			"%s AND BF.analysis IN ( '%s' )",
			sql,
			paste( analyses, collapse = "','" )
		)
	}

	sql = sprintf(
		sql, 
		opts$bfcolumn,
		opts$type_column,
		opts$info_column,
		opts$bftable,
		opts$typetable,
		range$chromosome,
		range$start,
		range$end
	)
	
	if( verbose ) {
		print( sql )
	}

	data = dbGetQuery( db, sql ) ;
	if( verbose ) {
		cat( "Data looks like:\n" )
		print( data[1:min(10,nrow(data)),] )
	}

	data$mean_bf = as.numeric( data$mean_bf )
	data = data[which( !is.na( data$mean_bf ) & data$mean_bf < Inf ), ]
	data$typed = data$typed
	
	data$typed[ data$typed == 0 ] = "Imputed variant"
	data$typed[ data$typed == 1 ] = "Omni 2.5M variant"
	data$typed[ data$analysis == 'hlaimp' ] = "Imputed HLA allele"
	
	{
		a = data$mean_bf / max( data$mean_bf, na.rm = T )
		sum_of_a = sum( a, na.rm = T )
		data$maller_et_al_posterior = a / sum_of_a
	}
	
	return( data )
}

get.focal.snp <- function( data, region = NULL, focal_snp_rsid = NULL ) {
	if( !is.null( focal_snp_rsid )) {
	   focal_snp = which( data$rsid == focal_snp_rsid ) ;
	   if( length( focal_snp ) < 1 ) {
	       cat( sprintf( "!! Focus SNP %s is not in the data!\n", focal_snp_rsid ) )
	       stop() ;
	   }
	   if( length( focal_snp ) > 1 ) {
	       cat( "!! More than one SNP matches focal SNP:\n" )
	       print( data[ focal_snp, ])
	       cat( "...I will take the first.\n" )
	       focal_snp = focal_snp[1,]
	   }
	} else {
	    bfs = data$mean_bf
		if( !is.null( region )) {
			bfs[ data$position < region[1] | data$position > region[2] ] = -1
		}
	    focal_snp = which.max( bfs )
	    rm( bfs )
	}
	hit = data[ focal_snp, ]
	return(hit)
}

annotate.ld <- function( data, haps, ld.focus, regularise = TRUE, verbose = FALSE ) {
	result = matrix( NA, ncol = 2, nrow = nrow( data ), dimnames = list( data$rsid, c( 'r_squared', 'absolute_Dprime' )))

	if( verbose ) {
		print( head( haps$data[, 1:10]))
	}
	if( length( ld.focus ) > 1 ) {
		cat( "More than one haplotype SNP is at the focal position, I'll use the first.\n" )
		ld.focus = ld.focus[1]
	}
	if( length( ld.focus ) ==  1 ) {
		cat( "Computing r^2 and D'...\n" )
		# For LD computation, we prevent bad behaviour by adding dummy columns to the haplotype file to
		# ensure all four haplotypes are represented at least once.
		# (Note this can change LD considerably at rare SNPs.)
		if( regularise ) {
			numberOfHaps = ncol( haps$data )
			additionalHaps = cbind( matrix( 0, ncol = 2, nrow = nrow( haps$data )), matrix( 1, ncol = 2, nrow = nrow( haps$data )))
			additionalHaps[ld.focus, c( 2, 4 )] = c( 1, 0 )
			haps$data = cbind( haps$data, additionalHaps )
		}
		ld.matrices = compute.ld( haps$data, focus = ld.focus )
		ld = data.frame(
			haps$variant,
			r = as.numeric( ld.matrices$R ),
			absolute_Dprime = as.numeric( ld.matrices$Dprime )
		) ;
		if( nrow( ld ) > 0 ) {
			result[,'r_squared'] = ld$r[ match( data$position, ld$position ) ]^2
			result[,'absolute_Dprime'] = ld$absolute_Dprime[ match( data$position, ld$position ) ]
		} else {
			result[,'r_squared'] = NA
			result[,'absolute_Dprime'] = NA
		}
	} else {
		result[,'r_squared'] = NA
		result[,'absolute_Dprime'] = NA
		ld = data.frame()
	}
	return( result )
}

query <- function(
	db,
	select,
	from,
	join = NULL,
	where = NULL,
	limit = NULL,
	explain = FALSE,
	verbose = FALSE
) {
	escape = sprintf( "`%s`", select )
	w = grep( "[*.]", select )
	if( length( w ) > 0 ) {
		escape[ w ] = select[w]
	}
	
	select.sql = sprintf(
		"SELECT %s",
		paste( escape, sep = ', ' )
	)
	from.sql = sprintf( "FROM `%s`", from )
	join.sql = ''
	where.sql = ''
	limit.sql = ''
	explain.sql = ''
	if( !is.null( join.sql )) {
		join.sql = join
	}
	if( !is.null( where.sql )) {
		where.sql = where
	}
	if( !is.null( limit )) {
		limit.sql = sprintf( "LIMIT %d", limit )
	}
	if( explain ) {
		explain.sql = "EXPLAIN QUERY PLAN"
	}
	sql = paste(
		sep = ' ',
		explain.sql,
		select.sql,
		from.sql,
		join.sql,
		where.sql,
		limit.sql
	)
	if( verbose ) {
		cat( "query=\"", sql, "\".\n", sep = '' )
	}
	dbGetQuery( db, sql )
}

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
 
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
 
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
 
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
 
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
 
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
 
  x1 <- switch(pos,
    topleft     =x[1] + sw, 
    left        =x[1] + sw,
    bottomleft  =x[1] + sw,
    top         =(x[1] + x[2])/2,
    center      =(x[1] + x[2])/2,
    bottom      =(x[1] + x[2])/2,
    topright    =x[2] - sw,
    right       =x[2] - sw,
    bottomright =x[2] - sw)
 
  y1 <- switch(pos,
    topleft     =y[2] - sh,
    top         =y[2] - sh,
    topright    =y[2] - sh,
    left        =(y[1] + y[2])/2,
    center      =(y[1] + y[2])/2,
    right       =(y[1] + y[2])/2,
    bottomleft  =y[1] + sh,
    bottom      =y[1] + sh,
    bottomright =y[1] + sh)
 
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
 
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


# Utility function for creating plot layouts
# Given an existing layout matrix (layout)
# And a single panel in that layout
# This expands the layout so that that panel has the layout of sub,
# Consecutive numbering is taken care of.
expand.layout = function(
	layout,
	panel,
	sublayout
) {
	sub = sublayout
	w = which( layout == panel, arr.in = T )
	x = as.integer(w[,1])
	y = as.integer(w[,2])
#	print(w)
	add.x = nrow(sub)
	add.y = ncol(sub)
	stopifnot( nrow(w) == 1 )
	new.layout = matrix(
		NA,
		nrow = nrow( layout ) + nrow( sub ) - 1,
		ncol = ncol( layout ) + ncol( sub ) - 1
	)
	wGreater = which( layout > panel )
	# Make space in the numbering
	if( length( wGreater ) > 0 ) {
		layout[wGreater] = layout[wGreater] + max( sub ) - 1
	}

	M.left = matrix(
		0,
		nrow = nrow( layout ) + nrow( sub ) - 1,
		ncol = nrow( layout )
	)
	M.right = matrix(
		0,
		nrow = ncol( layout ),
		ncol = ncol( layout ) + ncol( sub ) - 1
	)

	if(x>1) {
		diag(M.left)[1:x] = 1
	}
	if(y>1) {
		diag(M.right)[1:y] = 1
	}
	M.left[x:(x+add.x-1),x] = 1
	M.right[y,y:(y+add.y-1)] = 1
    if (x < nrow(layout)) {
        M.left[matrix(c((x + add.x):nrow(M.left), (x + 1):nrow(layout)), ncol = 2)] = 1
    }
    if (y < ncol(layout)) {
        M.right[matrix(c((y + 1):ncol(layout), (y + add.y):ncol(M.right)), ncol = 2)] = 1
    }
    result = M.left %*% layout %*% M.right
    result[x:(x + add.x - 1), y:(y + add.y - 1)] = sub + panel - 1
    result
}

load.qc.result <- function( callset, snp.db = NULL ) {
	if( is.null( snp.db )) {
		library( dplyr )
		library( dbplyr )
		snp.db = DBI::dbConnect( RSQLite::SQLite(), sprintf( "%s/qc/include_exclude/detail/snp_qc_metrics_aligned.sqlite", ANALYSISDIR ))
	}
	pop = gsub( "_GWAS-2.5M_b37", "", callset, fixed = T )
	
	cat( sprintf( "Loading SNP QC results for %s...\n", pop ))
	
	sql = "
	SELECT chromosome, position, rsid, alleleA, alleleB,
	CPR.platform_status AS platform_status,
	`[POP]_freq` AS alleleB_frequency,
	`[POP]_missing` AS missingness,
	`[POP]_plate` AS plate_test_pvalue,
	HWE.`[POP]` AS HWE_test_pvalue,
	RC.`[POP]` AS recall_test_pvalue,
	CASE WHEN R.`[POP]_result` IS NULL THEN NULL WHEN R.`[POP]_result` == 'ok' THEN 'ok' ELSE 'excluded:' || R.`[POP]_result` END AS result,
	CASE WHEN CPR.duplicated_status == 'duplicate' THEN 1 ELSE 0 END AS excluded_as_duplicate,
	CASE WHEN CPR.cluster_plot_status == 'ok' THEN 0 ELSE 1 END AS `superset:bad_cluster_plot`,
	CASE WHEN CPR.result == 'ok' AND CPR.platform_status == 'both' THEN 1 ELSE 0 END AS `superset:included`
	FROM Result R
	INNER JOIN CrossPopResultView CPR ON CPR.variant_id == R.variant_id
	INNER JOIN Plate P ON P.variant_id == R.variant_id
	INNER JOIN HWE ON HWE.variant_id == R.variant_id
	INNER JOIN RecallComparison RC ON RC.variant_id == R.variant_id
	"
	qc.result = DBI::dbGetQuery( snp.db, gsub( '[POP]', pop,  fixed = T, sql ))

	if( callset == "Kenya_GWAS-2.5M_b37" ) {
		qc.result$result[ qc.result$platform_status == 'v8_only' ] = NA
	}

	return( qc.result )
}

load.sample.stats <- function( type = "superset1_snps" ) {
	library( RSQLite )
	if( type != "superset1_snps" ) {
		filename = sprintf( "%s/qc/stats/sample_stats.sqlite", ANALYSISDIR )
		sql = "SELECT * FROM SampleStatsUnalignedView S
		INNER JOIN Analysis A ON A.id == S.analysis_id
		WHERE analysis NOT LIKE '%-Kumasi%' AND analysis NOT LIKE '%-Navrongo%'
		AND A.chunk == 'autosomes'"
		
	} else {
		filename = sprintf( "%s/qc/stats/stats.sqlite", ANALYSISDIR )
		sql = "SELECT * FROM SampleStatsView S"
		
	}
	stats.db = dbConnect( dbDriver( "SQLite" ), filename )
	stats = dbGetQuery( stats.db, sql )
	stats$pop = sapply( stats$analysis, function(a) { strsplit( a, split = "_", fixed = T )[[1]][1] } )
	stats$pop = factor(
		stats$pop,
		levels = c(
			"Gambia", "Mali", "BurkinaFaso", "Ghana", "Nigeria", "Cameroon",
			"Malawi", "Tanzania", "Kenya",
			"Vietnam", "PNG"
		)
	)
	props = dbGetQuery( stats.db, sprintf( "SELECT * FROM AnalysisPropertyView WHERE analysis_id == %d", stats$analysis_id[1] ))
	print( props )
	return( stats )
}
