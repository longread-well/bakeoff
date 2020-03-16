import csv, os
from glob import glob
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

files = {
	"10X": glob( "%s/projects/wgs_10x/data/raw/11589MCpool01-N__JK_HV31_3/11589MCpool01-N__JK_HV31_3_S8_*.fastq.gz" % ROOT ),
	"PB-CCS": glob( "%s/projects/pacbio/data/raw/ccs/combined/JK_HV31.ccs.fastq.gz" % ROOT ),
	# CLR is not working, presumably because there's no base quality score in the FASTQ file.
	# "PB-CLR": glob( "%s/projects/pacbio/data/raw/clr/CLR_190625/m64016_190625_101024.subreads.fastq.gz" % ROOT ),
	"ONT": glob( "%s/projects/nanopore/data/raw/WTON00039*/trimmed/fastq-trimmed-all.fastq.gz" % ROOT )
}

# kmer sizes:
# Flye uses 15 in standard mode, 17 in CCS mode and 31 in contig assembly mode
# Jia-Yuan used 18 in his histograms
ks = [ 15, 17, 18, 23, 31 ]

# We ignore bases with lower quality than this:
# NB. Illumina base qualities are binned, the values are 
# J   F   A   <   7   2   -   (   #
# 41  37  32  27  22  17  12  7   2
minBaseQualityChars = [       # See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
#	"I",                  # quality=40, or P(error) = 0.0001.  (NB. 41 is the max Illumina outputs)
#	"?",                  # quality=30, or P(error) = 0.001
	"5",                  # quality=20, or P(error) = 0.01
#	"+"                   # quality=10, or P(error) = 0.1
]

outDir = "%s/shared/analysis/kmer/distribution" % ROOT

rule all:
	input:
		[ 
		"{dir}/{platform}/{platform}.k={k}.bq={bq}.histogram".format( dir = outDir, platform = platform, k = k, bq = bq )
		for platform in files.keys()
		for k in ks
		for bq in minBaseQualityChars
	] + [ 
		"{dir}/images/{platform}.k={k}.bq={bq}.histogram.pdf".format( dir = outDir, platform = platform, k = k, bq = bq )
		for platform in files.keys()
		for k in ks
		for bq in minBaseQualityChars
	]

rule countKmers:
	input:
		fastq = lambda w: files[ w.platform ]
	output:
		jf = "%s/{platform}/{platform}.k={k}.bq={bq}.jf" % outDir
	params:
		inputs = lambda w: ' '.join( files[w.platform] )
	shell: """
		zcat {params.inputs} | jellyfish count -m {wildcards.k} -s 3G --min-qual-char="{wildcards.bq}" -C -o {output.jf} /dev/fd/0
	"""

rule computeHistogram:
	input:
		"%s/{platform}/{platform}.k={k}.bq={bq}.jf" % outDir
	output:
		"%s/{platform}/{platform}.k={k}.bq={bq}.histogram" % outDir
	shell: """
		jellyfish histo {input} -o {output}
	"""

rule plotHistogram:
	input:
		"%s/{platform}/{platform}.k={k}.bq={bq}.histogram" % outDir
	output:
		"%s/images/{platform}.k={k}.bq={bq}.histogram.pdf" % outDir
	shell: """
		echo 'X = read.table( "{input}", hea=F, as.is=T ); \
pdf( file = "{output}", width = 8, height = 4 ); \
par( mfrow = c( 1, 2 ))
plot( pmin( X[,1], 100 ), X[,2], type = "l", xlab = "Coverage", ylab = "Count" ); \
grid(); \
legend( "topright", c( sprintf( "Platform: %s", "{wildcards.platform}" ), sprintf( "k: %s", "{wildcards.k}" ), sprintf( "bq >= %s", {wildcards.bq} )), bty = "n" ) \
m = max( X[ X[,1] > 10 ) \
plot( pmin( X[,1], 100 ), X[,2], type = "l", xlab = "Coverage", ylab = "Count (zoomed)", ylim = c( 0, m * 1.5 )) ; \
grid(); \
dev.off()' \
| /apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/R --vanilla
	"""

