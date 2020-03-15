import csv, os
from glob import glob
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

files = {
	"10X": glob( "%s/projects/wgs_10x/data/raw/11589MCpool01-N__JK_HV31_3/11589MCpool01-N__JK_HV31_3_S8_*.fastq.gz" % ROOT ),
	"PB-CCS": glob( "%s/projects/pacbio/data/raw/ccs/combined/JK_HV31.ccs.fastq.gz" % ROOT )
	"PB-CLR": glob( "%s/projects/pacbio/data/raw/clr/CLR_190625/m64016_190625_101024.subreads.fastq.gz" % ROOT ),
	"ONT": glob( "/projects/nanopore/data/raw/WTON00039*/trimmed/fastq-trimmed-all.fastq.gz" % ROOT )
}

ks = [ 15, 17, 23, 31 ]

# We ignore bases with lower quality than this:
# NB. Illumina base qualities are binned, the values are 
# J   F   A   <   7   2   -   (   #
# 41  37  32  27  22  17  12  7   2
minBaseQualityChar = "?" # quality=30, or P(error) = 0.001. See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm

outDir = "%s/shared/analysis/kmer/regions/distribution" % ROOT

rule all:
	input:
	[ 
		"{dir}/{key}.k={k}.histogram".format( dir = outDir, key = key, k = k )
		for key in files.keys()
		for k in ks
	]

rule countKmers:
	input:
		fastq = lambda w: files[ w.key ]
	output:
		jf = "%s/{key}.k={k}.jf" % outDir
	params:
		minQualChar = minBaseQualityChar
	shell: """
		jellyfish -m {wildcards.k} -s 3.2G --min-qual-char="{params.minQualChar}" --both-strands -o {output.jf}
	"""

rule computeHistogram:
	input:
		"%s/{key}.k={k}.jf" % outDir
	output:
		"%s/{key}.k={k}.histogram" % outDir
	shell: """
		jellyfish histo {input} -o {output}
	"""

