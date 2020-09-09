import csv, os
from glob import glob
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

tools = {
	"Rscript": "/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript",
	"genomescope": "~/Projects/Software/3rd_party/genomescope2.0/genomescope.R"
}

files = {
	"NA12878": glob( "%s/shared/data/NA12878/reads/ERR3239334_[12].fastq.gz" % ROOT ),
	"10X": (
		glob( "%s/projects/wgs_10x/data/raw/11589MCpool01-N__JK_HV31_3/11589MCpool01-N__JK_HV31_3_S8_L00*R1_001.nobarcode.fastq.gz" % ROOT ) +
		glob( "%s/projects/wgs_10x/data/raw/11589MCpool01-N__JK_HV31_3/11589MCpool01-N__JK_HV31_3_S8_L00*R2_001.fastq.gz" % ROOT )
	),
	"10X:JKOBB008": glob( "%s/projects/wgs_10x/data/raw/11589MCpool01-N__JKOBB008/11589MCpool01-N__JKOBB008_S5_L00*R[12]_001.fastq.gz" % ROOT ),
	"PB-CCS": glob( "%s/projects/pacbio/data/raw/ccs/combined/JK_HV31.ccs.fastq.gz" % ROOT ),
	"PB-CLR": glob( "%s/projects/pacbio/data/raw/clr/CLR_190625/m64016_190625_101024.subreads.fastq.gz" % ROOT ),
	"ONT": glob( "%s/projects/nanopore/data/raw/WTON00039*/trimmed/fastq-trimmed-all.fastq.gz" % ROOT ),
	"MGI": glob( "%s/projects/mgi/HV31/L*/V300037759_L0[1234]_[12]_[12].fq.gz" % ROOT ),
	"MGI:NA12878": glob( "%s/projects/mgi/NA12878/L*/V300037759_L0[1234]_[34]_[12].fq.gz" % ROOT ),
	"stLFR": (
		glob( "%s/projects/stLFR/HV31/L0[12]/V300037763_L0[12]_read_1.fq.gz" % ROOT ) +
		glob( "%s/projects/stLFR/HV31/L0[12]/V300037763_L0[12]_read_2_nobarcode.fq.gz" % ROOT )
	),
	"stLFR:NA12878": (
		glob( "%s/projects/stLFR/NA12878/L0[34]/V300037763_L0[34]_read_1.fq.gz" % ROOT ) +
		glob( "%s/projects/stLFR/NA12878/L0[34]/V300037763_L0[34]_read_2_nobarcode.fq.gz" % ROOT )
	),
	"illumina_pcr_free": glob( "%s/projects/illumina/data/raw/pcr_free/200515_A00711_0200_AHLTJ3DRXX/WTCHG_788803_[12].fastq.gz" % ROOT ),
	"CoolMPS": glob( "%s/projects/mgi_coolmps/HV31/L0[1234]/V300037630_L0[1234]_[12]_[12].fq.gz" % ROOT ),
	"CoolMPS:NA12878": glob( "%s/projects/mgi_coolmps/NA12878/L0[1234]/V300037630_L0[1234]_[34]_[12].fq.gz" % ROOT )
}

merges = {}
merges[ '10X+MGI+illumina_pcr_free+stLFR' ] = [ '10X', 'MGI', 'illumina_pcr_free', 'stLFR' ]
merges[ '10X+MGI+illumina_pcr_free+stLFR+CoolMPS' ] = [ '10X', 'MGI', 'illumina_pcr_free', 'stLFR', 'CoolMPS' ]
merges[ 'MGI+illumina_pcr_free+stLFR+CoolMPS' ] = [ '10X', 'MGI', 'illumina_pcr_free', 'stLFR', 'CoolMPS' ]
merges[ 'PB-CCS+MGI+illumina_pcr_free+stLFR+CoolMPS' ] = [ 'PB-CCS', 'MGI', 'illumina_pcr_free', 'stLFR', 'CoolMPS' ]
merges[ 'MGI+illumina_pcr_free' ] = [ 'MGI', 'illumina_pcr_free' ]

merges[ 'PB-CCS+10X+MGI+illumina_pcr_free+stLFR' ] = [ 'PB-CCS', '10X', 'MGI', 'illumina_pcr_free', 'stLFR' ]
merges[ 'PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS' ] = [ 'PB-CCS', '10X', 'MGI', 'illumina_pcr_free', 'stLFR', 'CoolMPS' ]
merges[ 'PB-CCS+10X+MGI' ] = [ 'PB-CCS', '10X', 'MGI' ]
merges[ 'PB-CCS+MGI' ] = [ 'PB-CCS', 'MGI' ]
merges[ '10X+MGI' ] = [ '10X', 'MGI' ]

wildcard_constraints:
	platform = '|'.join( files.keys() ),
	mergedPlatform = '|'.join( merges.keys() ).replace( '+', '[+]' )

for name in files.keys():
	if len( files[name] ) == 0:
		raise Exception( "Expected at least one fastq file for \"%s\"." % name )


print( "files: %s" % ' '.join( files.keys() ))
print( "merges: %s" % ' '.join( merges.keys() ))

# kmer sizes:
# Flye uses 15 in standard mode, 17 in CCS mode and 31 in contig assembly mode
# Jia-Yuan used 18 in his histograms
ks = [ 14, 16, 17, 18, 22, 23, 31 ]

# We ignore bases with lower quality than this:
# NB. Illumina base qualities are binned, the values are 
# J   F   A   <   7   2   -   (   #
# 41  37  32  27  22  17  12  7   2
minBaseQualityChars = { # See https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
	"40": "I",            # quality=40, or P(error) = 0.0001.  (NB. 41 is the max Illumina outputs)
	"30": "?",            # quality=30, or P(error) = 0.001
	"20": "5",            # quality=20, or P(error) = 0.01
	"10": "+"             # quality=10, or P(error) = 0.1
}

# These thresholds are based on looking at plots and locate a set of
# kmers frequencies that should include most 'solid' (i.e. not error, not > 2 copy number) kmers
# the intention here is to make these sets large enough to contain the true solid kmers, so they can be later further filtered
import math
solidKmerThresholds = {
	'PB-CCS+10X+MGI.k=16.bq=20': [ 20, math.floor( 80*1.5 ) ],
	'PB-CCS+10X+MGI.k=17.bq=20': [ 20, math.floor( 78*1.5 ) ],
	'PB-CCS+10X+MGI.k=18.bq=20': [ 19, math.floor( 76*1.5 ) ],
	'PB-CCS+10X+MGI.k=22.bq=20': [ 16, math.floor( 71*1.5 ) ],
	'PB-CCS+10X+MGI.k=23.bq=20': [ 16, math.floor( 70*1.5 ) ],
	'PB-CCS+10X+MGI.k=31.bq=20': [ 12, math.floor( 61*1.5 ) ],
	'PB-CCS+MGI.k=16.bq=20': [ 7, math.floor( 54*1.5 ) ],
	'PB-CCS+MGI.k=17.bq=20': [ 7, math.floor( 53*1.5 ) ],
	'PB-CCS+MGI.k=18.bq=20': [ 6, math.floor( 52*1.5 ) ],
	'PB-CCS+MGI.k=22.bq=20': [ 5, math.floor( 49*1.5 ) ],
	'PB-CCS+MGI.k=23.bq=20': [ 5, math.floor( 48*1.5 ) ],
	'PB-CCS+MGI.k=31.bq=20': [ 4, math.floor( 42*1.5 ) ],
	'10X+MGI.k=16.bq=20': [ 19, math.floor( 68*1.5 ) ],
	'10X+MGI.k=17.bq=20': [ 17, math.floor( 67*1.5 ) ],
	'10X+MGI.k=18.bq=20': [ 16, math.floor( 66*1.5 ) ],
	'10X+MGI.k=22.bq=20': [ 14, math.floor( 61*1.5 ) ],
	'10X+MGI.k=23.bq=20': [ 13, math.floor( 60*1.5 ) ],
	'10X+MGI.k=31.bq=20': [ 11, math.floor( 52*1.5 ) ]
} ;

def makeMinQualOption( bq ):
	if bq == "any":
		return ""
	else:
		return '--min-qual-char "%s"' % minBaseQualityChars[bq]

outDir = "distribution"

rule all:
	input:
		[ 
		"{dir}/images/{platform}.k={k}.bq={bq}.histogram.pdf".format( dir = outDir, platform = platform, k = k, bq = bq )
		for platform in list( files.keys() ) + list( merges.keys())
		for k in ks
		for bq in list( minBaseQualityChars.keys() ) + [ 'any' ]
	] + [
		"{dir}/{platform}/{platform}.k={k}.bq={bq}.genomescope/{platform}.k={k}.bq={bq}_summary.txt".format( dir = outDir, platform = platform, k = k, bq = bq )
		for platform in list( files.keys() ) + list( merges.keys())
		for k in ks
		for bq in list( minBaseQualityChars.keys() ) + [ 'any' ]
	] + [
		"%s/histograms.csv" % outDir
	] + [
		"%s/PB-CCS+10X+MGI/PB-CCS+10X+MGI.k=22.bq=20.jf" % outDir,
		"%s/PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS/PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS.k=22.bq=20.jf" % outDir,
		"%s/PB-CCS+10X+MGI/PB-CCS+10X+MGI.k=31.bq=20.jf" % outDir,
		"%s/PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS/PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS.k=31.bq=20.jf" % outDir
	]

rule countKmers:
	input:
		fastq = lambda w: files[ w.platform ]
	output:
		jf = temp( "%s/{platform}/{platform}.k={k}.bq={bq}.jf" % outDir )
	params:
		inputs = (lambda w: ' '.join( files[w.platform] )),
		minQualOption = (lambda w: makeMinQualOption( w.bq )),
		threads = 10
	shell: """
		echo {params.inputs}
		zcat {params.inputs} | jellyfish count --threads {params.threads} -m {wildcards.k} -s 3G {params.minQualOption} -C -o {output.jf}.tmp /dev/fd/0
		mv {output.jf}.tmp {output.jf}
	"""

rule countMergedKmers:
	input:
		fastq = lambda w: [ item for platform in merges[ w.mergedPlatform ] for item in files[ platform ]]
	output:
		jf = temp( "%s/{mergedPlatform}/{mergedPlatform}.k={k}.bq={bq}.jf" % outDir )
	params:
		inputs = (lambda w, output, input: ' '.join( input.fastq )),
		minQualOption = (lambda w: makeMinQualOption( w.bq )),
		threads = 10
	shell: """
		echo {params.inputs}
		zcat {params.inputs} | jellyfish count --threads {params.threads} -m {wildcards.k} -s 4G {params.minQualOption} -C -o {output.jf}.tmp /dev/fd/0
		mv {output.jf}.tmp {output.jf}
	"""


# Merging is not always working because the hashes get different sizes.
# So I reverted to directly computing as above.
#rule mergeKmers:
#	input:
#		jf = lambda w: [ "{dir}/{platform}/{platform}.k={k}.bq={bq}.jf".format(
#				dir = outDir,
#				platform = platform,
#				k = w.k,
#				bq = w.bq
#			) for platform in merges[ w.mergedPlatform ]
#		]
#	output:
#		jf = temp( "%s/{mergedPlatform}/{mergedPlatform}.k={k}.bq={bq}.jf" % outDir )
#	shell: """
#		jellyfish merge -o {output.jf} {input.jf}
#	"""

rule computeHistogram:
	input:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.jf" % outDir
	output:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.histogram" % outDir
	shell: """
		jellyfish histo {input} -o {output}.tmp
		mv {output}.tmp {output}
	"""

def getSampleName( platform ):
	bits = platform.split( ':' )
	if len( bits ) == 2:
		return { 'sample': bits[1], 'platform': bits[0] }
	elif platform == "NA12878":
		return { 'sample': "NA12878", 'platform': "Novaseq" }
	else:
		return { 'sample': 'HV31', 'platform': platform }

rule plotHistogram:
	input:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.histogram" % outDir
	output:
		"%s/images/{platform_or_merge}.k={k}.bq={bq}.histogram.pdf" % outDir
	params:
		script = "pipelines.symlink/scripts/plot_kmer_histogram.R",
		sample = lambda w: getSampleName( w.platform_or_merge )['sample'],
		platformName = lambda w: getSampleName( w.platform_or_merge )['platform']
	shell: """
		/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript --vanilla {params.script} {input} {output} {params.sample} {params.platformName} {wildcards.k} {wildcards.bq}
	"""

rule collectHistograms:
	input:
		histograms = [ ( "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.histogram" % outDir ).format(
				platform_or_merge = p,
				k = k,
				bq = bq
			)
			for p in [ '10X', 'MGI', 'CoolMPS', 'stLFR', 'illumina_pcr_free', 'PB-CCS' ]
			for k in [ 23, 31 ]
			for bq in list( minBaseQualityChars.keys() ) + [ 'any' ]
		]
	output:
		csv = "%s/histograms.csv" % outDir
	params:
		platform_or_merges = [ '10X', 'MGI', 'CoolMPS', 'stLFR', 'illumina_pcr_free', 'PB-CCS' ],
		ks = [ 23, 31 ],
		bqs = list( minBaseQualityChars.keys() ) + [ 'any' ],
		tmp = "%s/histograms.csv.tmp" % outDir
	run:
		shell("""rm -f {output.csv}""" )
		shell("""rm -f {params.tmp}""" )
		shell( """echo sample,platform,k,bq,coverage,count > {params.tmp}""" )
		for platform_or_merge in params.platform_or_merges:
			sample = getSampleName( platform_or_merge )
			platformName = sample['platform']
			sample = sample['sample']
			for k in params.ks:
				for bq in params.bqs:
					filename = "%s/%s/%s.k=%s.bq=%s.histogram" % ( outDir, platform_or_merge, platform_or_merge, k, bq )
					shell( """cat """ + filename + """ | awk "{{printf( \\"%s,%s,%s,%s,%s,%s\\n\\", \\"{sample}\\", \\"{platform}\\", \\"{k}\\", \\"{bq}\\", \$1, \$2 )}}" >> {params.tmp}""" )
		shell( """mv {params.tmp} {output.csv}""" )

rule run_genomescope:
	input:
		histogram = "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.histogram" % outDir
	output:
		summary = "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.genomescope/{platform_or_merge}.k={k}.bq={bq}_summary.txt" % outDir
	params:
		dir = "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.genomescope" % outDir
	shell: """
		{tools[Rscript]} --vanilla {tools[genomescope]} -i {input.histogram} -o {params.dir} -p 2 -k {wildcards.k} -n {wildcards.platform_or_merge}.k={wildcards.k}.bq={wildcards.bq}
	"""

rule extractSolidKmers:
	input:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.jf" % outDir
	output:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_kmer_superset.txt.gz" % outDir
	params:
		lower = lambda w: solidKmerThresholds.get("%s.k=%s.bq=%s" % (w.platform_or_merge, w.k, w.bq), [0,1000] )[0],
		upper = lambda w: solidKmerThresholds.get("%s.k=%s.bq=%s" % (w.platform_or_merge, w.k, w.bq), [0,1000] )[1],
		output = "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_kmer_superset.txt" % outDir
	shell: """
		jellyfish dump -c --lower-count {params.lower} --upper-count {params.upper} -o {params.output} {input}
		gzip {params.output}
	"""

rule extractSolidAndHighCopyKmers:
	input:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.jf" % outDir
	output:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_and_high_copy_kmer_superset.txt.gz" % outDir
	params:
		lower = lambda w: solidKmerThresholds.get( "%s.k=%s.bq=%s" % (w.platform_or_merge, w.k, w.bq), [0] )[0],
		output = "%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_and_high_copy_kmer_superset.txt" % outDir
	shell: """
		jellyfish dump -c --lower-count {params.lower} -o {params.output} {input}
		gzip {params.output}
	"""

rule extractKmersForMhap:
	input:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_kmer_superset.txt.gz" % outDir
	output:
		"%s/{platform_or_merge}/{platform_or_merge}.k={k}.bq={bq}.solid_kmer_superset_mhap.txt.gz" % outDir
	shell: """
		./pipelines.symlink/scripts/make_mhap_kmer_file {input} | gzip -c > {output}
	"""
