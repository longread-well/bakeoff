import csv, os
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

def loadRegions( filename = '../resources/regions.tsv' ):
	regions = []
	with open( filename, newline = '' ) as file:
		regionReader = csv.DictReader( file, delimiter = '\t' )
		for region in regionReader:
			print( region )
			regions.append( dict( region ))
	return regions

regions = loadRegions()

basedir = "%s/projects" % ROOT
builds = {
	"GRCh37": '%s/reference/GRCh37/hs37d5/hs37d5.fa' % basedir,
	"GRCh38": '%s/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa' % basedir
}
tools = {
	'samtools': '/apps/well/samtools/1.8-gcc5.4.0/bin/samtools',
	'clustalo': '/apps/well/clustal-omega/1.2.1/bin/clustalo',
	'mummer': '/apps/well/mummer/3.23-2-gcc4.7.2',
	'selfmap': '/users/kwiatkowski/gav/Projects/Software/usr/bin/selfmap_v2.1-dev',
	'selfmap.R': '%s/shared/analysis/bakeoff-scripts/analysis/kmer/scripts/selfmap.R' % ROOT,
	'Rscript': '/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript'
}

def getRegionFastaPath( region, build ):
	return "regions/sequence/{region}.{build}.fa".format(
		region = region['acronym'],
		build = build
	)

def getRegionDefinition( acronym, build ):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	return result[0]

assemblyMethods = [
	'canu-racon-medaka',
	'wtdbg2-racon-medaka',
	'flye-racon-medaka'
]

rule all:
	input:
		fasta = [ getRegionFastaPath( region, build ) for build in builds.keys() for region in regions ],
		delta = [ "regions/comparison/{acronym}.GRCh37-GRCh38.delta".format( acronym = region['acronym']) for region in regions ],
		selfmap = [ "regions/comparison/{acronym}.GRCh37-GRCh38.k=100.selfmap".format( acronym = region['acronym']) for region in regions ],
		dotplot = [ "regions/images/{acronym}.GRCh37-GRCh38.k=100.dotplot.pdf".format( acronym = region['acronym']) for region in regions ],
		assemblyDelta= [ "regions/comparison/{acronym}.{build}-{method}.delta".format( build = build, method = method, acronym = region['acronym']) for build in builds.keys() for method in assemblyMethods for region in regions if region['acronym'] == 'TRA' ],
		assemblyDotplot = [ "regions/images/{acronym}.{build}-{method}.k=100.dotplot.pdf".format( build = build, method = method, acronym = region['acronym']) for build in builds.keys() for method in assemblyMethods for region in regions if region['acronym'] == 'TRA' ]

wildcard_constraints:
	build='G[A-Za-z0-9]*'

rule extractFasta:
	input:
		fasta = lambda w: builds[ w.build ]
	output:
		fasta = "regions/sequence/{acronym}.{build}.fa",
		md5sum = "regions/sequence/{acronym}.{build}.fa.md5"
	params:
		chromosome = lambda x: getRegionDefinition( x.acronym, x.build )['chromosome'],
		start = lambda x: getRegionDefinition( x.acronym, x.build )['start'],
		end = lambda x: getRegionDefinition( x.acronym, x.build )['end']
	shell:
		"""
			samtools faidx {input.fasta} {params.chromosome}:{params.start}-{params.end} > {output.fasta}
			echo {output.fasta} `tail -n +2 {output.fasta} | md5sum | cut -d' ' -f1` > {output.md5sum}
		"""

	
rule nucmerCompare:
	input:
		b37 = "regions/sequence/{acronym}.GRCh37.fa",
		b38 = "regions/sequence/{acronym}.GRCh38.fa"
	output:
		delta = "regions/comparison/{acronym}.GRCh37-GRCh38.delta",
		coords = "regions/comparison/{acronym}.GRCh37-GRCh38.coords",
		tiling = "regions/comparison/{acronym}.GRCh37-GRCh38.tiling",
		snps = "regions/comparison/{acronym}.GRCh37-GRCh38.snps",
		dotplot = "regions/comparison/{acronym}.GRCh37-GRCh38.png"
	params:
		stub = "regions/comparison/{acronym}.GRCh37-GRCh38"
	shell:
		"""
			{tools[mummer]}/nucmer -o -p {params.stub} {input.b37} {input.b38}
			{tools[mummer]}/show-tiling {output.delta} > {output.tiling}
			{tools[mummer]}/show-snps {output.delta} > {output.snps}
			{tools[mummer]}/mummerplot {output.delta} -p {params.stub} -t png -l -S
		"""

rule selfmapCompare:
	input:
		b37 = "regions/sequence/{acronym}.GRCh37.fa",
		b38 = "regions/sequence/{acronym}.GRCh38.fa"
	output:
		selfmap100 = "regions/comparison/{acronym}.GRCh37-GRCh38.k=100.selfmap"
	shell:
		"""
			{tools[selfmap]} -sequence GRCh37={input.b37} GRCh38={input.b38} -o {output.selfmap100} -kmer-size 100
		"""

rule selfmapPlot:
	input:
		b37 = "regions/sequence/{acronym}.GRCh37.fa",
		b38 = "regions/sequence/{acronym}.GRCh38.fa"
	output:
		dotplot = "regions/images/{acronym}.GRCh37-GRCh38.k=100.dotplot.pdf"
	params:
		chromosome_b37 = lambda x: getRegionDefinition( x.acronym, 'GRCh37' )['chromosome'],
		chromosome_b38 = lambda x: getRegionDefinition( x.acronym, 'GRCh37' )['chromosome']
	shell:
		"""
			{tools[Rscript]} --vanilla {tools[selfmap.R]} \
			--s1 'GRCh37={input.b37}' \
			--s2 'GRCh38={input.b38}' \
			--k 100 \
			--chromosome1 {params.chromosome_b37} \
			--chromosome2 {params.chromosome_b38} \
			--output {output.dotplot}
		"""

rule vsNanoporeSelfmapPlot:
	input:
		build = "regions/sequence/{acronym}.{build}.fa",
		assembled = "%s/projects/nanopore/{acronym}_assembly_test/{method}/consensus.fasta" % ROOT
	output:
		dotplot = "regions/images/{acronym}.{build}-{method}.k=100.dotplot.pdf"
	params:
		chromosome = lambda x: getRegionDefinition( x.acronym, x.build )['chromosome'],
	shell:
		"""
			{tools[Rscript]} --vanilla {tools[selfmap.R]} \
			--s1 '{wildcards.build}={input.build}' \
			--s2 '{wildcards.method}={input.assembled}' \
			--k 100 \
			--chromosome1 {params.chromosome} \
			--chromosome2 {params.chromosome} \
			--output {output.dotplot}
		"""

rule vsNanoporeNucmerCompare:
	input:
		reference = "regions/sequence/{acronym}.{build}.fa",
		assembled = "%s/shared/analysis/hroberts/nanopore/{acronym}_assembly_test/{method}/consensus.fasta" % ROOT
	output:
		delta = "regions/comparison/{acronym}.{build}-{method}.delta",
		coords = "regions/comparison/{acronym}.{build}-{method}.coords",
		tiling = "regions/comparison/{acronym}.{build}-{method}.tiling",
		snps = "regions/comparison/{acronym}.{build}-{method}.snps",
		dotplot = "regions/comparison/{acronym}.{build}-{method}.png"
	params:
		stub = "regions/comparison/{acronym}.{build}-{method}"
	shell:
		"""
			{tools[mummer]}/nucmer -o -p {params.stub} {input.reference} {input.assembled}
			{tools[mummer]}/show-tiling {output.delta} > {output.tiling}
			{tools[mummer]}/show-snps {output.delta} > {output.snps}
			{tools[mummer]}/mummerplot {output.delta} -p {params.stub} -t png -l -S
		"""
