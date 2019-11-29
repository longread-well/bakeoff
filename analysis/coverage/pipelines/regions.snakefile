import os
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

tools = {
	'samtools': '/apps/well/samtools/1.8-gcc5.4.0/bin/samtools'
}

data = {
        'JK_HV31': {
                "pacbio-ccs": "%s/projects/pacbio/data/aligned/pbmm2/ccs/combined/GRCh38.allCCSreads.bam" % ROOT,
		"pacbio-clr":  "%s/projects/pacbio/data/aligned/pbmm2/clr/CLR_190625/m64016_190625_101024.GRCh38.bam" % ROOT,	
                "nanopore": "%s/projects/nanopore/data/aligned/minimap2/combined/minimap2_all-build38.bam" % ROOT,
		"10x": "%s/projects/wgs_10x/data/aligned/JK_HV31_3.bam" % ROOT
        }
}

margin = 500000

def loadRegions(
	filename = '%s/shared/analysis/bakeoff-scripts/resources/regions.tsv' % ROOT,
	builds = [ 'GRCh38' ]
):
	import csv
	regions = {}
	with open( filename, newline = '' ) as file:
		regionReader = csv.DictReader( file, delimiter = '\t' )
		for region in regionReader:
			print( region )
			if region['build'] in builds:
				regions[region['acronym']] = dict( region )
	return regions

regions = loadRegions()
print( regions )
rule all:
	input:
		[
			"results/{sample}/{region}/{sample}_{region}_{platform}_coverage_{size}kb.bed".format(
				sample = sample,
				region = region,
				platform = platform,
				size = size
			)
			for sample in data.keys()
			for region in regions.keys()
			for platform in data[sample].keys()
			for size in [ 10, 25, 100 ]
		]


rule computeCoverage:
	input:
		bed = "results/{sample}/{region}/intervals_{size}kb.bed",
		bam = lambda w: data[w.sample][w.platform]
	output:
		bed = "results/{sample}/{region}/{sample}_{region}_{platform}_coverage_{size}kb.bed"
	shell:
		"""
		{tools[samtools]} bedcov -Q20 {input.bed} {input.bam} > {output.bed}
	"""

rule makeBeds:
	output:
		bed = "results/{sample}/{region}/intervals_{size}kb.bed"
	run:
		chromosome = regions[wildcards.region]['chromosome']
		size_in_bp = int( wildcards.size ) * 1000 
		start = int(regions[wildcards.region]['start'])
		end = int(regions[wildcards.region]['end'])
		mid = (start+end)/2.0
		margin = 0.25
		start = mid + (start-mid)*1.25
		end = mid + (end-mid)*1.25
		start = int( round( start / size_in_bp ) * size_in_bp )
		end = int( round( end / size_in_bp ) * size_in_bp )
		print( wildcards.region, chromosome, start, end )
		output = open( output.bed, 'w' )
		for i in range( start, end, size_in_bp ):
			output.write(
				"%s\t%s\t%s\n" % ( chromosome, i, i + size_in_bp )
			)
		output.close()

