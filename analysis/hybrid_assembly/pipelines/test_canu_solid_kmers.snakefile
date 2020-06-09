import os, pprint, re
ROOT = os.environ.get( 'BAKEOFF_ROOT' )

configfile: "%s/users/akl399/bakeoff/config/master.config.yaml" % ROOT

def readRegions( filename = "%s/shared/analysis/bakeoff-scripts/resources/regions.tsv" % ROOT ):
	import pandas
	X = pandas.read_csv( filename, sep='\t' )
	X = X[ X.build == 'GRCh38' ]
	result = {}
	for index, row in X.iterrows():
		result[ row['acronym'] ] = {
			"acronym": row['acronym'],
			"build": row['build'],
			"chromosome": row['chromosome'],
			"start": row['start'],
			"end": row['end']
		}
	return result

dataSources = {
	"ONT": { "bam": "%s/projects/nanopore/data/aligned/minimap2/combined/minimap2_all-build38.bam" % ROOT },
	"CCS": { "bam": "%s/projects/pacbio/data/aligned/pbmm2/ccs/combined/GRCh38.allCCSreads.bam" % ROOT },
	"CLR": { "bam": "%s/projects/pacbio/data/aligned/pbmm2/clr/CLR_190625/m64016_190625_101024.GRCh38.bam" % ROOT }
}

tools = {
	"pbsim": "/users/kwiatkowski/gav/Projects/Software/3rd_party/PBSIM-PacBio-Simulator/src/pbsim"
}

def get_region_length(acronym):
	region = regions[acronym]
	length = region['end'] - region['start']
	return length

data = {
	"CLR": "%s/projects/pacbio/data/raw/clr/combined/JK_HV31.clr.fastq" % ROOT,
	"ONT": "%s/projects/nanopore/data/raw/combined/JK_HV31.ont.fastq" % ROOT
}
regions = readRegions()

pp = pprint.PrettyPrinter(width=120, compact=True)
pp.pprint( regions )


kmers = {
	"CLR": "%s/shared/analysis/kmer/distribution/PB-CCS+10X+MGI/PB-CCS+10X+MGI.k=16.bq=20.solid_kmer_superset_mhap.txt" % ROOT,
	"ONT": "%s/shared/analysis/kmer/distribution/PB-CCS+10X+MGI/PB-CCS+10X+MGI.k=16.bq=20.solid_kmer_superset_mhap.txt" % ROOT,
	"pbsim-CLR": "%s/shared/analysis/assembly/regions/data/reference/GRCh38_unique_16mers.mhap.txt" % ROOT,
	"pbsim_sourcereads-CLR": "%s/shared/analysis/assembly/regions/data/reference/GRCh38_unique_16mers.mhap.txt" % ROOT
}

print( kmers )

tools = {
	"canu": "/users/kwiatkowski/gav/Projects/Software/canu/Linux-amd64/bin/canu"
}

rule all:
	input:
		[
			"data/reads/{platform}/{region}/{platform}_{region}.fastq.gz".format(
				platform = platform, region = region
			)
			for region in regions.keys()
			for platform in dataSources.keys()
		],
		[ "data/reference/GRCh38_{region}.fa".format( region = region ) for region in regions.keys() ],
		[
			"data/reads/pbsim-CLR/{region}/pbsim-CLR_{region}.source_reads.gz".format(
				region = region
			) for region in regions.keys()
		],
		[
			"data/reads/pbsim_sourcereads-CLR/{region}/pbsim_sourcereads-CLR_{region}.overlaps.txt.gz".format(
				region = region
			)
			for region in regions.keys()
		],
		[
			"results/canu/correction/{method}/{platform}/{region}/canu.correctedReads.fasta.gz".format(
				method = method,
				platform = platform,
				region = region
			)
			for method in [ 'standard', 'solid' ]
			for region in regions.keys()
			for platform in list( data.keys() ) + [ 'pbsim-CLR', 'pbsim_sourcereads-CLR' ]
		],
		[
			"results/canu/correction/{method}/{platform}/{region}/canu.correctedReads.realigned.sam.gz".format(
				method = method,
				platform = platform,
				region = region
			)
			for method in [ 'standard', 'solid' ]
			for region in regions.keys()
			for platform in  [ 'pbsim-CLR', 'pbsim_sourcereads-CLR' ]
		]

def getSolidKmerOption( wildcards ):
	if wildcards.method == 'standard':
		return ''
	elif wildcards.method == 'solid':
		return 'corMhapSolidKmers=%s' % kmers[wildcards.platform]

rule realign:
	input:
		corrected = "results/canu/correction/{method}/{platform}/{region}/canu.correctedReads.fasta.gz",
		ref = "data/reads/pbsim-CLR/{region}/pbsim-CLR_{region}.source_reads.gz"
	output:
		aligned = "results/canu/correction/{method}/{platform}/{region}/canu.correctedReads.realigned.sam.gz"
	shell: """
		minimap2 -x map-pb -a --cs=long {input.ref} {input.corrected} | gzip -c > {output.aligned}
	"""

rule correct:
	input:
                reads = "data/reads/{platform}/{region}/{platform}_{region}.fastq.gz",
		kmers = lambda w: kmers[w.platform]
	output:
                reads = "results/canu/correction/{method}/{platform}/{region}/canu.correctedReads.fasta.gz"
	params:
		output_folder = "results/canu/correction/{method}/{platform}/{region}",
		tech_flag = lambda wildcards: {"ONT": "nanopore", "CLR": "pacbio", "pbsim-CLR": "pacbio", "pbsim_sourcereads-CLR": "pacbio" }[wildcards.platform],
		length = lambda wildcards: get_region_length(wildcards.region),
		kmerOption = getSolidKmerOption
		
	threads: 4
	shell: """
		rm -rf {params.output_folder}
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{tools[canu]} \
		-correct \
		-d {params.output_folder} \
		-p canu \
		genomeSize={params.length} \
		useGrid=false \
		maxThreads=4 \
		{params.kmerOption} \
		-{params.tech_flag} \
		{input.reads}
	"""

rule extractReads:
	input:
		reads = lambda w: dataSources[w.platform]['bam']
	output:
		"data/reads/{platform}/{region}/{platform}_{region}.fastq.gz"
	wildcard_constraints:
		platform = "[A-Z]+"
	params:
		chromosome = lambda w: regions[w.region]['chromosome'],
		start = lambda w: regions[w.region]['start'],
		end = lambda w: regions[w.region]['end']
	shell: """
		samtools view -b {input.reads} {params.chromosome}:{params.start}-{params.end} \
		| samtools bam2fq - | bgzip -c > {output}.tmp
		mv {output}.tmp {output}
	"""

rule extractRef:
	input:
		reference = "%s/projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa" % ROOT
	output:
		"data/reference/GRCh38_{region}.fa"
	params:
		chromosome = lambda w: regions[w.region]['chromosome'],
		start = lambda w: regions[w.region]['start'],
		end = lambda w: regions[w.region]['end']
	shell: """
		samtools faidx {input.reference} {params.chromosome}:{params.start}-{params.end} > {output}.tmp
		mv {output}.tmp {output}
		samtools faidx {output}
	"""

rule createSourceSimulationReadSets:
	input:
		"data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.source_reads.gz"
	output:
		"data/reads/pbsim_sourcereads-{platform}/{region}/pbsim_sourcereads-{platform}_{region}.fastq.gz"
	wildcard_constraints:
		platform = "[A-Z]+"
	shell: """
		~/Projects/Software/3rd_party/seqtk/seqtk seq -F '?' {input} | gzip -c > {output}
	"""

rule indexSourceSimulationReads:
	input:
		"data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.source_reads.gz"
	output:
		"data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.source_reads.mmi"
	shell: """
		minimap2 -x map-pp -d {output} {input}
	"""

rule computeActualOverlaps:
	input:
		maf = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.maf.gz",
		fasta = "data/reference/GRCh38_{region}.fa"
	output:
		"data/reads/pbsim_sourcereads-{platform}/{region}/pbsim_sourcereads-{platform}_{region}.overlaps.txt.gz"
	params:
		tmp = lambda w, output: [ elt + "._under_construction" for elt in output ]
	run:
		whitespace = re.compile(r"\s+")
		header = open( input.fasta, 'r' ).readline()
		header = header.strip( "> \t\n" )
		shell( "rm -f {params.tmp}" )
		count = 1
		regions = []
		print( "Loading reads...\n" )
		for line in shell( "zcat {input.maf} | cut -c1-500 | grep '^s' | grep '%s' | cut -d' ' -f1-" % header, iterable = True ):
			line = line.strip()
			line = whitespace.sub( " ", line )
			elts = line.split()
			chromosome = header
			start = int( elts[2] )
			length = int( elts[3] )
			end = start + length
			regions.append( ( "S1_%d" % count, start, end, length ) )
			count = count + 1
		print( "Computing overlaps...\n" )
		result = open( params.tmp[0], 'w' )
		result.write( "r1_id\tr1_start\tr1_end\tr2_id\tr2_start\tr2_end\tlength_of_overlap\n" )
		for i in range( 0, len( regions )):
			count = 0
			r1 = regions[i]
			for j in range( 0, len( regions )):
				r2 = regions[j]
				overlap = min( max( r1[2] - r2[1], 0 ), max( r2[2] - r1[1], 0 ))
				if overlap > 0:
					result.write(
						"%s\t%d\t%d\t%s\t%d\t%d\t%d\n" % (
							r1[0], r1[1], r1[2],
							r2[0], r2[1], r2[2],
							overlap
						)
					)
					count = count + 1
			print( "Read %s has %d overlapping other reads.\n" % ( r1[0], count - 1 ))
		result.close()
		shell( """gzip -c {params.tmp} > {output}; rm {params.tmp}""" )

rule extractSourceSimulationReads:
	input:
		maf = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.maf.gz",
		fasta = "data/reference/GRCh38_{region}.fa"
	output:
		"data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.source_reads.gz"
	params:
		tmp = lambda w, output: [ elt + "._under_construction" for elt in output ]
	run:
		whitespace = re.compile(r"\s+")
		header = open( input.fasta, 'r' ).readline()
		header = header.strip( "> \t\n" )
		shell( "rm -f {params.tmp}" )
		count = 1
		for line in shell( "zcat {input.maf} | cut -c1-500 | grep '^s' | grep '%s' | cut -d' ' -f1-" % header, iterable = True ):
			line = line.strip()
			line = whitespace.sub( " ", line )
			elts = line.split()
			chromosome = header
			print( line )
			print( elts )
			start = int( elts[2] )
			length = int( elts[3] )
			end = start + length
			cmd = "samtools faidx {input.fasta} %s:%d-%d  | sed -e 's/>/>%s /' >> {params.tmp}" % ( chromosome, start, end, "S1_%d" % count )
			print( cmd )
			shell( cmd )
			count = count + 1
		shell( "mv {params.tmp} {output}" )
		shell( "gzip {output}" )
	
rule simulateCLR:
	input:
		reads = "data/reads/{platform}/{region}/{platform}_{region}.fastq.gz",
		fasta = "data/reference/GRCh38_{region}.fa"
	output:
		reads = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.fastq.gz",
		ref = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.ref.gz",
		maf = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}.maf.gz"
	params:
		prefix = "data/reads/pbsim-{platform}/{region}/pbsim-{platform}_{region}",
		pbsim = "/users/kwiatkowski/gav/Projects/Software/3rd_party/PBSIM-PacBio-Simulator/src/pbsim",
		model = "/users/kwiatkowski/gav/Projects/Software/3rd_party/PBSIM-PacBio-Simulator/data/model_qc_clr"
	shell: """
		#{params.pbsim} --data-type {wildcards.platform} --depth 20 --sample-fastq {input.reads} {input.fasta} --prefix {params.prefix}
		{params.pbsim} --data-type {wildcards.platform} --depth 20 --model_qc {params.model} {input.fasta} --prefix {params.prefix}
		mv {params.prefix}_0001.fastq {params.prefix}.fastq
		mv {params.prefix}_0001.ref {params.prefix}.ref
		mv {params.prefix}_0001.maf {params.prefix}.maf
		gzip {params.prefix}.ref
		gzip {params.prefix}.maf
		bgzip {params.prefix}.fastq
	"""

rule makeMhapKmerFile:
	input:
		txt = "data/reference/GRCh38_unique_16mers.txt.gz"
	output:
		#mhap = "data/reference/GRCh38_unique_16mers.mhap.txt"
		mhap = kmers['pbsim-CLR']
	shell: """
		../../kmer/pipelines.symlink/scripts/make_mhap_kmer_file {input.txt} > {output.mhap}
	"""

rule findGRCh38UniqueKmers:
	input:
		reference = "%s/projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa" % ROOT
	output:
		jf = "data/reference/GRCh38_16mers.jf",
		txt = "data/reference/GRCh38_unique_16mers.txt.gz"
	params:
		output = "data/reference/GRCh38_unique_16mers.txt"
	shell: """
		jellyfish count -m 16 -s 100M {input.reference} -C -o {output.jf}
		jellyfish dump -c --lower-count 1 --upper-count 1 -o {params.output} {output.jf}
		gzip {params.output}
	"""

