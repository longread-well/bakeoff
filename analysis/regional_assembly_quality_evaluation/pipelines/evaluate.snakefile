import csv, os
from Bio import SeqIO
from glob import glob
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

tools = {
        "jellyfish": "jellyfish",
	"Rscript": "/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript",
	"RepeatMasker": "/apps/well/RepeatMasker/open-4-0-9-p2/RepeatMasker"
}

regions = [	
	"GYP", "HLA", "IGH", "IGK", "IGL", "KIR", "TRA", "TRB", "TRG"
]

# We make a version of the validation plot
# for different sets of short read data
validationPlatforms = [
	'PB-CCS+10X+MGI',
	'PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS'
]

# Expected diploid (kmer) coverage is hard-coded for now
# (I got these values from the Jellyfish kmer distribution plots)
expectedDiploidCoverage = {
	'PB-CCS+10X+MGI': {
		'k=22.bq=20': 71.0,
		'k=31.bq=20': 61.0
	},
	'PB-CCS+10X+MGI+illumina_pcr_free+stLFR+CoolMPS': {
		'k=22.bq=20': 177.0,
		'k=31.bq=20': 152.0
	}
}

# Work with these kmer sizes.  Relevant Jellyfish .jf files must already have been
# computed in the path given below
ks = [ '22', '31' ]

# Locations of data files.
data = {
	"whole genome assembly": "%s/users/akl399/bakeoff/resources/data/whole_genome_assembly/CCS/canu19/canu.contigs.fasta" % ROOT,
	"unscaffolded contigs": "%s/users/akl399/bakeoff/analysis/@bam_extract/@wga_contigs/CCS/canu19/data/{region}/contigs.fasta" % ROOT,
	"regional assemblies": "%s/users/akl399/bakeoff/priority.symlink/data/{region}/assembly.symlink.fasta" % ROOT,
	"kmer counts": "%s/shared/analysis/kmer/distribution/{platform}/{platform}.k={k}.bq={bq}.jf" % ROOT
}


rule all:
	input:
		"results/assembly/stitched.fasta",
		jf = [ "results/assembly/stitched.k={k}.jf".format( k = k ) for k in [ 22, 31 ] ],
		multiplicities = [ "results/{region}/{region}.{platform}.k={k}.bq={bq}.txt".format( region = r, platform = p, k = k, bq = '20' ) for r in regions for p in validationPlatforms for k in ks ],
		plots = [ "results/images/{region}.{platform}.k={k}.bq={bq}.png".format( region = r, platform = p, k = k, bq = '20' ) for r in regions for p in validationPlatforms for k in ks ],
		repeatmasker = [ "results/{region}/repeatmasker/assembly.symlink.fasta.out".format( region = r ) for r in regions ]


rule create_stitched_assembly:
	# This rule takes the whole-genome CCS contig
	# removes contigs that went into scaffolding
	# then inserts the scaffolded regional assemblies in their place.

	input:
		wga = data[ "whole genome assembly" ],
		unscaffolded = [ data["unscaffolded contigs"].format( region = region ) for region in regions ],
		finalised = [ data["regional assemblies"].format( region = region ) for region in regions ]
	output:
		"results/assembly/stitched.fasta"
	run:
		contigIDs = set()
		for file in input.unscaffolded:
			for record in SeqIO.parse( file, "fasta" ):
				contigIDs.add( record.id )

		finalisedContigs = []
		for ( regionIndex, file ) in enumerate( input.finalised ):
			print( file )
			for ( contigIndex, record ) in enumerate( SeqIO.parse( file, "fasta" )):
				record.id = "%s_%d" % ( regions[regionIndex], contigIndex )
				finalisedContigs.append( record )

		with open( output[0], 'w' ) as handle:
			SeqIO.write( finalisedContigs, handle, "fasta" )
			sequences = SeqIO.parse( input.wga, "fasta" )
			# the following is a generator, i.e. lazy view of the above
			filtered_seq_iterator = (record for record in sequences if not record.id in contigIDs)
			SeqIO.write( filtered_seq_iterator, handle, "fasta" )

rule count_stitched_assembly_kmers:
	input:
		assembly = "results/assembly/stitched.fasta"
	output:
		jf = "results/assembly/stitched.k={k}.jf"
	params:
		threads = 4
	shell: """
		jellyfish count -t {params.threads} -m {wildcards.k} -s 3G -C -o {output.jf} {input.assembly}
	"""

rule count_regional_assembly_kmers:
	input:
		finalised = data["regional assemblies"]
	output:
		jf = "results/{region}/{region}.assembly.k={k}.jf"
	shell: """
		jellyfish count -m {wildcards.k} -s 100M -C -o {output.jf} {input.finalised}
	"""

rule compute_whole_genome_assembly_kmer_multiplicity:
	# Output multiplicity of each kmer appearing in the regional assembly
	# as computed in the whole-genome stitched assembly
	input:
		assembly = lambda w: data["regional assemblies"].format( region = w.region ),
		assembly_jf = "results/assembly/stitched.k={k}.jf"
	output:
		assembly = "results/{region}/{region}.stitched_assembly.k={k}.txt"
	shell: """
		jellyfish query {input.assembly_jf} -s {input.assembly} -o {output.assembly}
	"""

rule compute_regional_kmer_multiplicity:
	# Output multiplicity of each kmer appearing in the regional assembly
	# as computed in the regional assembly itself
	input:
		assembly = lambda w: data["regional assemblies"].format( region = w.region ),
		assembly_jf = lambda w: "results/{region}/{region}.assembly.k={k}.jf".format( region = w.region, k = w.k )
	output:
		assembly = "results/{region}/{region}.assembly.k={k}.txt"
	shell: """
		jellyfish query {input.assembly_jf} -s {input.assembly} -o {output.assembly}
	"""
		
rule count_platform_kmer_multiplicity:
	# Output multiplicity of each kmer appearing in the regional assembly
	# as computed in the reads from the given platform (e.g. short-read data)
	input:
		assembly = lambda w: data["regional assemblies"].format( region = w.region ),
		platform_jf = lambda w: data["kmer counts"].format( platform = w.platform, k = w.k, bq = w.bq )
	output:
		platform = "results/{region}/{region}.{platform}.k={k}.bq={bq}.txt"
	shell: """
		jellyfish query {input.platform_jf} -s {input.assembly} -o {output.platform}
	"""

rule mask_repeats:
	input:
		assembly = data["regional assemblies"]
	output:
		"results/{region}/repeatmasker/assembly.symlink.fasta.out"
	shell: """
		{tools[RepeatMasker]}  -dir {output} {input.assembly}
	"""

rule plot:
	input:
		regional = "results/{region}/{region}.assembly.k={k}.txt",
		whole_genome = "results/{region}/{region}.stitched_assembly.k={k}.txt",
		validation = "results/{region}/{region}.{platform}.k={k}.bq={bq}.txt",
		repeats = "results/{region}/repeatmasker/assembly.symlink.fasta.out"
	output:
		"results/images/{region}.{platform}.k={k}.bq={bq}.png"
	params:
		script = srcdir( "scripts/plot_kmers_across_region.R" ),
		expected = lambda w: expectedDiploidCoverage[ w.platform ][ 'k=%s.bq=%s' % ( w.k, w.bq ) ] / 2
	shell: """
		{tools[Rscript]} --vanilla {params.script} \
		--name {wildcards.region} \
		--k {wildcards.k} \
		--expected_haploid {params.expected} \
		--regional {input.regional} \
		--whole_genome {input.whole_genome} \
		--validation {input.validation} \
		--repeats {input.repeats} \
		--output {output}
	"""

