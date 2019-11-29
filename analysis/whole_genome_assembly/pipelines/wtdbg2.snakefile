import os
ROOT = os.environ[ 'BAKEOFF_ROOT' ]

tools = {
	'wtdbg2': '/users/kwiatkowski/gav/Projects/Software/3rd_party/wtdbg2/2.5-20190621',
}

data = {
	'JK_HV31': {
		"pacbio-ccs": [
			"%s/projects/pacbio/data/raw/ccs/data/ccs/HiFi_190625/m64016_190625_230427.ccs.fastq.gz" % ROOT,
			"%s/projects/pacbio/data/raw/ccs/data/ccs/HiFi_190628/m64016_190628_125725.ccs.fastq.gz" % ROOT,
			"%s/projects/pacbio/data/raw/ccs/data/ccs/HiFi_190701/m64016_190701_082735.ccs.fastq.gz" % ROOT
		],
		"pacbio-clr": [ "%s/projects/pacbio/data/raw/clr/CLR_190625/m64016_190625_101024.subreads.fastq.gz" % ROOT ],
		"nanopore": [
			"%s/projects/nanopore/data/raw/WTON000392/trimmed/fastq-trimmed-all.fastq.gz"% ROOT,
			"%s/projects/nanopore/data/raw/WTON000393/trimmed/fastq-trimmed-all.fastq.gz"% ROOT
		]
	}
}

flags = {
	"pacbio-ccs": 'ccs',
	"pacbio-clr": 'sq',
	"nanopore": 'ont'
}

samples = [ 'JK_HV31' ]
datatypes = [ 'pacbio-ccs', 'pacbio-clr', 'nanopore' ]

threads = 40

rule all:
	input:
		[
			'wtdbg2/{datatype}/{sample}.{datatype}.ctg.lay.gz'.format(
				sample = sample,
				datatype = datatype
			)
			for sample in samples
			for datatype in datatypes
		],
		[
			'wtdbg2/{datatype}/{sample}.{datatype}.ctg.fa'.format(
				sample = sample,
				datatype = datatype
			)
			for sample in samples
			for datatype in datatypes
		]

rule make_consensus:
	input:
		layout = 'wtdbg2/{datatype}/{sample}.{datatype}.ctg.lay.gz'
	output:
		contigs = 'wtdbg2/{datatype}/{sample}.{datatype}.ctg.fa'
	params:
		threads = threads
	shell:
		"""
		{tools[wtdbg2]}/wtpoa-cns -t {params.threads} -i {input.layout} -fo {output.contigs}
	"""

rule make_layout:
	input:
		fastq = lambda w: data[w.sample][w.datatype]
	output:
		layout = 'wtdbg2/{datatype}/{sample}.{datatype}.ctg.lay.gz'
	params:
		preset = lambda w: flags[ w.datatype ],
		threads = threads,
		outputPrefix = 'wtdbg2/{datatype}/{sample}.{datatype}'
	shell:
		"""
		{tools[wtdbg2]}/wtdbg2 -x {params.preset} -t {params.threads} -g 3.2g -fo {params.outputPrefix} -i {input.fastq}
	"""
