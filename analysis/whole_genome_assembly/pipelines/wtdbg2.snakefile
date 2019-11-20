tools = {
	'wtdbg2': '/users/kwiatkowski/gav/Projects/Software/3rd_party/wtdbg2/2.5-20190621',
}

data = {
	"ccs": "data/ccs/{sample}.ccs.fastq"
}

flags = {
	"pacbio-ccs": 'ccs',
	"pacbio-clr": 'sq',
	"nanopore": 'ont'
}

samples = [ 'JK_HV31' ]
datatypes = [ 'pacbio-ccs' ] #, 'pacbio-clr', 'nanopore' ]

threads = 24

rule all:
	input:
		[
			'wtdbg2/{sample}.{datatype}.ctg.lay.gz'.format(
				sample = sample,
				datatype = datatype
			)
			for sample in samples
			for datatype in datatypes
		],
		[
			'wtdbg2/{sample}.{datatype}.ctg.fa'.format(
				sample = sample,
				datatype = datatype
			)
			for sample in samples
			for datatype in datatypes
		]

rule make_consensus:
	input:
		layout = 'wtdbg2/{sample}.{datatype}.ctg.lay.gz'
	output:
		contigs = 'wtdbg2/{sample}.{datatype}.ctg.fa'
	params:
		threads = threads
	shell:
		"""
		{tools[wtdbg2]}/wtpoa-cns -t {params.threads} -i {input.layout} -fo {output.contigs}
	"""

rule make_layout:
	input:
		fastq = lambda w: data["ccs"].format( sample = w.sample )
	output:
		layout = 'wtdbg2/{sample}.{datatype}.ctg.lay.gz'
	params:
		preset = lambda w: flags[ w.datatype ],
		threads = threads,
		outputPrefix = 'wtdbg2/{sample}.{datatype}'
	shell:
		"""
		{tools[wtdbg2]}/wtdbg2 -x {params.preset} -t {params.threads} -g 3.2g -i {input.fastq} -fo {params.outputPrefix}
	"""
