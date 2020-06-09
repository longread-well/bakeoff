include: "initialize.snakefile"

canu_correct = Operation("canu_correct", ["tech"])
canu_trim = Operation("canu_trim", ['tech'])

rule canu_correct:
	input:
		reads = join("{prefix}", "data", canu_correct.wildcards.tech, "{region}", "raw_reads.fastq"),
	output:
		reads = "{prefix}" / canu_hybrid("shared") / "data" / canu_correct(...) / "{region}" / "canu.correctedReads.fasta.gz",
		symlink = "{prefix}" / canu_hybrid("shared") / "data" / canu_correct(...) / "{region}" / "reads.symlink.fasta.gz",
	params:
		output_folder = "{prefix}" / canu_hybrid("shared") / "data" / canu_correct(...) / "{region}",
		tech_flag = lambda wildcards: {"ONT": "nanopore", "CLR": "pacbio"}[wildcards.canu_correct_tech],
		length = lambda wildcards: get_region_length(wildcards.region),
	threads: 4
	shell:
		"""
		rm -rf {params.output_folder}
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{config[tools][canu20]} -correct -d {params.output_folder} -p canu genomeSize={params.length} 'useGrid=false' -{params.tech_flag} {input.reads}
		ln -rs {output.reads} {output.symlink}
		"""

rule canu_trim:
	input:
		reads = "{prefix}" / canu_hybrid("shared") / "data" / canu_correct("{tech}") / "{region}" / "reads.symlink.fasta.gz",
	output:
		reads = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim("{tech}") / "{region}" / "canu.trimmedReads.fasta.gz",
		symlink = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim("{tech}") / "{region}" / "reads.symlink.fasta.gz",
	params:
		output_folder = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim("{tech}") / "{region}",
		tech_flag = lambda wildcards: {"ONT": "nanopore", "CLR": "pacbio"}[wildcards.tech],
		length = lambda wildcards: get_region_length(wildcards.region),
	threads: 4
	wildcard_constraints:
		tech = "CLR|ONT"
	shell:
		# -stopOnLowCoverage=1 was included to avoid Canu from raising an error
		# useGrid is disabled for regional assembly
		"""
		rm -rf {params.output_folder}
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{config[tools][canu20]} -trim -d {params.output_folder} -p canu genomeSize={params.length} 'useGrid=false' -corrected -{params.tech_flag} {input.reads}
		ln -rs {output.reads} {output.symlink}
		"""

rule canu_trim_CCS:
	input:
		reads = join("{prefix}", "data", canu_trim.wildcards.tech, "{region}", "raw_reads.fastq"),
	output:
		reads = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim(tech = ...) / "{region}" / "canu.trimmedReads.fasta.gz",
		symlink = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim(tech = ...) / "{region}" / "reads.symlink.fasta.gz",
	params:
		output_folder = "{prefix}" / canu_hybrid("shared") / "data" / canu_trim(tech = ...) / "{region}",
		length = lambda wildcards: get_region_length(wildcards.region),
	threads: 4
	wildcard_constraints:
		canu_trim_tech = "CCS"
	shell:
		"""
		rm -rf {params.output_folder}
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{config[tools][canu20]} -trim -d {params.output_folder} -p canu genomeSize={params.length} 'useGrid=false' -pacbio-hifi {input.reads} -stopOnLowCoverage=1
		ln -rs {output.reads} {output.symlink}
		"""

def get_canu_assemble_input(wildcards):
	input_files = dict()
	techs = ['ONT', 'CCS', 'CLR']
	region = str(wildcards.region)
	prefix = str(wildcards.prefix)
	for tech in techs:
		if tech in wildcards.canu_hybrid_tech.split('+'):
			reads = prefix / canu_hybrid("shared") / "data" / canu_trim(tech) / region / "reads.symlink.fasta.gz"
			input_files[tech] = reads
	return input_files

def get_canu_assemble_parameters(wildcards, input):
	parameters = ""
	techs = wildcards.canu_hybrid_tech.split('+')
	if 'ONT' in techs:
		parameters += " -trimmed -corrected -nanopore " + input.ONT
	if 'CCS' in techs:
		parameters += " -trimmed -corrected -pacbio " + input.CCS
	if 'CLR' in techs:
		parameters += " -trimmed -corrected -pacbio " + input.CLR
	return parameters



rule canu_assemble:
	input:
		unpack(get_canu_assemble_input)
	output:
		contigs = "{prefix}" / canu_hybrid(...) / 'data' / "{region}" / "canu.contigs.fasta",
		symlink = "{prefix}" / canu_hybrid(...) / 'data' / "{region}" / "assembly.symlink.fasta",
	params:
		output_folder = "{prefix}" / canu_hybrid(...) / 'data' / "{region}",
		length = lambda wildcards: get_region_length(wildcards.region),
		canu_input = get_canu_assemble_parameters
	threads: 8
	wildcard_constraints:
		canu_hybrid_tech = r"(CCS\+CLR\+ONT|CCS\+ONT|CCS\+CLR|CLR\+ONT)"
	shell:
		# -stopOnLowCoverage=1 was included to avoid Canu from raising an error
		# useGrid is disabled for regional assembly
		"""
		rm -rf {params.output_folder}
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{config[tools][canu20]} -assemble -d {params.output_folder} -p canu genomeSize={params.length} 'useGrid=false' {params.canu_input}
		ln -rs {output.contigs} {output.symlink}
		"""
