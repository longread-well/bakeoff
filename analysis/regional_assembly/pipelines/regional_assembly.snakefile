import csv, os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	return result[0]

def get_output_path(tech, build, acronym, method):
	return ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/".format(tech=tech, build=build, acronym=acronym, method=method)

def get_output_file(tech, build, acronyms, method):
	return get_output_path(tech, build, acronyms, method) + "Output.symlink.fasta"

def get_bam_file(tech, build):
	assert build == "GRCh38"
	if tech == "ONT":
		return "/well/longread/projects/nanopore/data/aligned/minimap2/combined/minimap2_all-build38.bam"
	elif tech == "PB-CLR":
		return "/well/longread/projects/pacbio/data/aligned/pbmm2/clr/CLR_190625/m64016_190625_101024.GRCh38.bam"
	elif tech == "PB-CCS":
		return "/well/longread/projects/pacbio/data/aligned/pbmm2/ccs/combined/GRCh38.allCCSreads.bam"
	else:
		raise Exception("Parameter error")

def get_fastq_file(tech):
	# TODO
	if tech == "ONT":
		return ROOT + "/resources/data/ONT/JK_HV31.ont.fastq"
	elif tech == "PB-CCS":
		return ROOT + "/resources/data/PB-CCS/JK_HV31.ccs.fastq"
	elif tech == "PB-CLR":
		return ROOT + "/resources/data/PB-CLR/JK_HV31.clr.fastq"
	else:
		raise Exception("Parameter error")

wildcard_constraints:
	build='GRCh3[7|8]',
	acronym=r'\w+'

rule Extract_reads:
	input:
		bam = lambda wildcards: get_bam_file(wildcards.tech, wildcards.build),
		fastq = lambda wildcards: get_fastq_file(wildcards.tech),
		region_definition = ROOT + '/resources/regions.tsv'

	output:
		mapped_bam = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/mapped_reads.bam",
		mapped_fastq = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/mapped_reads.fastq",
		raw_fastq = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq",

	params:
		chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
		start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
		end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
		samtools = tools['samtools'],
		extract_raw_reads = tools['extract_raw_reads'],
		bakeoff_env = tools['bakeoff_env']

	shell:
		"""
		{params.samtools} view -b {input.bam} {params.chromosome}:{params.start}-{params.end} -o {output.mapped_bam}
		{params.samtools} fastq {output.mapped_bam} > {output.mapped_fastq}
		set +u; source activate {params.bakeoff_env}; set -u
		python {params.extract_raw_reads} -mapped {output.mapped_fastq} -all {input.fastq} -out {output.raw_fastq}
		"""

rule Flye:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq"

	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Flye/assembly.fasta",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Flye/Output.symlink.fasta"

	params:
		flye = tools['flye'],
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Flye"),
		tech_flag = lambda wildcards: {"ONT": "nano-raw", "PB-CLR": "pacbio-raw", "PB-CCS": "pacbio-corr"}[wildcards.tech],
		length = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["length"],

	shell:
		"""
		{params.flye} --{params.tech_flag} {input.reads} -g {params.length} -o {params.output_folder} -t 4 --asm-coverage 40
		ln -s {output.consensus} {output.symlink}
		"""

rule Medaka:
	input:
		scaffold = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq"
	output:
		assembly = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}_Medaka/consensus.fasta",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}_Medaka/Output.symlink.fasta"
	params:
		medaka = tools['medaka'],
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, wildcards.method + "_Medaka"),

	shell:
		# TODO: confirm medaka -m parameter. See https://nanoporetech.github.io/medaka/installation.html
		"""
		# set +u and set -u are necessary: see https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
		set +u; source {params.medaka}; set -u
		medaka_consensus -i {input.reads} -d {input.scaffold} -o {params.output_folder} -t 4 -m r941_prom_high
		ln -s {output.assembly} {output.symlink}
		"""

rule Canu:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq"

	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Canu/canu.contigs.fasta",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Canu/Output.symlink.fasta"

	params:
		canu = tools['canu'],
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Canu"),
		tech_flag = lambda wildcards: {"ONT": "nanopore-raw", "PB-CLR": "pacbio-raw", "PB-CCS": "pacbio-hifi"}[wildcards.tech],
		length = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["length"],

	shell:
		# -stopOnLowCoverage=1 was included to avoid Canu from raising an error
		# useGrid is disabled for regional assembly
		"""
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{params.canu} -d {params.output_folder} -p canu genomeSize={params.length} 'useGrid=false' -{params.tech_flag} {input.reads} -stopOnLowCoverage=1
		ln -s {output.consensus} {output.symlink}
		"""

rule Wtdbg2:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq"
	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2/Wtdbg2.ctg.fa",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2/Output.symlink.fasta"

	params:
		wtdbg2 = tools['wtdbg2'],
		output_prefix = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Wtdbg2") + "Wtdbg2",
		tech_flag = lambda wildcards: {"ONT": "ont", "PB-CLR": "sq", "PB-CCS": "ccs"}[wildcards.tech],
		length = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["length"],
	shell:
		# TODO: check wtdbg2 -x option (RSII or Squell for CLR?) See https://github.com/ruanjue/wtdbg2
		"""
		{params.wtdbg2}/wtdbg2 -x {params.tech_flag} -g {params.length} -t 4 -i {input.reads} -fo {params.output_prefix}
		{params.wtdbg2}/wtpoa-cns -t 4 -i {params.output_prefix}.ctg.lay.gz -fo {params.output_prefix}.ctg.fa

		ln -s {output.consensus} {output.symlink}
		"""

rule Falcon:
	# TODO
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq"
	output:
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Falcon/Output.symlink.fasta"
	params:
		fofn = ROOT + "/analysis/regional_assembly/resources/falcon_input_fasta.fofn",
		cfg = ROOT + "/analysis/regional_assembly/resources/falcon_{tech}.cfg",
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Falcon"),
	shell:
		"""
		touch {params.fofn}
		echo {input.reads} >> {params.fofn}
		source activate pb-assembly
		rm -rf {params.output_folder}
		mkdir {params.output_folder}
		cd {params.output_folder}
		fc_run {params.cfg}
		rm {params.fofn}
		"""
rule SDA:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq",
		contigs = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
	output:
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}_SDA/Output.symlink.fasta"
	params:
		sda = tools['sda'],
		tech_flag = lambda wildcards: {"ONT": "ont", "PB-CCS": "ccs", "PB-CLR": "subread"}[wildcards.tech],
		output_path = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, wildcards.method + "_SDA"),
		samtools = tools['samtools']
	shell:
		"""
		{params.samtools} faidx {input.contigs}
		echo "{input.reads}" > {params.output_path}reads.fofn
		{params.sda} denovo --fofn {params.output_path}/reads.fofn --species human --ref {input.contigs} --dir {params.output_path} --prefix sda --assemblers canu --drmaa " -l mfree={{resources.mem}}G -pe serial {{threads}} -l h_rt=128:00:00 -V -cwd -j y -S /bin/bash -N SDA -P todd.prjc -q himem.qh -o qsub_output.log -e qsub_error.log"
		"""

output_files = []
for tech in ["ONT", "PB-CCS", "PB-CLR"]:
	for build in ["GRCh38"]:
		for acronym in acronyms:
			if tech == "ONT":
				methods = ["Flye", "Flye_Medaka", "Canu", "Canu_Medaka", "Wtdbg2", "Wtdbg2_Medaka"]
			elif tech == "PB-CCS" or tech == "PB-CLR":
				methods = ["Flye", "Canu", "Wtdbg2"]
			for method in methods:
				output_files.append(get_output_file(tech, build, acronym, method))


rule All:
	input:
		output_files
