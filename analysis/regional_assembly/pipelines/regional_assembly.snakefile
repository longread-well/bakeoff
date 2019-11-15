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
	build_number = build[4:]
	if tech == "ONT":
		return "/well/longread/projects/nanopore/data/minimap2_all-build{build_number}.bam".format(build_number = build_number)
	elif tech == "PB-CLR":
		return "/well/longread/projects/pacbio/data/CLR_190625/m64016_190625_101024.GRCh38.bam"
	elif tech == "PB-CCS":
		return "/well/longread/projects/pacbio/data/GRCh38.allCCSreads.bam"

wildcard_constraints:
	build='GRCh3[7|8]',
	acronym=r'\w+'


rule Extract_target_region:
	input:
		bam = lambda wildcards: get_bam_file(wildcards.tech, wildcards.build),
		region_definition = ROOT + '/resources/regions.tsv'

	output:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"

	params:
		chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
		start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
		end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
		samtools = tools['samtools'],

	shell:
		"{params.samtools} view -b {input.bam} {params.chromosome}:{params.start}-{params.end} | {params.samtools} fastq - > {output.reads}"

rule Flye:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"

	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Flye/assembly.fasta",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Flye/Output.symlink.fasta"

	params:
		flye = tools['flye'],
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Flye"),
		# TODO: confirm tech_flag for PB-CCS
		tech_flag = lambda wildcards: {"ONT": "nano-raw", "PB-CLR": "pacbio-raw", "PB-CCS": "pacbio-corr"}[wildcards.tech]

	shell:
		"""
		{params.flye} --{params.tech_flag} {input.reads} -g 1m -o {params.output_folder} -t 4 --asm-coverage 40
		ln -s {output.consensus} {output.symlink}
		"""

rule Medaka:
	input:
		scaffold = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"
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
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"

	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Canu/canu.contigs.fasta",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Canu/Output.symlink.fasta"

	params:
		canu = tools['canu'],
		output_folder = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Canu"),
		tech_flag = lambda wildcards: {"ONT": "nanopore-raw", "PB-CLR": "pacbio-raw", "PB-CCS": "pacbio-hifi"}[wildcards.tech]

	shell:
		# -stopOnLowCoverage=1 was included to avoid Canu from raising an error
		# useGrid is disabled for regional assembly
		"""
		module load gcc/5.4.0
		module load java/1.8.0_latest
		{params.canu} -d {params.output_folder} -p canu genomeSize=1m 'useGrid=false' -{params.tech_flag} {input.reads} -stopOnLowCoverage=1
		ln -s {output.consensus} {output.symlink}
		"""

rule Wtdbg2:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"
	output:
		consensus = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2/Wtdbg2.ctg.fa",
		symlink = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2/Output.symlink.fasta"

	params:
		wtdbg2 = tools['wtdbg2'],
		output_prefix = lambda wildcards: get_output_path(wildcards.tech, wildcards.build, wildcards.acronym, "Wtdbg2") + "Wtdbg2",
		tech_flag = lambda wildcards: {"ONT": "ont", "PB-CLR": "sq", "PB-CCS": "ccs"}[wildcards.tech]
	shell:
		# TODO: check wtdbg2 -x option (RSII or Squell for CLR?) See https://github.com/ruanjue/wtdbg2
		"""
		{params.wtdbg2}/wtdbg2 -x {params.tech_flag} -g 1m -t 4 -i {input.reads} -fo {params.output_prefix}
		{params.wtdbg2}/wtpoa-cns -t 4 -i {params.output_prefix}.ctg.lay.gz -fo {params.output_prefix}.ctg.fa

		ln -s {output.consensus} {output.symlink}
		"""

rule Falcon:
	input:
		reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Extract_target_region/reads.fastq"
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
	threads: 8
