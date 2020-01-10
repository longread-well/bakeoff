import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_file = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{assembly_type}/{tech}/{build}/{acronym}/{method}.bam"

def get_input_fasta(assembly_type, tech, build, acronym, method):
	fasta = ROOT + "/analysis/{assembly_type}/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta"
	return fasta.format(assembly_type = assembly_type, tech = tech, build = build, acronym = acronym, method = method)

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	return result[0]

rule Align_WGA_contigs:
	input:
		fasta = ROOT + "/analysis/whole_genome_assembly/data/{tech}/{method}/Output.symlink.fasta",
		ref = lambda wildcards: REF[wildcards.build]
	output:
		bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/whole_genome_assembly/{tech}/{build}/whole_genome/{method}.bam"
	params:
		minimap2 = tools['minimap2'],
		samtools = tools['samtools']
	shell:
		"""
		{params.minimap2} -ax asm5 {input.ref} {input.fasta} --secondary=no | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
		"""

rule Subset_WGA_alignment:
	input:
		bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/whole_genome_assembly/{tech}/{build}/whole_genome/{method}.bam"
	output:
		bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/whole_genome_assembly/{tech}/{build}/{acronym}/{method}.bam"
	params:
		samtools = tools['samtools'],
		chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
		start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
		end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
	wildcard_constraints:
		acronym = ".{3,4}"
	shell:
		"""
		{params.samtools} view -b {input.bam} {params.chromosome}:{params.start}-{params.end} -o {output.bam}
		"""

rule Align_regional_contigs:
	input:
		fasta = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
		ref = lambda wildcards: REF[wildcards.build]
	output:
		bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/regional_assembly/{tech}/{build}/{acronym}/{method}.bam"
	params:
		minimap2 = tools['minimap2'],
		samtools = tools['samtools']
	wildcard_constraints:
		acronym = ".{3,4}"
	shell:
		"""
		{params.minimap2} -ax asm5 {input.ref} {input.fasta} --secondary=no | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
		"""

output_files = []
assembly_type = "regional_assembly"
for tech in ["ONT", "PB-CCS", "PB-CLR"]:
	for build in ["GRCh38"]:
		for acronym in ['IGH']:
			if tech == "ONT":
				methods = ["Flye", "Flye_Medaka", "Canu", "Canu_Medaka", "Wtdbg2", "Wtdbg2_Medaka", "Canu_Purge"]
			elif tech == "PB-CCS" or tech == "PB-CLR":
				methods = ["Flye", "Canu", "Wtdbg2", "Canu_Purge"]
			for method in methods:
				output_files.append(output_file.format(assembly_type = assembly_type, tech = tech, build = build, acronym = acronym, method = method))
assembly_type = "whole_genome_assembly"
for tech in ["ONT", "PB-CCS", "PB-CLR", "10X"]:
	for build in ["GRCh38"]:
		for acronym in acronyms:
			break
			if tech == "ONT":
				methods = ["Wtdbg2"]
			elif tech == "PB-CCS":
				methods = ["Canu", "Wtdbg2"]
			elif tech == "PB-CLR":
				methods = ['Wtdbg2']
			elif tech == "10X":
				methods = ["Supernova"]
			for method in methods:
				output_files.append(output_file.format(assembly_type = assembly_type, tech = tech, build = build, acronym = acronym, method = method))

rule All:
	input:
		output_files

localrules: All, Align_regional_contigs, Subset_WGA_alignment

# To run:
# snakemake -s pipelines/Alignment.snakefile All --cluster "qsub -P todd.prjc -q himem.qh -pe shmem 4 -N alignment -cwd -j y -o qsub_output.log -e qsub_error.log" --jobs 2
