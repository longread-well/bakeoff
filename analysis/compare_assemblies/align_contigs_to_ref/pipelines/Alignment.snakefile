import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_file = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/{method}.bam"

rule Align_contigs_to_ref:
	input:
		fasta = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
		ref = lambda wildcards: REF[wildcards.build]
	output:
		bam = output_file
	params:
		minimap2 = tools['minimap2'],
		samtools = tools['samtools']
	shell:
		"""
		{params.minimap2} -a {input.ref} {input.fasta} -g 500000 | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
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
				output_files.append(output_file.format(tech = tech, build = build, acronym = acronym, method = method))

rule All:
	input:
		output_files
