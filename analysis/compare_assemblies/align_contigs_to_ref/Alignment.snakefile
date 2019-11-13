import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_input_fasta(method):
	prefix = ROOT + "/analysis/de_novo_assembly/data/{tech}/{build}/{acronym}/"

	if method in ("flye-racon-medaka", "canu-racon-medaka", "wtdbg2-racon-medaka"):
		filename = method + "/Output.symlink.fasta"
	elif method in ("flye-racon", "canu-racon", "wtdbg2-racon"):
		filename = method + "-medaka" + "/3xRacon/racon-3.ctg.fa.gz"
	elif method in ("flye", "canu", "wtdbg2"):
		filename = method + "-racon-medaka" + "/Racon_input.symlink.fasta"

	return prefix + filename

output_file = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/{method}.bam"

rule Align_contigs_to_ref:
	input:
		fasta = lambda wildcards: get_input_fasta(wildcards.method),
		ref = lambda wildcards: REF[wildcards.build]
	output:
		bam = output_file
	params:
		minimap2 = tools['minimap2'],
		samtools = tools['samtools']
	shell:
		"""
		{params.minimap2} -a {input.ref} {input.fasta} | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
		"""

rule All:
	input:
		[ output_file.format(tech=tech, build=build, acronym=acronym, method=method)
			for tech in TECHNOLOGIES
			for build in ["GRCh38"]
			for acronym in acronyms
			for method in ["flye-racon-medaka", "canu-racon-medaka", "flye", "canu"]
			]
