import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

REF = dict(
	# Generated with: /well/ont/apps/minimap2-v2.14/minimap2/minimap2 -d GRCh37_REF.mmi /well/longread/projects/reference/GRCh37/10X/refdata-b37-2.1.0/fasta/genome.fa
	GRCh37 = ROOT + "/analysis/shared/data/reference/GRCh37_REF.mmi",
	# Symlink to mmi file enerated by Hannah Roberts
	GRCh38 = ROOT + "/analysis/shared/data/reference/GRCh38_REF.mmi",
	)

output_file = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{build}/{acronym}/{method}.bam"

rule Align_contigs_to_ref:
	input:
		fasta = ROOT + "/analysis/de_novo_assembly/data/{build}/{acronym}/{method}/Output.symlink.fasta",
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
		[ output_file.format(build=build, acronym=acronym, method=method)
			for build in ["GRCh37", "GRCh38"]
			for acronym in acronyms
			for method in methods
			]
