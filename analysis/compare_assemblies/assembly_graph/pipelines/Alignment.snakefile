import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

tools = {
	"minimap2": "/well/ont/apps/minimap2-v2.14/minimap2/minimap2",
	"samtools": "/apps/well/samtools/1.4.1/bin/samtools"
}

REF_37 = ""
REF_38 = "/well/longread/projects/reference/GRCh38/ont/genome.mmi"


def get_output_file(build, acronym, method):
	return ROOT + "/analysis/compare_assemblies/data/{build}/{acronym}/{method}/alignment.bam".format(build=build, acronym=acronym, method=method)

rule Align_contigs_to_ref:
	input:
		fasta = ROOT + "/analysis/de_novo_assembly/data/{build}/{acronym}/{method}/Output.symlink.fasta",
		ref = lambda wildcards: REF_38 if wildcards.build == "GRCh38" else REF_37
	output:
		bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{build}/{acronym}/{method}/{method}.bam",
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
		output_file = [ ROOT + "/analysis/compare_assemblies/data/{build}/{acronym}/{method}/{method}.bam".format(build=build, acronym=acronym, method=method) 
			for build in ["GRCh38"] 
			for acronym in ['HLA', 'TRA']
			for method in ['flye-racon-medaka', 'canu-racon-medaka', 'wtdbg2-racon-medaka']]


