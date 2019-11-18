import csv, os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	chromosome = result[0]['chromosome']
	start = result[0]['start']
	end = result[0]['end']
	region = chromosome + ":" + start + "-" + end
	return region

rule Align_10X_contigs_to_ref:
	input:
		fasta = "/well/longread/projects/wgs_10x/asm/JK_HV31_3/pseudohap.fasta.gz",
		ref = lambda wildcards: REF[wildcards.build]
	output:
		bam = ROOT + "/analysis/compare_assemblies/extract_10X_contigs/data/align_10X_contigs_to_ref/{build}_10X_pseudohap.bam"
	params:
		minimap2 = tools['minimap2'],
		samtools = tools['samtools']
	shell:
		"""
		{params.minimap2} -a {input.ref} {input.fasta} | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
		"""

output_file = ROOT + "/analysis/compare_assemblies/extract_10X_contigs/data/{build}/10X_{acronym}.bam"

rule Subset_10x_contigs:
	input:
		bam = ROOT + "/analysis/compare_assemblies/extract_10X_contigs/data/align_10X_contigs_to_ref/{build}_10X_pseudohap.bam"
	output:
		bam = output_file
	params:
		region = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions),
		samtools = tools['samtools']
	shell:
		"""
		{params.samtools} view -b {input.bam} {params.region} -o {output.bam}
		{params.samtools} index {output.bam}
		"""

rule All:
	input:
		[output_file.format(build = build, acronym = acronym) for build in ['GRCh38'] for acronym in acronyms]
