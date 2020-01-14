import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_output_file(tech, build, acronym):
	return ROOT + "/analysis/compare_assemblies/map_reads_to_reference/data/{tech}/{build}/{acronym}.bam".format(tech = tech, build = build, acronym = acronym)

rule Map_reads_to_contigs:
    input:
        reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq",
        ref = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta"
    output:
        bam = ROOT + "/analysis/compare_assemblies/map_reads_to_reference/data/{tech}/{build}/{acronym}.bam",
        coverage = ROOT + "/analysis/compare_assemblies/map_reads_to_reference/data/{tech}/{build}/{acronym}.tsv"
    params:
        minimap2 = tools['minimap2'],
        samtools = tools['samtools'],
        tech_flag = lambda wildcards: {"PB-CLR": "map-pb", "PB-CCS": "asm20", "ONT": "map-ont"}[wildcards.tech]
    shell:
        """
		{params.minimap2} -t 4 -ax {params.tech_flag} {input.ref} {input.reads} --secondary=no | {params.samtools} sort -m 1G -o {output.bam} -T tmp.ali
        {params.samtools} index {output.bam}
        {params.samtools} depth {output.bam} > {output.coverage}
		"""

output_files = []
for tech in ["ONT", "PB-CCS", "PB-CLR"]:
	for build in ["GRCh38"]:
		for acronym in acronyms:
			output_files.append(get_output_file(tech, build, acronym))

rule All:
    input:
        output_files
