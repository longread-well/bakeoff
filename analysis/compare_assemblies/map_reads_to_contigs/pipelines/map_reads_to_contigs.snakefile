import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_output_file(tech, build, acronym, method):
	return ROOT + "/analysis/compare_assemblies/map_reads_to_contigs/data/{tech}/{build}/{acronym}/{method}.bam".format(tech=tech, build=build, acronym=acronym, method=method)

rule Map_reads_to_contigs:
    input:
        reads = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq",
        contigs = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta"
    output:
        bam = ROOT + "/analysis/compare_assemblies/map_reads_to_contigs/data/{tech}/{build}/{acronym}/{method}.bam",
        coverage = ROOT + "/analysis/compare_assemblies/map_reads_to_contigs/data/{tech}/{build}/{acronym}/{method}.tsv"
    params:
        minimap2 = tools['minimap2'],
        samtools = tools['samtools'],
        tech_flag = lambda wildcards: {"PB-CLR": "map-pb", "PB-CCS": "asm20", "ONT": "map-ont"}[wildcards.tech]
    shell:
        """
		{params.minimap2} -t 4 -ax {params.tech_flag} {input.contigs} {input.reads} --secondary=no | {params.samtools} sort -m 1G -o {output.bam} -T tmp.ali
        {params.samtools} index {output.bam}
        {params.samtools} depth {output.bam} > {output.coverage}
		"""

output_files = []
for tech in ["ONT", "PB-CCS", "PB-CLR", 'CLR+ONT']:
	for build in ["GRCh38"]:
		for acronym in acronyms:
			if tech == "ONT":
				methods = ["Flye", "Flye_Medaka", "Canu", "Canu_Medaka", "Wtdbg2", "Wtdbg2_Medaka", "Canu_Purge"]
			elif tech == "PB-CCS":
				methods = ["Flye", "Canu", "Wtdbg2", "Canu_Purge"]
			elif tech == "PB-CLR":
				methods = ["Flye", "Canu", "Wtdbg2", "Canunc", "Canu_Purge"]
			elif tech == 'CLR+ONT':
				methods = ['Canu'] # TODO
			for method in methods:
				output_files.append(get_output_file(tech, build, acronym, method))

rule All:
    input:
        output_files
