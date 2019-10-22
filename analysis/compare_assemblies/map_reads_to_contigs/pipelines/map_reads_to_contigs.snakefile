import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_file = ROOT + "/analysis/compare_assemblies/map_reads_to_contigs/data/{build}/{acronym}/{acronym}.bam"

rule Map_reads:
    input:
        reads = ROOT + "/analysis/de_novo_assembly/data/{build}/{acronym}/Extract_target_region/reads.fastq",
        ref = lambda wildcards: REF[wildcards.build]
    output:
        bam = output_file
    params:
        minimap2 = tools['minimap2'],
        samtools = tools['samtools']
    shell:
        """
		{params.minimap2} -a {input.ref} {input.reads} | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
		"""

rule All:
    input:
        [output_file.format(build = build, acronym = acronym)
            for build in builds
            for acronym in acronyms]
