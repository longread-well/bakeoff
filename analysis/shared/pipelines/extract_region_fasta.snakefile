import os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

BUILD = dict(
    GRCh37 = "/well/longread/projects/reference/GRCh37/10X/refdata-b37-2.1.0/fasta/genome.fa",
    GRCh38 = "/well/longread/projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa"
)

output_file = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta"

get_region_defination = lambda acronym, build: [region for region in regions if region['acronym'] == acronym and region['build'] == build][0]

rule Extract_region_fasta:
    input:
        build = lambda wildcards: BUILD[wildcards.build]
    output:
        fasta = output_file
    params:
        chromosome = lambda wildcards: get_region_defination(wildcards.acronym, wildcards.build)['chromosome'],
        start = lambda wildcards: get_region_defination(wildcards.acronym, wildcards.build)['start'],
        end = lambda wildcards: get_region_defination(wildcards.acronym, wildcards.build)['end'],
        samtools = tools['samtools']
    shell:
        """
        {params.samtools} faidx {input.build} {params.chromosome}:{params.start}-{params.end} > {output.fasta}
        """

rule All:
    input:
        fasta = [output_file.format(build = build, acronym = acronym)
        for build in ['GRCh37', 'GRCh38']
        for acronym in acronyms
        ]
