import os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

BUILD = dict(
    GRCh37 = ROOT + "/well/longread/projects/reference/GRCh37/10X/refdata-b37-2.1.0/fasta/genome.fa",
    GRCh38 = ROOT + "/well/longread/projects/reference/GRCh38/10X/refdata-GRCh38-2.1.0/fasta/genome.fa"
)

output_path = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta"

rule Extract_region_fasta:
    input:
        build = lambda wildcards: BUILD[wildcards.build]
    output:
        fasta = ROOT + "/analysis/
