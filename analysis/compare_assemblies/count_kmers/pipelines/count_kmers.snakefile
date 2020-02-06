import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

K = 16

def get_fastq_file(tech):
	if tech == "ONT":
		return ROOT + "/resources/data/ONT/JK_HV31.ont.fastq"
	elif tech == "PB-CCS":
		return ROOT + "/resources/data/PB-CCS/JK_HV31.ccs.fastq"
	elif tech == "PB-CLR":
		return ROOT + "/resources/data/PB-CLR/JK_HV31.clr.fastq"
	elif tech == "10X":
		return "/well/longread/projects/wgs_10x/data/raw/11589MCpool01-N__JK_HV31_3/*.fastq.gz"
	else:
		raise Exception("Parameter error")

rule count_kmers_in_reads:
    input:
        fasta = lambda wildcards: get_fastq_file(wildcards.tech),
    output:
        png = ROOT + "/analysis/compare_assemblies/count_kmers/data/sequencing_reads/{tech}.png"
    params:
        count_kmers = tools['count_kmers'],
        k = K,
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.count_kmers} -i {input.fasta} -o {output.png} -k {params.k} --canonical True --size 150G
        """

rule count_kmers_in_reference:
    input:
        fasta = ROOT + "/analysis/shared/data/reference/{build}.fasta",
    output:
        png = ROOT + "/analysis/compare_assemblies/count_kmers/data/reference_genome/{build}.png"
    params:
        count_kmers = tools['count_kmers'],
        k = K,
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.count_kmers} -i {input.fasta} -o {output.png} -k {params.k} --canonical True --size 3G
        """


output_files = []
for build in ['GRCh38', 'GRCh37']:
    output_files.append(ROOT + "/analysis/compare_assemblies/count_kmers/data/reference_genome/{build}.png".format(
        build = build
    ))
for tech in ['PB-CCS', 'PB-CLR', 'ONT']:
    output_files.append(ROOT + "/analysis/compare_assemblies/count_kmers/data/sequencing_reads/{tech}.png".format(
        tech = tech
    ))

rule All:
    input:
        output_files
