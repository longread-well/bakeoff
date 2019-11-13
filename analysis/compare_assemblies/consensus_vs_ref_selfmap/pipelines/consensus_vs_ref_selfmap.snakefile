import os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_path = ROOT + "/analysis/compare_assemblies/consensus_vs_ref_selfmap/data/{tech}/{build}/{acronym}/{method}/"

rule Selfmap:
    input:
        ref = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta",
        consensus = ROOT + "/analysis/de_novo_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta"
    output:
        k_100_pdf = output_path + "selfmap_k_100.pdf",
    params:
        selfmap_plot = tools['selfmap_plot'],
        Rscript = tools["Rscript"],
        pdftoppm = tools['pdftoppm'],
        chromosome = lambda wildcards: [region for region in regions if region['build'] == wildcards.build and region['acronym'] == wildcards.acronym][0]['chromosome']
    shell:
        """
        # Include selfmap path in $PATH
        export PATH=$PATH:/well/longread/users/akl399/bin/
		{params.Rscript} --vanilla {params.selfmap_plot} \
            --s1 '{wildcards.build}-{wildcards.acronym}={input.ref}' \
            --s2 '{wildcards.method}={input.consensus}' \
            --k 100 \
            --chromosome1 {params.chromosome} \
            --chromosome2 {params.chromosome} \
            --output {output.k_100_pdf}

        {params.pdftoppm} {output.k_100_pdf} {output.k_100_pdf} -png -rx 300 -ry 300
		"""

rule All:
    input:
        [output_path.format(tech = tech, build = build, acronym = acronym, method = method) + "selfmap_k_100.pdf"
            for tech in ["PB-CCS", "PB-CLR"]
            for build in ['GRCh38']
            for acronym in acronyms
            for method in ['flye-racon-medaka', 'canu-racon-medaka']
            ]
