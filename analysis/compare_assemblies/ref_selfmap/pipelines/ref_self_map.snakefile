import os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_path = ROOT + "/analysis/compare_assemblies/ref_selfmap/data/{build}/{acronym}/"

rule Selfmap:
    input:
        fasta = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta"
    output:
        k_100_pdf = output_path + "selfmap_k_100.pdf"
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
            --s1 'ref1={input.fasta}' \
            --s2 'ref2={input.fasta}' \
            --k 100 \
            --chromosome1 {params.chromosome} \
            --chromosome2 {params.chromosome} \
            --output {output.k_100_pdf}

        {params.pdftoppm} {output.k_100_pdf} {output.k_100_pdf} -png -rx 300 -ry 300
		"""

rule All:
    input:
        [output_path.format(build = build, acronym = acronym) + "selfmap_k_100.pdf"
            for build in builds
            for acronym in acronyms
            ]
