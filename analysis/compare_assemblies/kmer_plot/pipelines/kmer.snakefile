import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

K = 50

output_file = ROOT + "/analysis/compare_assemblies/kmer_plot/data/{assembly_type}/{tech}/{build}/{acronym}/{method}.png"

rule Kmer_plot:
    input:
        ref = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta",
        contigs = ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta",
        ref_coverage = ROOT + "/analysis/compare_assemblies/map_reads_to_reference/data/{tech}/{build}/{acronym}.tsv",
        contigs_coverage = ROOT + "/analysis/compare_assemblies/map_reads_to_contigs/data/{tech}/{build}/{acronym}/{method}.tsv"
    output:
        png = output_file
    params:
        kmer_plot = tools['kmer_plot'],
        k = K,
        xlabel = lambda wildcards: '"Ref: %s"' % wildcards.acronym,
        ylabel = lambda wildcards: '"{tech} + {method}"'.format(tech = wildcards.tech, method = wildcards.method),
        title = '"K-mer sharing plot (k = %s)"' % K,
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.kmer_plot} -x {input.ref} -y {input.contigs} -k {params.k} -o {output.png} --x_coverage {input.ref_coverage} --y_coverage {input.contigs_coverage} --xlabel {params.xlabel} --ylabel {params.ylabel} --title {params.title}
        """

rule Ref_self_plot:
    input:
        ref = ROOT + "/analysis/shared/data/regions/sequence/{build}/{acronym}.fasta",
    output:
        png = ROOT + "/analysis/compare_assemblies/kmer_plot/data/Reference/{build}/{acronym}.png"
    params:
        kmer_plot = tools['kmer_plot'],
        k = K,
        xlabel = lambda wildcards: '"{build}: {acronym}"'.format(build = wildcards.build, acronym = wildcards.acronym),
        ylabel = lambda wildcards: '"{build}: {acronym}"'.format(build = wildcards.build, acronym = wildcards.acronym),
        title = '"K-mer sharing plot (k = %s)"' % K,
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.kmer_plot} -x {input.ref} -y {input.ref} -k {params.k} -o {output.png} --xlabel {params.xlabel} --ylabel {params.ylabel} --title {params.title}
        """

rule GRCh37_vs_GRCh38:
    input:
        ref37 = ROOT + "/analysis/shared/data/regions/sequence/GRCh37/{acronym}.fasta",
        ref38 = ROOT + "/analysis/shared/data/regions/sequence/GRCh38/{acronym}.fasta",
    output:
        png = ROOT + "/analysis/compare_assemblies/kmer_plot/data/Reference/GRCh37_vs_GRCh38/{acronym}.png"
    params:
        kmer_plot = tools['kmer_plot'],
        k = K,
        xlabel = lambda wildcards: '"GRCh38: %s"' % wildcards.acronym,
        ylabel = lambda wildcards: '"GRCh37: %s"' % wildcards.acronym,
        title = '"K-mer sharing plot (k = %s)"' % K,
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.kmer_plot} -x {input.ref38} -y {input.ref37} -k {params.k} -o {output.png} --xlabel {params.xlabel} --ylabel {params.ylabel} --title {params.title}
        """


output_files = []
for assembly_type in ['regional_assembly']:
    for tech in ['ONT', 'PB-CCS', 'PB-CLR', 'CLR+ONT']:
        if tech == "ONT":
            methods = ["Flye", "Flye_Medaka", "Canu", "Canu_Medaka", "Wtdbg2_Medaka", "Canu_Purge"]
        elif tech == "PB-CCS":
            methods = ["Flye", "Canu", "Wtdbg2", "Canu_Purge"]
        elif tech == "PB-CLR":
            methods = ["Flye", "Canu", "Wtdbg2", "Canu_Purge", "Canunc", "Canunc_Purge"]
        elif tech == 'CLR+ONT':
            methods = ['Canu']

        for build in ['GRCh38']:
            for acronym in acronyms:
                for method in methods:
                    output_files.append(output_file.format(assembly_type = assembly_type, tech = tech, build = build, acronym = acronym, method = method))

ref_output = ROOT + "/analysis/compare_assemblies/kmer_plot/data/Reference/{build}/{acronym}.png"
for build in ['GRCh37', 'GRCh38', 'GRCh37_vs_GRCh38']:
    for acronym in acronyms:
        output_files.append(ref_output.format(build = build, acronym = acronym))

rule All:
    input:
        output_files
