import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_file = ROOT + "/analysis/compare_assemblies/visualize_bam_file/data/{tech}/{build}/{acronym}/assemblies.png"

def get_bam_input(tech, build, acronym, method):
    if tech in ["ONT", "PB-CCS", "PB-CLR"]:
        bam_file = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/regional_assembly/{tech}/{build}/{acronym}/{method}.bam".format(
        tech = tech, build = build, acronym = acronym, method = method)
    if tech == "10X":
        bam_file = ROOT + "/analysis/compare_assemblies/extract_10X_contigs/data/{build}/10X_{acronym}.bam".format(
        build = build, acronym = acronym)

    return bam_file

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	return result[0]

rule Group_bam_files:
    input:
        bam = lambda wildcards: get_bam_input(wildcards.tech, wildcards.build, wildcards.acronym, wildcards.method),
        bai = lambda wildcards: get_bam_input(wildcards.tech, wildcards.build, wildcards.acronym, wildcards.method) + '.bai',
    output:
        bam = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/{tech} + {method}.bam",
        bai = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/{tech} + {method}.bam.bai",
    shell:
        """
        ln -s "{input.bam}" "{output.bam}"
        ln -s "{input.bai}" "{output.bai}"
        """

# Modify as needed
def get_summary_input():
    input_files = []
    template = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/{tech} + {method}.bam"
    for tech in ['ONT', 'PB-CCS', 'PB-CLR', '10X']:
        for build in ['GRCh38']:
            for acronym in acronyms:
                if tech in ['ONT']:
                    methods = ['Flye_Medaka', 'Canu_Medaka']
                elif tech in ['PB-CCS', 'PB-CLR']:
                    methods = ['Flye', 'Canu']
                elif tech in ['10X']:
                    methods = ['Supernova']

                if tech == 'PB-CCS' and acronym in ['GYP', 'IL27']:
                    methods = ['Flye', 'Canu', 'Canu_Purge']
                for method in methods:
                    input_files.append(template.format(build = build, acronym = acronym, tech = tech, method = method))

    return input_files

rule Summarize:
    input:
        bam = lambda wildcards: get_summary_input(),
    output:
        fig = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/summary/assemblies.png",
    params:
        visualize = tools['visualize_bam_file'],
        input_dir = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/",
        output_dir = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/summary/",
        chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
        start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
        end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
    shell:
        """
        set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
        python {params.visualize} -i {params.input_dir} -o {params.output_dir} -r {params.chromosome}:{params.start}-{params.end}
        """

def get_all_input():
    input_files = []
    template = ROOT + "/analysis/compare_assemblies/regional_summary/data/{build}/{acronym}/summary/assemblies.png"
    for build in ['GRCh38']:
        for acronym in acronyms:
            input_files.append(template.format(build = build, acronym = acronym))
    return input_files

rule All:
    input: lambda wildcards: get_all_input()
