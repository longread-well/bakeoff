import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

output_file = ROOT + "/analysis/compare_assemblies/visualize_bam_file/data/{tech}/{build}/{acronym}/assemblies.png"

def get_input_dir(tech, build, acronym):
    return ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/".format(tech = tech, build = build, acronym = acronym)

def get_output_dir(tech, build, acronym):
    return ROOT + "/analysis/compare_assemblies/visualize_bam_file/data/{tech}/{build}/{acronym}/".format(tech = tech, build = build, acronym = acronym)

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	return result[0]


rule visualize_bam_file:
    output:
        png = output_file
    params:
        input_dir = lambda wildcards: get_input_dir(wildcards.tech, wildcards.build, wildcards.acronym),
        output_dir = lambda wildcards: get_output_dir(wildcards.tech, wildcards.build, wildcards.acronym),
        chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
        start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
        end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
        visualize = tools['visualize_bam_file'],

    shell:
        """
        set +u; source activate /users/todd/akl399/bin/miniconda/envs/Longread; set -u
        python3 {params.visualize} -i {params.input_dir} -o {params.output_dir} -r {params.chromosome}:{params.start}-{params.end}
        """

rule All:
    input:
        [output_file.format(tech = tech, build = build, acronym = acronym)
        for tech in ["ONT", "PB-CCS", "PB-CLR"]
        for build in ["GRCh38"]
        for acronym in acronyms]
