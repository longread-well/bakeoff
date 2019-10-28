import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_input_file(tech, build, acronym, method):
	prefix = ROOT + "/analysis/de_novo_assembly/data/{tech}/{build}/{acronym}/{method}/".format(tech=tech, build=build, acronym=acronym, method=method)
	if method == "flye-racon-medaka":
		input_file = "Flye/assembly_graph.gfa"
	elif method == "canu-racon-medaka":
		input_file = "Canu/canu.contigs.gfa"
	else:
		raise Exception("Only methods flye-racon-medaka and canu-racon-medaka are supported.")
	return prefix + input_file


output_file = ROOT + "/analysis/compare_assemblies/assembly_graph/data/{tech}/{build}/{acronym}/{method}.svg"

rule Bandage:
	input:
		gfa = lambda wildcards: get_input_file(wildcards.tech, wildcards.build, wildcards.acronym, wildcards.method)
	output:
		svg = output_file
	params:
		bandage = tools['bandage']
	shell:
		"""
		{params.bandage} image {input.gfa} {output.svg} --scope entire
		"""


rule All:
	input:
		output_file = [ output_file.format(tech=tech, build=build, acronym=acronym, method=method)
			for tech in TECHNOLOGIES
			for build in ["GRCh38"]
			for acronym in acronyms
			for method in ['flye-racon-medaka', 'canu-racon-medaka']]
