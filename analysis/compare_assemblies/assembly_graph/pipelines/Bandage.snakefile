import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

tools = {
	"bandage": "/well/longread/users/akl399/bin/Bandage_CentOS_static_v0_8_1/Bandage"
}

def get_input_file(build, acronym, method):
	prefix = ROOT + "/analysis/de_novo_assembly/data/{build}/{acronym}/{method}/".format(build=build, acronym=acronym, method=method)
	if method == "flye-racon-medaka":
		input_file = "Flye/assembly_graph.gfa"
	elif method == "canu-racon-medaka":
		input_file = "Canu/canu.contigs.gfa"
	else:
		raise Exception("Only methods flye-racon-medaka and canu-racon-medaka are supported.")
	return prefix + input_file

def get_output_file(build, acronym, method):
	return ROOT + "/analysis/compare_assemblies/data/{build}/{acronym}/{method}/assembly_graph.svg".format(build=build, acronym=acronym, method=method)

rule Bandage:
	input:
		gfa = lambda wildcards: get_input_file(wildcards.build, wildcards.acronym, wildcards.method)
	output:
		svg = ROOT + "/analysis/compare_assemblies/data/{build}/{acronym}/{method}/assembly_graph.svg"
	params:
		bandage = tools['bandage']
	shell:
		"""
		{params.bandage} image {input.gfa} {output.svg} --scope entire
		"""


rule All:
	input:
		output_file = [ get_output_file(build, acronym, method) 
			for build in ["GRCh37", "GRCh38"] 
			for acronym in acronyms
			for method in ['flye-racon-medaka', 'canu-racon-medaka']]


