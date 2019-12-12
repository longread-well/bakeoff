import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

def get_region_definition(acronym, build, regions):
	result = [ region for region in regions if region['build'] == build and region['acronym'] == acronym ]
	if len( result ) != 1:
		raise Exception( "No region with  acronym %s and build %s can be found." % ( acronym, build ) )
	chromosome = result[0]['chromosome']
	start = result[0]['start']
	end = result[0]['end']
	region = chromosome + ":" + start + "-" + end
	return region

def get_fasta_input(tech, build, method):
    assert build == "GRCh38" and method == "Wtdbg2" and tech in ['ONT', 'PB-CCS', 'PB-CLR']
    prefix = "/well/longread/shared/analysis/assembly/whole_genome/wtdbg2/"
    if tech == "ONT":
        tech_str = "nanopore"
    elif tech == "PB-CCS":
        tech_str = "pacbio-ccs"
    elif tech == "PB-CLR":
        tech_str = "pacbio-clr"
    return prefix + tech_str + "/JK_HV31." + tech_str + ".ctg.fa"

rule Map_WGA_contigs_to_reference:
    input:
        fasta = lambda wildcards: get_fasta_input(wildcards.tech, wildcards.build, wildcards.method),
        ref = lambda wildcards: REF[wildcards.build]
    output:
        bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/map_WGA_contigs_to_reference/{tech}/{build}/{method}.bam"
    params:
        minimap2 = tools['minimap2'],
        samtools = tools['samtools']
    shell:
        """
        {params.minimap2} -a {input.ref} {input.fasta} | {params.samtools} view -bS - | {params.samtools} sort -o {output.bam} -
		{params.samtools} index {output.bam}
        """

rule Subset_WGA_contigs:
    input:
        bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/map_WGA_contigs_to_reference/{tech}/{build}/{method}.bam"
    output:
        bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}_WGA.bam"
    params:
        region = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions),
        samtools = tools['samtools'],
    shell:
        """
		{params.samtools} view -b {input.bam} {params.region} -o {output.bam}
		{params.samtools} index {output.bam}
		"""

rule Fetch_regional_assemblies:
    input:
        bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/{method}.bam",
        bai = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/{method}.bam.bai"
    output:
        bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}_regional.bam",
        bai = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}_regional.bam.bai"
    shell:
        """
        cp {input.bam} {output.bam}
        cp {input.bai} {output.bai}
        """

rule Visualize_assemblies:
    input:
        WGA = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2_WGA.bam",
        regional = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2_regional.bam"
    output:
        png = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/summary/assemblies.png"
    params:
        visualize = tools['visualize_bam_file'],
        input_dir = lambda wildcards: ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/".format(tech=wildcards.tech, build = wildcards.build, acronym = wildcards.acronym),
        output_dir = lambda wildcards: ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/summary/".format(tech=wildcards.tech, build = wildcards.build, acronym = wildcards.acronym),
        region = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions),
    shell:
        """
        set +u; source activate /users/todd/akl399/bin/miniconda/envs/Longread; set -u
        python {params.visualize} -i {params.input_dir} -o {params.output_dir} -r {params.region}
        """

output_file = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/summary/assemblies.png"
rule All:
    input:
        [output_file.format(tech=tech, build=build, acronym=acronym) for tech in ['PB-CCS', 'PB-CLR'] for build in ['GRCh38'] for acronym in acronyms]
