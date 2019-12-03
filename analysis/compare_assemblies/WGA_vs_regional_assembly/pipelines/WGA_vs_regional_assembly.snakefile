import os, sys

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
include: ROOT + "/analysis/shared/scripts/initialize.py"

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

rule Symlink_regional_assemblies:
    input:
        bam = ROOT + "/analysis/compare_assemblies/align_contigs_to_ref/data/{tech}/{build}/{acronym}/{method}.bam"
    output:
        symlink = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}_regional.symlink.bam"
    shell:
        """
        ln -s {input.bam} {output.symlink}
        """

rule Visualize_assemblies:
    input:
        WGA_bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/Wtdbg2_WGA.bam",
        regional_bam = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}/regional_assembly.symlink.bam"
    output:
        png = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/summary/assemblies.png"
    params:
        visualize = tools['visualize_bam_file'],
        input_dir = os.path.join(ROOT, "analysis/compare_assemblies/WGA_vs_regional_assembly/data/WGA_vs_regional_assembly/{tech}/{build}/{acronym}/"),
        output_dir = os.path.join(ROOT, "analysis/compare_assemblies/WGA_vs_regional_assembly/data/WGA_vs_regional_assembly/{tech}/{build}/{acronym}/summary/"),
        chromosome = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["chromosome"],
        start = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["start"],
        end = lambda wildcards: get_region_definition(wildcards.acronym, wildcards.build, regions)["end"],
    shell:
        """
        set +u; source activate /users/todd/akl399/bin/miniconda/envs/Longread; set -u
        python {params.visualize} -i {params.input_dir} -o {params.output_dir} -r {params.chromosome}:{params.start}-{params.end}
        """

output_file = ROOT + "/analysis/compare_assemblies/WGA_vs_regional_assembly/data/{tech}/{build}/{acronym}/{method}_WGA.bam"
rule All:
    input:
        [output_file.format(tech=tech, build=build, acronym=acronym, method=method) for tech in TECHNOLOGIES for build in ['GRCh38'] for acronym in acronyms] for method in ['Wtdbg2']
