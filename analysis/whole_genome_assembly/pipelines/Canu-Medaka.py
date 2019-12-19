#!/usr/bin/env python
import os
import argparse

ROOT = os.environ['BAKEOFF_ROOT']
canu = "/well/longread/users/akl399/bin/canu-1.9/Linux-amd64/bin/canu"
medaka = "/well/ont/apps/medaka/venv/bin/activate"
cluster_template = ROOT + "/analysis/whole_genome_assembly/scripts/cluster_template.sge.sh"

def get_fastq_file(tech):
	# TODO
	if tech == "ONT":
		return ROOT + "/resources/data/ONT/JK_HV31.ont.fastq"
	elif tech == "PB-CCS":
		return ROOT + "/resources/data/PB-CCS/JK_HV31.ccs.fastq"
	elif tech == "PB-CLR":
		return ROOT + "/resources/data/PB-CLR/JK_HV31.clr.fastq"
	else:
		raise Exception("Parameter error")

def generate_cluster_script(template, output_file, job_name, project_name, target_queue, slots, shell_command):
	template_file = open(template, 'r')
	assert template_file.mode == 'r'
	script = template_file.read()
	script = script.format(job_name = job_name, project_name = project_name, target_queue = target_queue,
		slots = slots, shell_command = shell_command)
	with open(output_file, 'w+') as output_file:
		output_file.write(script)

def Canu_WGA(fastq, output_path, tech, genome_size = '3.2g', use_grid = True):
    assert os.path.isfile(fastq)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    tech_flag = {"ONT": "nanopore-raw", "PB-CLR": "pacbio-raw", "PB-CCS": "pacbio-hifi"}[tech]

    command = """
    module load gcc/5.4.0
    module load java/1.8.0_latest
    """

    if use_grid:
        command += """
        {canu} -p canu -d {output_path} genomeSize={genome_size} 'gridOptions=-q himem.qh -P todd.prjc -S /bin/bash -N {tech}_Canu -e {tech}_Canu_error.log' 'gridEngineResourceOption=-pe shmem THREADS' -{tech_flag} {fastq}
        """.format(canu = canu, tech = tech, output_path = output_path, genome_size = genome_size, tech_flag = tech_flag, fastq = fastq)
    else:
        command += """
        {canu} -p {tech} -d {output_path} genomeSize={genome_size} 'useGrid=false' -{tech_flag} {fastq}
        """.format(canu = canu, tech = tech, output_path = output_path, genome_size = genome_size, tech_flag = tech_flag, fastq = fastq)

    command += """
    ln -rs {contigs} {symlink}
    """.format(contigs = os.path.join(output_path, "canu.contigs.fasta"), symlink = os.path.join(output_path, "Output.symlink.fasta"))

    print(command)
    os.system(command)

def Medaka_WGA(fastq, input_contigs, output_path):
	assert os.path.isfile(fastq)
	assert os.path.isfile(input_contigs)
	if not os.path.exists(output_path):
		os.makedirs(output_path)

	command = """
	set +u; source {medaka}; set -u
	medaka_consensus -i {fastq} -d {input_contigs} -o {output_path} -t 4 -m r941_prom_high
	ln -rs {output_contigs} {symlink}
	""".format(medaka = medaka, fastq = fastq, input_contigs = input_contigs, output_path = output_path,
	output_contigs = os.path.join(output_path, "consensus.fasta"),
	symlink = os.path.join(output_path, "Output.symlink.fasta"))

	print(command)
	cluster_script = os.path.join(output_path, "Medaka.sge.sh")
	generate_cluster_script(cluster_template, cluster_script, "Medaka_WGA", "todd.prjc", "himem.qh", 6, command)
	os.system("qsub %s" % cluster_script)

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Whole genome de novo assembly using Canu or Canu_Medaka (only for ONT)')
	parser.add_argument('--tech', required = True, help = "ONT | PB-CCS | PB-CLR")
	args = parser.parse_args()
	tech = args.tech
	assert tech in ['ONT', 'PB-CCS', 'PB-CLR']
	fastq = get_fastq_file(tech)
	canu_output_path = ROOT + "/analysis/whole_genome_assembly/data/{tech}/Canu/".format(tech = tech)
	canu_contigs = os.path.join(canu_output_path, "Output.symlink.fasta")

	if not os.path.isfile(canu_contigs):
		Canu_WGA(fastq, canu_output_path, tech)
	elif tech == "ONT":
		medaka_output_path = ROOT + "/analysis/whole_genome_assembly/data/{tech}/Canu_Medaka/".format(tech = tech)
		medaka_contigs = os.path.join(medaka_output_path, "Output.symlink.fasta")
		if not os.path.isfile(medaka_contigs):
			Medaka_WGA(fastq, canu_contigs, medaka_output_path)
		else:
			print("Nothing to do.")
