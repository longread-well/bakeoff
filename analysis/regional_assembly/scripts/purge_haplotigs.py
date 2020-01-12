import os
import csv

ROOT = os.environ[ 'BAKEOFF_ROOT' ]
ph = "/well/longread/users/akl399/env/purge_haplotigs_env"
os.system("source activate {ph}".format(ph = ph))
minimap2 = "/well/ont/apps/minimap2-v2.14/minimap2/minimap2"
samtools = "/apps/well/samtools/1.4.1/bin/samtools"

def get_contigs(tech, build, acronym, method):
    return ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}/Output.symlink.fasta".format(
        tech = tech, build = build, acronym = acronym, method = method
    )

def get_reads(tech, build, acronym):
    return ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/Raw_reads/raw_reads.fastq".format(
        tech = tech, build = build, acronym = acronym
    )

def get_output_path(tech, build, acronym, method):
    return ROOT + "/analysis/regional_assembly/data/{tech}/{build}/{acronym}/{method}_Purge/".format(
        tech = tech, build = build, acronym = acronym, method = method
    )

def map_reads_to_contigs(reads, contigs, tech, output_path):
    tech_flag = {"PB-CCS": "asm20", "PB-CLR": "map-pb", "ONT": "map-ont"}[tech]
    os.system("cd {output_path} && {minimap2} -t 4 -ax {tech_flag} {contigs} {reads} --secondary=no | {samtools} sort -m 1G -o aligned.bam -T tmp.ali && touch map_reads_to_contigs.done".format(
        output_path = output_path, minimap2 = minimap2, tech_flag = tech_flag, contigs = contigs, reads = reads, samtools = samtools
    ))

def ph_hist(contigs, mapped_reads, output_path):
    os.system("source activate {ph} && cd {output_path} && purge_haplotigs  hist  -b {mapped_reads}  -g {contigs} && touch ph_hist.done".format(
        ph = ph, output_path = output_path, contigs = contigs, mapped_reads = mapped_reads
    ))

def ph_cov(l, m, h, output_path):
    os.system("source activate {ph} && cd {output_path} && purge_haplotigs  cov  -i aligned.bam.gencov -l {l} -m {m} -h {h} && touch ph_cov.done".format(
        ph = ph, output_path = output_path, l = l, m = m, h = h
    ))

def ph_purge(contigs, mapped_reads, output_path):
    if os.system("source activate {ph} && cd {output_path} && purge_haplotigs  purge  -g {contigs} -c coverage_stats.csv -dotplot {mapped_reads}".format(ph = ph, output_path = output_path, contigs = contigs, mapped_reads = mapped_reads)) == 0:
        os.system("cd {output_path} && touch ph_purge.done && ln -s curated.fasta Output.symlink.fasta".format(output_path = output_path))
    else:
        # If no contigs flagged as either suspects or artefacts, purge_haplotigs will through an error. Here we mask that error and use input contigs as the output
        os.system("cd {output_path} && touch ph_purge.done && ln -rs {contigs} Output.symlink.fasta".format(output_path = output_path, contigs = contigs))

def load_regions( filename = ROOT + '/resources/regions.tsv' ):
	"""
	Returns a list of regions, each as a OrderedDict.
	Example:
		OrderedDict([('name', 'T Cell Receptor Alpha'),
			('acronym', 'TRA'), ('build', 'GRCh38'),
			('chromosome', 'chr14'), ('start', '21170441'),
			('end', '22134544'),
			('start_sequence', 'CCTACCCCACACACTTATCACCCAGCAGGGAACCCTCAGGTTGGGCCCACAGCCCCCATT'),
			('end_sequence', 'CCCTGGAAACCATTATTTTATTCTCTGCTTCTATGAATTTGACTATTTTAGATACCTCAAATAAGTGTAATCATGCATTATTTG')，
			('length', '964103')
			])
	"""
	regions = []
	with open( filename, newline = '' ) as file:
		regionReader = csv.DictReader( file, delimiter = '\t' )
		for region in regionReader:
			regions.append( dict( region ))
	return regions

regions = load_regions()
acronyms = [region['acronym'] for region in regions]

for tech in ['ONT', 'PB-CCS', 'PB-CLR']:
    for build in ['GRCh38']:
        for acronym in acronyms:
            if tech == "ONT":
                methods = ['Canu']
            elif tech == 'PB-CCS':
                methods = ['Canu']
            elif tech == 'PB-CLR':
                methods = ['Canu', 'Canunc']
            for method in methods:
                print(tech + "/" + build + '/' + acronym + '/' + method)
                reads = get_reads(tech, build, acronym)
                contigs = get_contigs(tech, build, acronym, method)
                output_path = get_output_path(tech, build, acronym, method)
                if not os.path.exists(output_path):
                    os.makedirs(output_path)
                if not os.path.isfile(os.path.join(output_path, "map_reads_to_contigs.done")):
                    map_reads_to_contigs(reads, contigs, tech, output_path)
                    assert os.path.isfile(os.path.join(output_path, "map_reads_to_contigs.done"))
                mapped_reads = os.path.join(output_path, "aligned.bam")
                cutoffs = os.path.join(output_path, "ph_hist.done")
                if not os.path.isfile(cutoffs):
                    ph_hist(contigs, mapped_reads, output_path)
                    assert os.path.isfile(cutoffs)
                # User manually enter three cutoff values into ph_hist.done
                if os.stat(cutoffs).st_size > 1 and not os.path.isfile(os.path.join(output_path, "ph_cov.done")):
                    l, m, h = open(cutoffs, 'r').read().split(" ")
                    l, m, h = int(l), int(m), int(h)
                    assert h > m and m > l
                    ph_cov(l, m, h, output_path)
                    assert os.path.isfile(os.path.join(output_path, "ph_cov.done"))
                if os.path.isfile(os.path.join(output_path, "ph_cov.done")) and not os.path.isfile(os.path.join(output_path, "ph_purge.done")):
                    ph_purge(contigs, mapped_reads, output_path)
                    assert os.path.isfile(os.path.join(output_path, "ph_purge.done"))
