import csv
import os

ROOT = os.environ[ 'BAKEOFF_ROOT' ]

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
builds = ["GRCh37", "GRCh38"]

tools = dict(
	bandage = "/well/longread/users/akl399/bin/Bandage_CentOS_static_v0_8_1/Bandage",
	flye = "/well/ont/apps/Flye-2.5/bin/flye",
	racon = "/well/ont/apps/racon/build/bin/racon",
	medaka = "/well/ont/apps/medaka/venv/bin/activate",
	canu = "/well/longread/users/akl399/bin/canu-1.9/Linux-amd64/bin/canu",
	wtdbg2 = "/well/ont/apps/wtdbg2",
	minimap2 = "/well/ont/apps/minimap2-v2.14/minimap2/minimap2",
	samtools = "/apps/well/samtools/1.4.1/bin/samtools",
	selfmap = "/well/longread/users/akl399/bin/selfmap_v2.1-dev",
	selfmap_plot = ROOT + "/analysis/shared/scripts/selfmap_plot.R",
	Rscript = '/apps/well/R/3.4.3-openblas-0.2.18-omp-gcc5.4.0/bin/Rscript',
	pdftoppm = "/users/todd/akl399/bin/miniconda/bin/pdftoppm",
	visualize_bam_file = ROOT + "/analysis/shared/scripts/visualize_bam_file_v2.py",
	extract_raw_reads = ROOT + "/analysis/shared/scripts/extract_raw_reads.py",
	bakeoff_env = "/well/longread/users/akl399/env/bakeoff",
	sda = "/well/longread/users/akl399/bin/SDA/SDA",
	kmer_plot = ROOT + "/analysis/shared/scripts/kmer_plot.py",
	purge_haplotigs = "/well/longread/users/akl399/env/purge_haplotigs_env",
	count_kmers = ROOT + "/analysis/shared/scripts/count_kmers.py",
)

REF = dict(
	# Generated with: /well/ont/apps/minimap2-v2.14/minimap2/minimap2 -d GRCh37_REF.mmi /well/longread/projects/reference/GRCh37/10X/refdata-b37-2.1.0/fasta/genome.fa
	GRCh37 = ROOT + "/analysis/shared/data/reference/GRCh37_REF.mmi",
	# Symlink to mmi file enerated by Hannah Roberts
	GRCh38 = ROOT + "/analysis/shared/data/reference/GRCh38_REF.mmi",
	)

TECHNOLOGIES = ['ONT', 'PB-CCS', 'PB-CLR']
