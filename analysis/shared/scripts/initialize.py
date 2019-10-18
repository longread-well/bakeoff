import csv

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
			('end_sequence', 'CCCTGGAAACCATTATTTTATTCTCTGCTTCTATGAATTTGACTATTTTAGATACCTCAAATAAGTGTAATCATGCATTATTTG')])
	"""
	regions = []
	with open( filename, newline = '' ) as file:
		regionReader = csv.DictReader( file, delimiter = '\t' )
		for region in regionReader:
			regions.append( dict( region ))
	return regions

regions = load_regions()
acronyms = [region['acronym'] for region in regions]
methods = ['flye-racon-medaka', 'canu-racon-medaka', 'wtdbg2-racon-medaka']

tools = {
	"minimap2": "/well/ont/apps/minimap2-v2.14/minimap2/minimap2",
	"samtools": "/apps/well/samtools/1.4.1/bin/samtools",
	
}

