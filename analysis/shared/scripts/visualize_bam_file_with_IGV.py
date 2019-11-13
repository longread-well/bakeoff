#!/usr/bin/env python

# This script generates a txt file containing IGV commands to visualize all bam files within the target dir.
# Absolute path will be used.

import os, csv, sys

if len(sys.argv) == 3:
    target_dir = sys.argv[1]
    regions_file = sys.argv[2]
else:
    print('Usage: python visualize_bam_file_with_IGV.py <target_dir> <regions_file>\nNote: Only absolute paths are recognized for target_dir.')
    exit()

def load_regions(regions_file):
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
	with open( regions_file, newline = '' ) as file:
		regionReader = csv.DictReader( file, delimiter = '\t' )
		for region in regionReader:
			regions.append( dict( region ))
	return regions

regions = load_regions(regions_file)
acronyms = [region['acronym'] for region in regions]

class BAM_file(object):
    def __init__(self, path, filename):
        self.path = path
        self.filename = filename

        global acronyms
        for acronym in acronyms:
            if "/" + acronym in self.path:
                self.acronym = acronym
        assert self.acronym

        if "/GRCh37" in self.path:
            self.build = "37"
        elif "/GRCh38" in self.path:
            self.build = "38"
        else:
            raise Exception("Invalid path")

        global regions
        for region in regions:
            if region['acronym'] == self.acronym and region['build'] == "GRCh" + self.build:
                self.chromosome, start, end = region['chromosome'], int(region['start']), int(region['end'])
                delta = end - start
                start = round(start - 0.1 * delta)
                end = round(end + 0.1 * delta)
                self.start = str(start)
                self.end = str(end)

                if "chr" not in self.chromosome:
                    self.chromosome = "chr" + self.chromosome
                self.position = self.chromosome + ":" + self.start + "-" + self.end

bam_files = []
for dirpath, dirnames, filenames in os.walk(target_dir):
    for filename in filenames:
        if filename.split('.')[-1] == 'bam':
            bam_files.append(BAM_file(dirpath, filename))


IGV_COMMAND = """
new
genome hg{build}
load {path}/{filename}
snapshotDirectory {path}/IGV_screenshots/
goto {position}
snapshot {filename}.png
"""

output = ""

for bam_file in bam_files:
    output += (IGV_COMMAND.format(build=bam_file.build, path=bam_file.path, filename=bam_file.filename,
                                    position=bam_file.position))

output_file = open("IGV_commands.txt", "w+")
output_file.write(output)
output_file.close()
