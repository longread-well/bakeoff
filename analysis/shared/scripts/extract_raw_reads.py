from Bio import SeqIO
import argparse
import os

# Parse command line arguments
parser = argparse.ArgumentParser(description='Extract raw read sequences of mapped reads from FASTQ file')
parser.add_argument('-mapped', required = True, help = "FASTQ file of mapped reads")
parser.add_argument('-all', required = True, help = "FASTQ file of all reads")
parser.add_argument('-out', required = True, help = "Output FASTQ file of raw read sequences of mapped reads")
args = parser.parse_args()

mapped_read_file = args.mapped
all_read_file = args.all
output_read_file = args.out

print("Mapped reads: {mapped}\nAll reads: {all}\nOutput: {out}".format(
    mapped = mapped_read_file, all = all_read_file, out = output_read_file
))

# Extract raw read sequences
print("Indexing raw reads ...")
all_reads_dict = SeqIO.index(all_read_file, "fastq")
print("Extracting mapped reads ...")
mapped_reads = [all_reads_dict[read.id] for read in SeqIO.parse(mapped_read_file, 'fastq')]
print("Writing output file ...")
SeqIO.write(mapped_reads, output_read_file, "fastq")
print("Check MD5 ...")
hash1 = os.system("md5sum %s" % mapped_read_file)
hash2 = os.system("md5sum %s" % output_read_file)
