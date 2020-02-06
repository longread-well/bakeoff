import pandas as pd
import matplotlib.pyplot as plt
import os
import tempfile
import argparse

JELLYFISH = "/well/longread/users/akl399/bin/Jellyfish/jellyfish-linux"
TEMP_DIR = "/well/longread/users/akl399/temp/"

def jellyfish(fasta, k, canonical = False, size = "10M", thread = 10):
    print("Counting kmers")
    count_file = tempfile.NamedTemporaryFile(dir = TEMP_DIR).name
    count_command = "{jellyfish} count -m {k} -s {size} -t {thread} {c} {fasta} -o {output_file}".format(
        jellyfish = JELLYFISH,
        k = k,
        size = size,
        thread = thread,
        c = "-C" if canonical else "",
        fasta = fasta,
        output_file = count_file
    )
    os.system(count_command)
    histo_file = tempfile.NamedTemporaryFile(dir = TEMP_DIR).name
    histo_command = "{jellyfish} histo {count_file} > {histo_file}".format(
        jellyfish = JELLYFISH,
        count_file = count_file,
        histo_file = histo_file
    )
    os.system(histo_command)
    if os.path.exists(count_file):
        os.remove(count_file)
    return histo_file

def histogram(histo_file, output_file):
    print("Drawing histogram")
    data = pd.read_csv(histo_file, sep = ' ', names = ['occurance', 'kmer_number'])
    plt.plot(data.occurance, data.kmer_number)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('k-mer occurance')
    plt.ylabel('Number of k-mers')
    plt.savefig(output_file, dpi=300)
    if os.path.exists(histo_file):
        os.remove(histo_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot histogram of k-mer occurance. A wrapper of Jellyfish')
    parser.add_argument('-i', required = True, help = "Input FASTA sequence")
    parser.add_argument('-o', required = True, help = "Output PNG histogram")
    parser.add_argument('-k', default = 16, type = int, help = "Length of unique k-mers to be counted")
    parser.add_argument('--canonical', default = False, type = bool, help = "Treat reverse complement of k-mers as identical. Use --canonical for sequencing reads but not for reference genome")
    parser.add_argument('--size', default = "10M", help = "Estimated total number of unique k-mers. Roughly equals to genome size for reference genome, or genome size + genome size * coverage * error rate * k for sequencing reads")
    args = parser.parse_args()
    print(args)


    histo_file = jellyfish(args.i, args.k, args.canonical)
    histogram(histo_file, args.o)
