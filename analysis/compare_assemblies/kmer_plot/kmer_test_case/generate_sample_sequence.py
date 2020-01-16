#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np

bases = ('A', 'T', 'C', 'G')
np.random.seed(2020)
elements = {k: ''.join(list(np.random.choice(bases, size = 5000, replace = True)))
            for k in ('A', 'B', 'C', 'D', 'X', 'Y', 'Z')}

def reverse_complement(seq):
    r = reversed(list(seq))
    complement_nucleotide = dict(A = "T", T = "A", C = "G", G = "C")
    rc = ''.join([complement_nucleotide[base] for base in r])
    return rc

elements['M'] = reverse_complement(elements['A'])

def generate_chimera(key):
    fasta = ">Seq\n"
    for k in key:
        fasta += elements[k]
    return fasta

keys = ['XABCY', 'XACY', 'XAAAY', "XAAY", "XAY", 'XMY', 'XAYAZ']
sequences = {key: generate_chimera(key) for key in keys}

for key, fasta in sequences.items():
    output_file = "sample_sequences/%s.fasta" % key
    if not os.path.isfile(output_file):
        with open(output_file, 'w+') as f:
            f.write(fasta)
