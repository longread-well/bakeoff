#!/usr/bin/env python
# coding: utf-8

import os

root = os.getcwd()
kmer_plot = "/well/longread/users/akl399/bakeoff/analysis/shared/scripts/kmer_plot.py"

cases = dict(deletion = ['XABCY', 'XACY'],
             insertion = ['XACY', 'XABCY'],
             inversion = ["XAY", "XMY"],
             one_more_copy = ['XAY', 'XAAY'],
             one_less_copy = ['XAAY', 'XAY'],
             two_copies = ['XAAY', 'XAAY'],
             three_copies = ['XAAAY', 'XAAAY'],
             unique = ['XAY', 'XAY'],
             repeat = ['XAYAZ', 'XAYAZ'])

for name, keys in cases.items():
    command = """
    set +u; source activate /well/longread/users/akl399/env/bakeoff; set -u
    python {kmer_plot} -x "{ref}" -y "{contig}" -k 50 -o "{png}" --figsize 4,4
    """.format(ref = os.path.join(root, 'sample_sequences', keys[0] + '.fasta'),
               contig = os.path.join(root, 'sample_sequences', keys[1] + '.fasta'),
               png = os.path.join(root, 'plots', name + '.png'),
               kmer_plot = kmer_plot)
    os.system(command)
