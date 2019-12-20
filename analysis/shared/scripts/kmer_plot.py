#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import pandas as pd
from random import sample
import argparse
import tempfile
import os

SELFMAP = "/well/longread/users/akl399/bin/selfmap/selfmap_v2.1-dev"

def hash_to_rgba(hash_value):
    h = hash_value % 1000 / 1000 # Hue range: 0-1
    s = hash_value % 1001 / 1001 * 0.2 + 0.65 # Saturation range: 0.65-0.85
    l = hash_value % 999 / 999 * 0.4 + 0.3 # Lightness range: 0.3-0.7
    a = 0.05

    def hsl_to_rgb(h, s, l):
        if s == 0:
            r = g = b = l
        else:
            def hue_to_rgb(p, q, t):
                if t < 0:
                    t += 1
                if t > 1:
                    t -= 1
                if t < 1/6:
                    return p + (q - p) * 6 * t
                if t < 1/2:
                    return q
                if t < 2/3:
                    return p + (q - p) * (2/3 - t) * 6
                return p

            q = l * (1 + s) if l < 0.5 else l + s - l * s
            p = 2 * l - q
            r = hue_to_rgb(p, q, h + 1/3)
            g = hue_to_rgb(p, q, h)
            b = hue_to_rgb(p, q, h - 1/3)

        return r, g, b

    r, g, b = hsl_to_rgb(h, s, l)
    return (r, g, b, a)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Plot shared k-mers between two sequences')
    parser.add_argument('-x', required = True, help = "Sequence to be plotted on x axis")
    parser.add_argument('-y', required = True, help = "Sequence to be plotted on x axis")
    parser.add_argument('-k', type = int, default = 50, help = "Minimum length of shared sequences")
    parser.add_argument('-o', default = "kmer_plot.png", help = "Output PNG file")
    parser.add_argument('--downsample', type = float, default = 1, help = "Downsample factor of scatter points")
    parser.add_argument('--xlabel', default = '', help = "X-axis label of k-mer plot")
    parser.add_argument('--ylabel', default = '', help = "Y-axis label of k-mer plot")
    parser.add_argument('--title', default = '', help = "Title of k-mer plot")

    args = parser.parse_args()
    selfmap_output = tempfile.NamedTemporaryFile().name

    # Run selfmap
    selfmap_command = "{selfmap} -sequence x={x} y={y} -kmer-size {k} -o {output}".format(
        selfmap = SELFMAP, x = args.x, y = args.y, k = args.k, output = selfmap_output
    )
    print("Running selfmap ...\n%s" % selfmap_command)
    os.system(selfmap_command)

    # Plot selfmap output
    print("Reading selfmap output ...")
    data = pd.read_csv(selfmap_output, sep = "\t", header = 8)
    points = []

    check_chromosomes = lambda row: str(row['chromosome']).split(':')[0] == 'x' and str(row['other_chromosome']).split(':')[0] == 'y'
    points = [dict(x = row['position'], y = row['other_position'], color = hash_to_rgba(hash(row['other_chromosome'])), contig_name = row['other_chromosome'])
              for index, row in data.iterrows() if check_chromosomes(row)]

    contig_names = set([point['contig_name'] for point in points]) # Select unique contig names
    y0 = {}
    for name in contig_names:
        y0[name] = min([point['y'] for point in points if point['contig_name'] == name]) # Find minimum y for each contig
    points = [dict(x = point['x'], y = point['y'] - y0[point['contig_name']], color = point['color']) for point in points] # Substract y_min from each contig to force y to start at 0

    points_downsampled = sample(points, round(len(points) * args.downsample))

    X = [point['x'] for point in points_downsampled]
    Y = [point['y'] for point in points_downsampled]
    C = [point['color'] for point in points_downsampled]

    print("Making k-mer plot (%s scatter points) ..." % len(points_downsampled))
    plt.figure(figsize = (8, 6))
    plt.gca().set_axisbelow(True)
    plt.grid(color='whitesmoke', linestyle='-', linewidth=1)
    plt.scatter(X, Y, s = 0.5, c = C)
    plt.xlim(min(X), max(X))
    plt.ylim(min(Y), max(Y))

    xl = args.xlabel if args.xlabel else os.path.basename(args.x)
    yl = args.ylabel if args.ylabel else os.path.basename(args.y)
    plt.xlabel(xl)
    plt.ylabel(yl)
    if args.title:
        plt.title(args.title)



    plt.savefig(args.o, dpi = 300)
    if os.path.exists(selfmap_output):
        os.remove(selfmap_output)
        print("Temp file removed")
