#!/usr/bin/env python
# coding: utf-8
import matplotlib.pyplot as plt
import pandas as pd
from random import sample
import argparse
import tempfile
import os
import numpy as np
from Bio import SeqIO

SELFMAP = "/well/longread/users/akl399/bin/selfmap/selfmap_v2.1-dev"

def color_hash(string):
    '''
    Assign a unique, random and nice-looking color (rgb format) to the input string.
    '''
    import hashlib
    hash_value = int(hashlib.md5(string.encode('utf-8')).hexdigest()[:4], base = 36)
    h = hash_value % 1009 / 1009 # Hue range: 0-1
    s = hash_value % 1021 / 1021 * 0.2 + 0.65 # Saturation range: 0.65-0.85
    l = hash_value % 1031 / 1031 * 0.4 + 0.3 # Lightness range: 0.3-0.7

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
    return (r, g, b)

def coverage_plot(ax, coverage_file, colors = None, bin_size = 10000, transpose = False, position_offset = 0):
    data = pd.read_csv(coverage_file, sep = '\t', header = None, names = ['sequence_name', 'position', 'coverage'])
    data['sequence_name'] = data['sequence_name'].apply(lambda s: s.split(':')[0])
    sequence_names = set(data.sequence_name)

    max_coverage = 0
    for sequence_name in sequence_names:
        subset = data.loc[data["sequence_name"] == sequence_name, :]
        if len(subset) < bin_size:
            continue
        coverage = list(np.convolve(subset.coverage, np.ones((bin_size,))/bin_size, mode='same'))
        position = list(subset.position + position_offset - 1)
        #print(len(coverage)); print(len(position)); print(sequence_name) # DEBUG
        max_coverage = max(np.max(coverage), max_coverage)
        color = colors[sequence_name] if colors else 'gray'
        plt.sca(ax)
        if not transpose:
            xy = position, coverage
            plt.ylabel("Coverage")
            plt.fill_between(position, coverage, 0, color = color, alpha = 0.2)
        else:
            xy = coverage, position
            plt.xlabel("Coverage")
            plt.fill_betweenx(position, coverage, 0, color = color, alpha = 0.2)
        plt.plot(*xy, color = color, alpha = 0.8)

    if not transpose:
        plt.ylim(0, max_coverage * 1.1)
    else:
        plt.xlim(0, max_coverage * 1.1)


def read_fasta_metadata(fasta_file):
    metadata = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence_name = record.id.split(":")[0]
        metadata[sequence_name] = dict(length = len(record.seq), offset = 0)
        if (record.id[:3] == 'chr' or record.id[0].isnumeric()) and ':' in record.id and "-" in record.id: # Check if record is part of a reference chromosome
            metadata[sequence_name]['offset'] = int(float(record.id.split(':')[1].split("-")[0]))
    print(metadata) # DEBUG
    return metadata

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Plot shared k-mers between two sequences')
    parser.add_argument('-x', required = True, help = "FASTA sequence to be plotted on x axis. Only one sequence can be included in the FASTA file.")
    parser.add_argument('-y', required = True, help = "FASTA sequence to be plotted on y axis. Multiple sequences can be included in the FASTA file. Scatter points will be colored according to their sequence name on y axis. ")
    parser.add_argument('-k', type = int, default = 50, help = "Minimum length of shared sequences")
    parser.add_argument('-o', default = "kmer_plot.png", help = "Output PNG file")
    parser.add_argument('--x_coverage', default = '', help = "Coverage file for x axis")
    parser.add_argument('--y_coverage', default = '', help = "Coverage file for y axis")
    parser.add_argument('--downsample', type = float, default = 1, help = "Downsample factor of scatter points")
    parser.add_argument('--xlabel', default = '', help = "X-axis label of k-mer plot")
    parser.add_argument('--ylabel', default = '', help = "Y-axis label of k-mer plot")
    parser.add_argument('--title', default = '', help = "Title of k-mer plot")
    parser.add_argument('--figsize', default = '6,6', help = "Figure size (width, height) in inches")
    parser.add_argument('--zoom', default = '', help = "Format: center_x,center_y,zoom_factor. If provided, the resulting k-mer plot will zoom in to centre at the selected position.")
    parser.add_argument('--alpha', type = float, default = '0.05', help = "Transparency of scatter points.")

    args = parser.parse_args()
    print(args)

    # Read sequence metadata
    x_metadata = read_fasta_metadata(args.x)
    y_metadata = read_fasta_metadata(args.y)
    x_offset = min(record['offset'] for record in x_metadata.values()) # selfmap outputs chromosome positions while samtools depth outputs 1-based positions, so offset values are read from FASTA headers to make positions consistant.
    y_offset = min(record['offset'] for record in y_metadata.values())

    # Run selfmap
    selfmap_output = tempfile.NamedTemporaryFile().name
    selfmap_command = "{selfmap} -sequence x={x} y={y} -kmer-size {k} -o {output}".format(
        selfmap = SELFMAP, x = args.x, y = args.y, k = args.k, output = selfmap_output
    )
    print("Running selfmap ...\n%s" % selfmap_command)
    os.system(selfmap_command)

    # Plot selfmap output
    print("Reading selfmap output ...")
    data = pd.read_csv(selfmap_output, sep = "\t", header = 8)
    if os.path.exists(selfmap_output):
        os.remove(selfmap_output)
        print("Temp file removed")

    points = []
    colors = {contig_name: color_hash(contig_name) for contig_name in y_metadata.keys()}
    print(colors) ## DEBUG

    check_chromosomes = lambda row: str(row['chromosome']).split(':')[0] == 'x' and str(row['other_chromosome']).split(':')[0] == 'y'
    points = [dict(x = row['position'], y = row['other_position'], contig_name = row['other_chromosome'].split(' ')[0].split(':')[1])
              for index, row in data.iterrows() if check_chromosomes(row)]

    points_downsampled = sample(points, round(len(points) * args.downsample))
    for point in points:
        point['color'] = colors[point['contig_name']]

    X = [point['x'] for point in points_downsampled]
    Y = [point['y'] for point in points_downsampled]
    C = [point['color'] for point in points_downsampled]

    print("Making k-mer plot (%s scatter points) ..." % len(points_downsampled))
    figsize = [int(s) for s in args.figsize.split(",")]
    if args.x_coverage and args.y_coverage:
        fig, ax_grid = plt.subplots(2, 2, gridspec_kw={'width_ratios': [4, 1], "height_ratios": [1, 4]}, figsize = figsize)
        main_axes, top_axes, right_axes = ax_grid[1,0], ax_grid[0,0], ax_grid[1,1]
        fig.delaxes(ax_grid[0,1])
    elif args.x_coverage and not args.y_coverage:
        fig, ax_grid = plt.subplots(2, 1, gridspec_kw={"height_ratios": [1, 4]}, figsize = figsize)
        main_axes, top_axes = ax_grid[1, 0], ax_grid[0, 0]
    elif args.y_coverage and not args.x_coverage:
        fig, ax_grid = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4, 1]}, figsize = figsize)
        main_axes, right_axes = ax_grid[0, 0], ax_grid[0, 1]
    else:
        fig = plt.figure(figsize = figsize)
        main_axes = plt.gca()


    main_axes.set_axisbelow(True)
    plt.sca(main_axes)
    plt.grid(color='whitesmoke', linestyle='-', linewidth=1)
    plt.scatter(X, Y, s = 0.5, c = C, alpha = args.alpha)


    x_min = min([contig['offset'] for contig in x_metadata.values()])
    x_max = max([contig['offset'] + contig['length'] for contig in x_metadata.values()])
    y_min = min([contig['offset'] for contig in y_metadata.values()])
    y_max = max([contig['offset'] + contig['length'] for contig in y_metadata.values()])

    if args.zoom:
        center_x, center_y, zoom_factor = [float(x) for x in args.zoom.split(",")]
        delta_x, delta_y = (x_max - x_min) / zoom_factor, (y_max - y_min) / zoom_factor
        x_min, x_max = center_x - delta_x, center_x + delta_x
        y_min, y_max = center_y - delta_y, center_y + delta_y

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)

    xl = args.xlabel if args.xlabel else os.path.basename(args.x)
    yl = args.ylabel if args.ylabel else os.path.basename(args.y)
    plt.xlabel(xl)
    plt.ylabel(yl)
    if args.title:
        plt.title(args.title)

    # Plot read coverage
    if args.x_coverage:
        print("Plotting read coverage on x axis ...")
        coverage_plot(top_axes, args.x_coverage, position_offset = x_offset)
        plt.sca(top_axes)
        plt.xticks([])
        plt.xlim(x_min, x_max)
    if args.y_coverage:
        print("Plotting read coverage on y axis ...")
        coverage_plot(right_axes, args.y_coverage, colors, transpose=True, position_offset = y_offset)
        plt.sca(right_axes)
        plt.yticks([])
        plt.ylim(y_min, y_max)

    fig.tight_layout()
    plt.savefig(args.o, dpi = 300)
