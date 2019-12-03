#!/usr/bin/env python
# coding: utf-8
import pysam
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import pandas as pd
import numpy as np
import pandas as pd
import argparse
from matplotlib import cm


ROOT = os.environ[ 'BAKEOFF_ROOT' ]
GENE_ANNOTATION_FILE = ROOT + "/resources/GRCh38_Genes.tsv"

# Customizable parameters
NM_PERCENTAGE_THRESHOLD = 0.2 # Contigs will be hidden and excluded from statistics if NM_percentage >= NM_PERCENTAGE_THRESHOLD; parameter used to exclude mismapped contigs
INSERTION_THRESHOLD = 10 # Insertions will hidden if length < INSERTION_THRESHOLD
INSERTION_THRESHOLD_2 = 100 # Insertions will be highlighted if length >= INSERTION_THRESHOLD_2
DELETION_THRESHOLD = 1 # Delections will not be displayed if length < DELETION_THRESHOLD
ALIGN_CONTINGENCY = 10 # Blocks will not be displayed if length < ALIGN_CONTINGENCY

RECTANGLE_HEIGHT = 0.2 # Hight of rectangles representing contigs and genes

class Canvas:
    def __init__(self, target_region):
        self.tracks = []
        self.y = 0 # Use y for variable values and y0 for fixed values
        self.target_region = target_region

    def add_track(self, track):
        self.tracks.append(track)
        track.y0 = self.y
        if track.type == Assembly:
            pass
        elif track.type == Raw:
            pass
        elif track.type == Gene_Annotation:
            pass
        self.y += track.height

    def draw(self, x_limit=None):

        # Determine x limit
        if x_limit:
            x_min, x_max = x_limit
        else:
            x_min = None
            x_max = None
            for track in self.tracks:
                if track.type in ["Assembly", "Raw"]:
                    if x_min is None or x_min > track.left:
                        x_min = track.left
                    if x_max is None or x_max < track.right:
                        x_max = track.right

            delta_x = (x_max - x_min) * 0.03
            x_max += delta_x
            x_min -= delta_x

            # Prevent x limit from getting too wide
            x_max = min(x_max, self.target_region.right * 1.2 - self.target_region.left * 0.2)
            x_min = max(x_min, self.target_region.left * 1.2 - self.target_region.right * 0.2)

        # Determine y limit
        y_min = 0
        y_max = self.y


        # Create figure
        figure_height = max(3, 0.5 * self.y)
        plt.figure(figsize = (10, figure_height))

        # Plot tracks
        for track in self.tracks:
            track.draw(x_min, x_max)

        # Draw a pair of vertical lines marking our region of interest
        kwargs = dict(linewidth=1, linestyle='--', color='lightgray', zorder = 0.1)
        plt.gca().axvline(x = self.target_region.left, **kwargs)
        plt.gca().axvline(x = self.target_region.right, **kwargs)

        # Hide ticks on y axis
        plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        # Misc
        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        plt.xlabel("Reference postition")
        plt.tight_layout()

    def __repr__(self):
        output = "Target region: {region}\n\n".format(region = self.target_region)
        for track in self.tracks:
            output += track.__repr__()
        return output

class Track:
    def __init__(self, data):
        self.data = data
        self.type = data.type # Assembly, Raw, Gene_Annotation
        self.left, self.right = data.left, data.right
        self.height = self.data.get_height()
        self.label = self.data.name

    def draw(self, x_min, x_max):
        if self.type == "Gene_Annotation":
            self.data.draw(self.y0, x_min, x_max)
        else:
            self.data.draw(self.y0)

        # Draw a horizontal seperating line
        plt.gca().axhline(y = self.y0, linewidth=1, color='lightgray')

        # Track label
        label_x_position = 0.97 * x_max + 0.03 * x_min
        plt.text(label_x_position, self.y0 + 0.05, self.label, ha = 'right', va = 'bottom')

    def __repr__(self):
        return self.data.__repr__()


class Assembly:
    def __init__(self, name):
        self.name = name
        self.contigs = []
        self.type = "Assembly"
        self.left = None
        self.right = None

        self.bottom_padding = 0.7
        self.top_padding = 0.5
        self.delta_y = 0.5

    def add_contig(self, contig):
        self.contigs.append(contig)
        self.n = len(self.contigs)
        if self.left is None or self.left > contig.left:
            self.left = contig.left
        if self.right is None or self.right < contig.left:
            self.right = contig.right

    def get_height(self):
        assigned_contigs = assign_rows(self.contigs)
        n_row = max([assigned_contig[0] for assigned_contig in assigned_contigs]) + 1 # Row number is 0-indexed
        self.height = (n_row - 1)* self.delta_y + self.bottom_padding + self.top_padding
        return self.height

    def get_ref_coverage(self, target_region):
        ref_covered = [False for position in range(target_region.length)]
        coverage = pd.DataFrame(columns = ("i", "contig_length", "covered_percentage"))
        i = 0
        self.contigs.sort(key = lambda contig: get_overlaped_length(contig.region, target_region), reverse = True)
        for contig in self.contigs:
            i += 1
            contig_start = contig.left - target_region.left
            contig_end = contig.right - target_region.left
            if contig_start < 0:
                contig_start = 0
            if contig_end > target_region.length - 1:
                contig_end = target_region.length - 1
            for j in range(contig_start, contig_end + 1):
                ref_covered[j] = True
            covered_percentage = sum(ref_covered) / target_region.length
            coverage = coverage.append(
                pd.DataFrame({"i":[i], "contig_length":get_overlaped_length(contig.region, target_region), "covered_percentage":[covered_percentage]}))
        return coverage

    def get_NG50(self, target_region):
        coverage = self.get_ref_coverage(target_region)
        NG50 = None
        for index, row in coverage.iterrows():
            if row.covered_percentage >= 0.5:
                NG50 = row.contig_length
                break
        return NG50

    def get_NG90(self, target_region):
        coverage = self.get_ref_coverage(target_region)
        NG90 = None
        for index, row in coverage.iterrows():
            if row.covered_percentage >= 0.9:
                NG90 = row.contig_length
                break
        return NG90

    def get_NM(self, target_region):
        # TODO: only calculate bases with the target_region?
        pass

    def draw(self, y0):
        # Assign y0 for each contig
        for i, contig in assign_rows(self.contigs):
            contig.draw(y0 + self.bottom_padding + i * self.delta_y)

    def calculate_statistics(self):
        pass # TODO

    def __repr__(self):
        output = "Assembly {name} ({n} contigs)\n".format(name=self.name, n=len(self.contigs))
        for contig in self.contigs:
            output += contig.__repr__()
        return output

class Contig:
    def __init__(self):
        self.region = None
        self.blocks = []
        self.insertions = []
        self.deletions = []
        self.left = None
        self.right = None
        self.NM = None
        self.NM_percentage = None

    def add_block(self, block):
        self.blocks.append(block)
        if self.left is None or self.left > block.left:
            self.left = block.left
        if self.right is None or self.right < block.left:
            self.right = block.right
        self.length = self.right - self.left + 1

    def draw(self, y0):
        # Determine contig color
        color = cm.viridis(self.NM_percentage / 0.05) # Brighter colors denote higher NM%

        # Draw contig spine
        spine_height = 0.05
        plt.gca().add_patch(patches.Rectangle((self.left, y0 - 0.5 * spine_height),
                                       self.length, spine_height, facecolor = color, edgecolor=None, zorder=0.8))

        # Draw blocks
        for block in self.blocks:
            block.draw(y0, color)

        # Draw insertions and deletions
        if self.insertions:
            for insertion in self.insertions:
                insertion.draw(y0)
        if self.deletions:
            for deletion in self.deletions:
                deletion.draw(y0)

    def __repr__(self):
        output = "Contig ({region}, {m} blocks, {n} insertions, {p} deletions) \n".format(
            region = self.region.__repr__(), m = len(self.blocks), n = len(self.insertions), p = len(self.deletions))
        return output

class Block:
    def __init__(self, region):
        self.region = region
        self.left, self.right = region.left, region.right
        self.length = self.right - self.left + 1

    def draw(self, y0, color):
        if self.right - self.left >= ALIGN_CONTINGENCY:
            plt.gca().add_patch(patches.Rectangle((self.left, y0 - 0.5 * RECTANGLE_HEIGHT), self.length, RECTANGLE_HEIGHT,
             zorder=1, facecolor=color, edgecolor=None))

class Insertion:
    def __init__(self, position, length):
        self.position = position
        self.length = length
    def draw(self, y):
        if self.length < INSERTION_THRESHOLD:
            return
        if not self.position or not y:
            return
        if self.length >= INSERTION_THRESHOLD_2:
            marker_color = 'red'
        else:
            marker_color = 'purple'

        def align_marker(marker, halign='center', valign='middle',):
            """
            create markers with specified alignment.
            From: https://stackoverflow.com/questions/26686722/align-matplotlib-scatter-marker-left-and-or-right
            """
            from matplotlib import markers
            from matplotlib.path import Path
            if isinstance(halign, (bytes, str)):
                halign = {'right': -1., 'middle': 0., 'center': 0., 'left': 1.,}[halign]
            if isinstance(valign, (bytes, str)):
                valign = {'top': -1., 'middle': 0., 'center': 0., 'bottom': 1.,}[valign]
            bm = markers.MarkerStyle(marker)
            m_arr = bm.get_path().transformed(bm.get_transform()).vertices
            m_arr[:, 0] += halign / 2
            m_arr[:, 1] += valign / 2
            return Path(m_arr, bm.get_path().codes)

        plt.plot(self.position, y - 0.5 * RECTANGLE_HEIGHT, marker = align_marker("^", valign="top"), color=marker_color, markersize=8, markeredgewidth=0)

    def __repr__(self):
        return "Insertion(position = {position}, length = {length})".format(
            position = self.position, length = self.length)

class Deletion():
    def __init__(self, left, right):
        self.left = left
        self.right = right
        self.width = self.right - self.left + 1
    def __repr__(self):
        output = "Deletion(width = {width}) ".format(width = self.width)
        return output
    def draw(self, y):
        if self.width < DELETION_THRESHOLD:
            return
        ax = plt.gca()
        ax.add_patch(patches.Rectangle((self.left, y - 0.5 * RECTANGLE_HEIGHT), self.width, RECTANGLE_HEIGHT, zorder=2, facecolor='red', edgecolor=None))

class Raw:
    def __init__(self):
        self.type = "Raw"
        self.reads = []
        # TODO

    def add_read(self, read):
        self.reads.append(read)
        self.n = len(self.reads)

class Gene_Annotation:
    def __init__(self, target_region):
        self.type = "Gene_Annotation"
        self.genes = []
        self.name = "Genes"
        self.region = target_region
        self.chromosome, self.left, self.right = target_region.chromosome, target_region.left, target_region.right
        self.delta_y = 0.5
        self.bottom_padding = 0.7
        self.top_padding = 0.5


    def add_gene(self, gene):
        self.genes.append(gene)

    def get_height(self):
        assigned_genes = assign_rows(self.genes)
        n_row = max([assigned_gene[0] for assigned_gene in assigned_genes]) + 1 # Row number is 0-indexed
        self.height = (n_row - 1)* self.delta_y + self.bottom_padding + self.top_padding
        return self.height

    def draw(self, y0, x_min, x_max):
        for i, gene in assign_rows(self.genes):
            gene.draw(y0 + self.bottom_padding + i * self.delta_y, x_min, x_max)

    def __repr__(self):
        output = "Gene Annotation ({n} genes)".format(n = len(self.genes))
        return output

class Gene:
    def __init__(self, region, name):
        self.name = name
        self.region = region
        self.chromosome, self.left, self.right = region.chromosome, region.left, region.right
        self.length = self.right - self.left + 1
        self.exons = [] # TODO
        self.protein_coding = None

    def draw(self, y0, x_min, x_max):
        # Draw rectangles
        c = np.random.rand(3,) * 0.8
        plt.gca().add_patch(patches.Rectangle((self.left, y0 - 0.5 * RECTANGLE_HEIGHT), self.length, RECTANGLE_HEIGHT, zorder=1,
                                       facecolor=c, edgecolor=None, alpha=0.8))

        # Label genes
        label_top_padding = 0.05
        if self.length >= 0.005 * (x_max - x_min):
            center = (self.left + self.right)/2
            plt.text(center, y0 - 0.5 * RECTANGLE_HEIGHT - label_top_padding, self.name, ha='center', va='top', color=c, fontsize=3, rotation = 90, clip_on=True)

    def __repr__(self):
        output = "Gene {name} ({region})".format(name = self.name, region = self.region.__repr__())
        return output

class Region():
    def __init__(self, chromosome, left, right):
        self.chromosome, self.left, self.right = chromosome, left, right
        self.length = self.right - self.left + 1

    def __repr__(self):
        output = "Region {chromosome}:{left}-{right}".format(chromosome = self.chromosome, left = self.left, right = self.right)
        return output




def region_overlap(region1, region2):
    # Determine if two regions overlap with each other
    if region1.chromosome != region2.chromosome:
        return False
    if region1.right >= region2.left and region1.right <= region2.right:
        return True
    if region1.left >= region2.left and region1.left <= region2.right:
        return True
    if region1.left <= region2.left and region1.right >= region2.right:
        return True
    return False

def get_overlaped_length(region1, region2):
    # Return the number of overlapped bases between two regions
    if region_overlap(region1, region2):
        overlapped_left = max(region1.left, region2.left)
        overlapped_right = min(region1.right, region2.right)
        overlapped_length = overlapped_right - overlapped_left + 1
        return overlapped_length
    else:
        return 0

def plot_coverage(coverage, axis, label = ""):
    # (0, l1), (c1, l1), (c1, l2), (c2, l2), (c2, l3), ... , (cn, ln)
    x = [0]
    y = [coverage.iloc[0, 1]]
    for i in range(len(coverage) - 1):
        x.append(coverage.iloc[i, 2])
        y.append(coverage.iloc[i, 1])

        x.append(coverage.iloc[i, 2])
        y.append(coverage.iloc[i + 1, 1])

    x.append(coverage.iloc[-1, 2])
    y.append(coverage.iloc[-1, 1])

    plt.sca(axis)
    plt.plot(x, y, label = label)
    plt.xlabel("x")
    plt.ylabel("NGx")
    plt.legend(frameon=False, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    # Set x ticks to percentage format
    ax = plt.gca()
    ax.set_xticklabels(['{:,.0%}'.format(x) for x in ax.get_xticks()])

    plt.tight_layout()

def assign_rows(blocks):
    # Assign blocks (e.g. contigs/reads/genes) to appropriate rows to avoid overlapping
    # Returns a list of tuples: [(y1, triangle1), (y2, triangle2), ...]

    class Row:
        def __init__(self, i):
            self.i = i
            self.blocks = []

        def add_block(self, new_block):
            overlapped = False
            for assigned_block in self.blocks:
                if region_overlap(new_block.region, assigned_block.region):
                    overlapped = True
                    break
            if not overlapped:
                self.blocks.append(new_block)
                return True
            else:
                return False

    i = 0
    rows = [Row(i)]

    blocks.sort(key=lambda block: block.length, reverse=True)

    for new_block in blocks:
        assigned = False
        for row in rows:
            if row.add_block(new_block):
                assigned = True
                break
        if not assigned:
            i += 1
            new_row = Row(i)
            new_row.add_block(new_block)
            rows.append(new_row)

    assigned_blocks = []
    for row in rows:
        for block in row.blocks:
            assigned_blocks.append((row.i, block))

    return assigned_blocks

def parse_region(region):
    # region str format: chr14:105094527-106881340
    chromosome, interval = region.split(":")
    left_limit, right_limit = interval.split("-")
    left_limit = int(left_limit)
    right_limit = int(right_limit)
    return Region(chromosome, left_limit, right_limit)

def read_bam_file(bam_file, type, target_region, name = None):
    if not name:
        # Use bamfile name without extension as the name
        base_name = os.path.basename(bam_file)
        name = os.path.splitext(base_name)[0]

    bam = pysam.AlignmentFile(bam_file, 'rb')

    def read_contig_data(contig_data):
        contig = Contig()
        chromosome = contig_data.reference_name
        blocks_data = contig_data.get_blocks()
        for block_data in blocks_data:
            left, right = block_data[0], block_data[1]
            block = Block(Region(chromosome, left, right))
            contig.add_block(block)
        contig.NM = contig_data.get_tag('NM')
        contig.NM_percentage = contig.NM / contig.length
        contig.region = Region(chromosome, contig.left, contig.right)



        # Detect insertions
        aligned_pairs = contig_data.get_aligned_pairs()
        length = 0
        position = 0
        insertions = []
        for query, reference in aligned_pairs:
            if query and not reference:
                if length == 0: # Start of a new insertion
                    pass
                else: # Continue current insertion
                    pass
                length += 1

            else:
                if length == 0: # Matched pair; record position
                    position = reference
                else: # End current insertion
                    insertions.append(Insertion(position, length))
                    position = 0
                    length = 0
        contig.insertions = insertions

        # Detect deletions
        deletions = []
        left = None
        right = None
        for query, reference in aligned_pairs:
            if reference and not query:
                if not left: # Start of a new deletion
                    left = reference
                    right = left + 1
                else: # Continue current deletion
                    right = reference
            else:
                if not left: # Matched pair; record position
                    pass
                else: # End current insertion
                    deletions.append(Deletion(left, right))
                    left = None
                    right = None
        contig.deletions = deletions

        return contig

    if type == "Assembly":
        assembly = Assembly(name=name)
        for contig_data in bam.fetch():
            contig = read_contig_data(contig_data)
            if region_overlap(contig.region, target_region):
                assembly.add_contig(contig)
        track = Track(data=assembly)
    elif type == "Raw":
        raw = Raw(name=name)
        for read_data in bam.fetch():
            read = read_contig_data(read_data)
            if region_overlap(read.region, target_region):
                raw.add_read(read)
        track = Track(data=raw)

    return track

def load_genes(gene_file, target_region):
    gene_data = pd.read_table(gene_file)
    gene_data["tx_length"] = gene_data['txEnd'] - gene_data['txStart']
    gene_data.sort_values(by = ['tx_length'], ascending = False)
    gene_data = gene_data.drop_duplicates(subset = ['name2'])
    genes = []
    for index, row in gene_data.iterrows():
        gene = Gene(name = row["name2"], region =
                    Region(chromosome = row['chrom'], left = row['txStart'], right = row['txEnd']))
        gene.protein_coding = not row["cdsStart"] == row["cdsEnd"]

        if region_overlap(gene.region, target_region) and gene.protein_coding:
            genes.append(gene)

    gene_annotation = Gene_Annotation(target_region)
    for gene in genes:
        gene_annotation.add_gene(gene)

    return Track(gene_annotation)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Visualize regional assemblies and calculate assembly statistics')
    parser.add_argument('-i', required = True, help = "Path containing bam files")
    parser.add_argument('-o', required = True, help = "Path to store output image and statictics")
    parser.add_argument('-r', required = True, help = "Genome region of interest; format: chr14:21170441-22134544")
    args = parser.parse_args()
    input_dir = args.i
    input_dir = os.path.join(input_dir, '') # add tailing "/" if missing
    output_dir = args.o
    output_dir = os.path.join(output_dir, '') # add tailing "/" if missing
    target_region = parse_region(args.r)
    print("Input path: {i}\nOutput path: {o}\nTarget region: {r}".format(i = input_dir, o = output_dir, r = target_region))

    # Load assemblies
    bamfiles = [os.path.join(input_dir, f) for f in os.listdir(input_dir)
                if os.path.isfile(os.path.join(input_dir, f)) and f.split('.')[-1] == 'bam' and f[0] != '.']
    bamfiles.sort() # Display bamfiles with similar names next to each other
    print("Found {n} bam files".format(n = len(bamfiles)))
    for b in bamfiles:
        print(b)
    assembly_tracks = [read_bam_file(bamfile, "Assembly", target_region) for bamfile in bamfiles]

    canvas = Canvas(target_region)
    for track in assembly_tracks:
        canvas.add_track(track)

    # Load genes
    gene_annotation_track = load_genes(GENE_ANNOTATION_FILE, target_region)
    canvas.add_track(gene_annotation_track)

    # Visualize
    canvas.draw()
    print(canvas)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    plt.savefig(os.path.join(output_dir, "assemblies.svg"))
    plt.savefig(os.path.join(output_dir, "assemblies.png"), dpi = 300)

    plt.figure(figsize = (8, 4))

    # Coverage plot
    plt.figure(figsize=(7,4))
    for track in canvas.tracks:
        if track.type == 'Assembly':
            assembly = track.data
            coverage = assembly.get_ref_coverage(target_region)
            plot_coverage(coverage, plt.gca(), assembly.name)
    plt.savefig(os.path.join(output_dir, "coverage_plot.png"), dpi=300)

    # Statistics
    statistics = pd.DataFrame(columns = ["name", "contig_number", "longest_contig_length",
                                         "NM_percentage", "NG50", "NG90", "max_coverage"])
    i = 0
    for track in canvas.tracks:
        if track.type == 'Assembly':
            assembly = track.data
            statistics = statistics.append(pd.DataFrame(dict(name = assembly.name,
                                               contig_number = len(assembly.contigs),
                                               longest_contig_length = max([contig.length for contig in assembly.contigs]),
                                               NM_percentage = -1,
                                               NG50 = assembly.get_NG50(target_region),
                                               NG90 = assembly.get_NG90(target_region),
                                               max_coverage = max([row.covered_percentage for index, row in assembly.get_ref_coverage(target_region).iterrows()])),
                                               index = [i]))
    print(statistics)
    statistics.to_csv(os.path.join(output_dir, "statistics.csv"))
