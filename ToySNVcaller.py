#!/usr/bin/env python3

import argparse
import gzip
import re
import pickle
from tqdm import tqdm
from collections import defaultdict
import logging
from Bio import SeqIO
from parse_annotation import parse_bed
from class_read import Read
from class_variant import Variant
from Bio import motifs
import matplotlib.pyplot as plt
import numpy as np


log = logging.getLogger("toy-logger")


# parsing command line arguments
parser = argparse.ArgumentParser(description='Toy aligner and variant caller for the code task from Inivata. '
                                             'In order to get gene names you need to prepare annotation first. '
                                             'For that - run parse_annotation.py and provide it with bed files.'
                                             'Dependencies: biopython')
parser.add_argument('-i', '--input', help='Input fasta files - one or more files containing reads',
                    metavar='amplicon.fa.gz', nargs='+', required=True)
parser.add_argument('-r', '--reference', help='Genome reference in fasta-format',
                    metavar='genome.fa.gz', required=True)
parser.add_argument('-b', '--bed', help='Reference bed file with amplicon coordinates',
                    metavar='amplicon_coordinates.bed', type=argparse.FileType(), required=True)
parser.add_argument('-o', '--output', help='Output filename for the table with '
                                           'Gene name, Frequency, N amplicons, amplicon seq, ref seq',
                    metavar='output.tsv', type=argparse.FileType('w'), default='output.tsv')
parser.add_argument('-bp', '--barplot', help='use if you want to get a bar plot with mutations data',
                    action='store_true')
parser.add_argument('-a', '--annotation', help='use if you want to get a bar plot with mutations data',
                    action='store_true')
args = parser.parse_args()


# there is this annotation parser (parse_annotation.py), however gene1 amplicon coordinates appeared to be
# outside of exon range of annotation, so I decided to hardcode the annotation just to save time for coding
annotation = {'gene1': [('55019278', '55205617')], 'gene2': [('140734770', '140924566')]}


# need to be uncommented in main function in order to read annotation from pkl-file
def read_annotation(annot_filename='data/annotation.pkl'):
    try:
        with open(annot_filename, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        print('There is probably no file  data/annotation.pkl. '
              'Please, generate it with parse_annotation.py!\n', e)


# read fasta and return as generator: fasta sequence and fasta ID
def parse_reads(input_file):
    for ref_file in input_file:
        with gzip.open(ref_file, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                yield str(record.seq), record.id


# use bed file in order to trim reference according to amplicon coordinates
def parse_ref(input_ref, input_bed):
    amplicons = defaultdict(tuple)
    with gzip.open(input_ref, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            for chromosome, start, end, amplicon in input_bed:
                amplicon = re.sub(r'[_]', '', amplicon)
                if record.id == chromosome:
                    start = int(start) - 1
                    amplicons[amplicon] = (record.seq[start: int(end)], (start, end))
    return amplicons


# count total read numer per amplicon coordinates
def count_total_read_number(coordinate, read_coordinates):
    total_reads = 0
    for coord, number in read_coordinates.items():
        if coordinate in range(coord[0], coord[1]):
            total_reads += number
    return total_reads


# if read contains mutations then this function returns alternative variants as a dict containing class Variant objects
def return_variants(mismatches, variants):
    for mismatch in mismatches:
        if mismatch[1] in variants:
            variants[mismatch[1]].add_variant(mismatch[0][1])
        else:
            variants[mismatch[1]] = Variant(mismatch[0][0], mismatch[1])
    return variants


# returns reference amplicon and mutated amplicon sequences
def mutate_amplicon(amplicons, coordinate, mutation):
    for amplicon, seq in amplicons.items():
        if coordinate in range(int(seq[1][0]), int(seq[1][1])):
            index = coordinate - int(seq[1][0])
            mutant_seq = str(seq[0])[:index] + mutation + str(seq[0])[index + 1:]
            return str(seq[0]), mutant_seq


# get gene name from the annotation
def get_gene_name(annotation, coordinate):
    for gene, coord in annotation.items():
        for c in coord:
            if coordinate in range(int(c[0]), int(c[1])):
                return gene


# bar plot for mutations detected
def plot_freqs(number_mut, number_total, legend):
    # creating the bar plot
    ind = np.arange(len(number_mut))
    p1 = plt.bar(ind, number_total)
    p2 = plt.bar(ind, number_mut, bottom=number_total)

    plt.xticks(ind, legend)
    plt.legend((p1[0], p2[0]), ('reference', 'mutation'))
    plt.ylabel("Frequency")
    plt.title("Frequency of detected mutations")
    plt.savefig('mutations_barplot.png')


def main():

    # prepare reference data
    if args.annotation:
        annotation = read_annotation() #this function was supposed to read annotation prepared by parse_annotation.py
    bed = parse_bed(args.bed)
    amplicons = parse_ref(args.reference, bed)

    # prepare to collect variants
    variants = defaultdict(Variant)
    read_coordinates = defaultdict(int)

    print('Watch out - the process might take some time to align all reads properly!')

    # parse reads and collect variant data
    cnt = 0
    for read, read_id in tqdm(parse_reads(args.input)):
        read = Read(read, read_id)
        read.align_score(amplicons)
        cnt += 1

        read.find_variants()
        if read.is_mutant():
            mismatches = read.return_mismatches()
            variants = return_variants(mismatches, variants)
        read_coordinates[read.read_coord] += 1

    print("%s reads have been processed!" % cnt)

    # prepare for bar plot
    if args.barplot:
        total_for_plot, muts_for_plot, legend_for_plot = [], [], []

    # count total read number per bed amplicon intervals and calculate variant frequency
    for coordinate, variant in variants.items():
        total_reads = count_total_read_number(coordinate, read_coordinates)
        ref, freqs = variant.calc_freq(total_reads)

        for freq in freqs:
            gene_name = get_gene_name(annotation, coordinate)
            ref_seq, mut_seq = mutate_amplicon(amplicons, coordinate, freq[0])
            args.output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene_name, ref, freq[0], freq[1], freq[2], ref_seq,
                                                                mut_seq))

        print("Variant %s in coordinate %s for %s has been written to %s file." % (freq[0], coordinate, gene_name,
                                                                                   args.output.name))

        # prepare data for bar plot
        if args.barplot:
            m, n, v = variant.to_plot_freqs()
            muts_for_plot.append(m)
            legend_for_plot.append(v)
            total_for_plot.append(n)

    # draw bar plot
    if args.barplot:
        plot_freqs(muts_for_plot, total_for_plot, legend_for_plot)

    print('All is DONE!')


if __name__ == '__main__':
    main()
