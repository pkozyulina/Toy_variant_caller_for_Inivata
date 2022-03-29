#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict
import pickle


# read bed file and return as list of tuples
def parse_bed(input_bed):
    content = []
    for line in input_bed:
        content.append(line.strip().split())
    return content


def reformat_bed(bed):
    if bed[2] < bed[3]:
        return tuple(bed[2:4])
    else:
        return (bed[3], bed[2])


def main():
    # parsing command line arguments
    parser = argparse.ArgumentParser(description='Run parse_annotation.py to generate annotation for coding test.')
    parser.add_argument('-i', '--input', help='Input bed files containing gene coordinates', metavar='gene.bed', nargs='+',
                        required=True)
    args = parser.parse_args()

    annotation = defaultdict(list)
    for bed_file in args.input:
        with open(bed_file, 'r') as file:
            bed = parse_bed(file)
            bed = list(map(reformat_bed, bed))
            gene_name = re.split(r"/", re.split(r'\.', bed_file)[0])[1]
            annotation[gene_name] = bed
            print(bed)

    with open('data/annotation.pkl', 'wb') as f:
        pickle.dump(annotation, f)


if __name__ == '__main__':
    main()
