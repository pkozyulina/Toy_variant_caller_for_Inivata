#!/usr/bin/env python3

from Bio import pairwise2
import re

class Read:

    def __init__(self, seq, read_id):
        self.id = read_id
        self.seq = seq
        self.ref = ""
        self.score = 0 # mapping distance score
        self.amplicon = re.split('\.', re.split('_', self.id)[1])[0]
        self.amplicon_coord = (0, 0)
        self.variants = []
        self.gene = re.split('_', self.id)[0]
        self.align_coordinates = (0, 0)
        self.read_coord = (0, 0)

    
    def align_score(self, ref):
        """
        Get the alignment score in order to deduce if the ref and the read are identical
        """
        for ref_name, reference in ref.items():

            for a in pairwise2.align.localms(reference[0], self.seq, 2, -1, -1, -0.5):
                al1, al2, score, begin, end = a

                if score > self.score:
                    self.score, self.align_coordinates = score, (begin, end)
                    self.amplicon = ref_name
                    self.amplicon_coord = reference[1]

        self.ref = ref[self.amplicon][0][self.align_coordinates[0]:self.align_coordinates[1]]
        self.read_coord = (self.align_coordinates[0] + self.amplicon_coord[0], self.align_coordinates[1] + self.amplicon_coord[0])


    def is_mutant(self):
        return self.score != len(self.seq)*2

    def mismatch_number(self):
        return len(self.seq)*2 - self.score
    
    def find_variants(self):
        """
        Find mismatch coordinate
        """
        for n in range(len(self.seq)):
            if self.seq[n] != self.ref[n]:
                self.variants.append(n)

    def mismatch(self, n):
        """
        Find mismatches
        """
        return self.ref[n], self.seq[n]

    def mismatches(self):
        return map(self.mismatch, self.variants)

    def __str__(self):
        return "%s: %s %s\n%s\n%s" % (self.score, self.amplicon, self.gene, self.seq, self.ref)

    def transfrom_coordinates(self, local_coord):
        return local_coord + self.align_coordinates[0] + self.amplicon_coord[0]

    def return_mismatches(self):
        mismatches = self.mismatches()
        return zip(mismatches, list(map(self.transfrom_coordinates, self.variants)))

