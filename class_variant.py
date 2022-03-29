#!/usr/bin/env python3

from collections import defaultdict


class Variant:

    def __init__(self, ref, coord):
        self.coordinate = coord
        self.ref = ref
        self.mut = defaultdict(int)
        self.freqs = []

    def add_variant(self, mut):
        self.mut[mut] += 1

    def calc_freq(self, total_reads):
        for mut, number in self.mut.items():
            freq = number / total_reads
            self.freqs.append([mut, freq, number])
        return self.ref, self.freqs

    def to_plot_freqs(self):
        return self.freqs[0][1], 1 - self.freqs[0][1], '%s\n%s->%s' % (self.coordinate, self.ref, self.freqs[0][0])
