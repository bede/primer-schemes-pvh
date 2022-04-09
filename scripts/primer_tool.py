#!/usr/bin/env python

# TODO: find a better name for this script

import argparse
from dataclasses import dataclass
from operator import attrgetter
import re
import sys

class FastaFile:
    # simple FASTA format parser included to avoid dependency on Biopython
    def __init__(self, filename):
        self.filename = filename
        self.id = None
        self.sequence = None

    def read(self):
        for line in open(self.filename):
            if line.startswith('>'):
                if self.id is not None:
                    yield self
                self.id = line.lstrip('>').rstrip().split()[0]
                self.sequence = ''
            else:
                self.sequence += line.strip()

@dataclass
class Primer:
    ref: str
    start: int
    end: int
    amplicon: str
    pool: str
    strand: str
    primer_str: str


def by_amplicon(primer):
    ampl_match = re.match(r'SARS-CoV-2_(\d+)_(LEFT|RIGHT)', primer.amplicon)
    ampl_number = float(ampl_match.group(1))
    side = ampl_match.group(2)
    if side == 'RIGHT':
        ampl_number += 0.5
    return ampl_number

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Utility script to generate different BED files, starting from ARTIC format .primer.bed and optionally .reference.fasta'
    )
    parser.add_argument('--output_format', default='bed', choices=['bed', 'insert', 'scheme', 'primers'])
    parser.add_argument('--reference_fasta', type=argparse.FileType())
    parser.add_argument('primer_file', type=argparse.FileType())
    args = parser.parse_args()

    primers = []
    sequence_id = None
    sequence = None
    if args.reference_fasta is not None:
        fasta_file = FastaFile(args.reference_fasta.name)
        fasta_file.read()
        sequence_id = fasta_file.id
        sequence = fasta_file.sequence

    for line in args.primer_file:
        parts = line.strip().split('\t')
        (ref, start_str, end_str, amplicon, pool, strand) = parts[:6]
        primer_str = None
        if len(parts) == 7:
            primer_str = parts[6]
            if sequence is not None:
                # This is here to test to see if the primers in the .primer.bed
                # match the sequence in the reference fasta file. In practice
                # you wouldn't use --reference_fasta with primer BED files that
                # contained the sequence already
                assert primer_str == sequence[start:end]
        elif len(parts) == 6 and ref == sequence_id and sequence is not None:
            start = int(start_str)
            end = int(end_str)
            primer_str = sequence[start:end]
        primer = Primer(ref, int(start_str), int(end_str), amplicon, pool, strand, primer_str)
        primers.append(primer)

    inserts = []
    insert_num = 1
    previous_primer = None
    for i, p in enumerate(sorted(primers, key=by_amplicon)):
        if args.output_format == 'bed' or args.output_format == 'scheme':
            print(p.ref, p.start, p.end, p.amplicon, end='\t', sep='\t')
            if args.output_format == 'bed':
                ampl_match = re.match(r'([^_]+)_(\d+)_(LEFT|RIGHT)(_alt\d+)?$', p.amplicon)
                if ampl_match is None:
                    exit(f"Failed to match against {p.amplicon}")
                if int(ampl_match.group(2)) % 2 != 0:
                    pool_name = 'SARS-CoV-2_1'
                else:
                    pool_name = 'SARS-CoV-2_2'
                print(pool_name, p.strand, sep='\t')
            else:
                print(p.pool, p.strand, sep='\t')
        elif args.output_format == 'insert':
            if previous_primer is not None and (i % 2) == 1:
                # skip the first primer - we are looking for inserts
                insert_start = previous_primer.end
                insert_end = p.start
                name = f'SARS-CoV-2_INSERT_{insert_num}'
                inserts.append(
                    Primer(p.ref, insert_start, insert_end, name, p.pool, p.strand, '')
                )
                insert_num += 1
            previous_primer = p
        elif args.output_format == 'primers':
            print(p.amplicon, p.start, p.end, p.pool, p.primer_str, sep='\t')
    if args.output_format == 'insert':
        for i, ins in enumerate(inserts):
            print(ins.ref, ins.start, ins.end, ins.amplicon, i % 2 + 1, '+', sep='\t')
