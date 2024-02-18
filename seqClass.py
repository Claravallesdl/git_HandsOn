#!/usr/bin/env python

import sys
import re
from argparse import ArgumentParser

def classify_sequence(sequence):
    sequence = sequence.upper()
    if re.search('^[ACGT]+$', sequence):
        return 'DNA'
    elif re.search('^[ACGU]+$', sequence):
        return 'RNA'
    else:
        return 'not DNA nor RNA'

parser = ArgumentParser(description='Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type=str, required=True, help="Input sequence")
parser.add_argument("-m", "--motif", type=str, required=False, help="Motif")

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

sequence_type = classify_sequence(args.seq)

print(f'The sequence is {sequence_type}')

if args.motif:
    motif = args.motif.upper()
    print(f'Motif search enabled: looking for motif "{motif}" in sequence "{args.seq}"... ', end='')
    if re.search(motif, args.seq.upper()):
        print("Found")
    else:
        print("NOT FOUND")

