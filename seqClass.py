#!/usr/bin/env python

<<<<<<< HEAD
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

=======
# Import necessary modules
import sys, re
from argparse import ArgumentParser

# Create an argument parser with options
>>>>>>> 5ff4be76e4a8913780e017be0d3a08d93add67c1
parser = ArgumentParser(description='Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type=str, required=True, help="Input sequence")
parser.add_argument("-m", "--motif", type=str, required=False, help="Motif")

# Check if there are any command-line arguments. If not, print the help message and exit.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

# Parse the command-line arguments
args = parser.parse_args()

<<<<<<< HEAD
sequence_type = classify_sequence(args.seq)

print(f'The sequence is {sequence_type}')
=======
# Convert the sequence to uppercase for case-insensitive comparison
args.seq = args.seq.upper()
>>>>>>> 5ff4be76e4a8913780e017be0d3a08d93add67c1

# Check if the sequence consists of only valid DNA or RNA characters
if re.search('^[ACGTU]+$', args.seq):
    # If 'T' is present, it's DNA; if 'U' is present, it's RNA
    if re.search('T', args.seq):
        print('The sequence is DNA')
    elif re.search('U', args.seq):
        print('The sequence is RNA')
    else:
        print('The sequence can be DNA or RNA')
else:
    print('The sequence is not DNA nor RNA')

# Check if a motif is provided
if args.motif:
<<<<<<< HEAD
    motif = args.motif.upper()
    print(f'Motif search enabled: looking for motif "{motif}" in sequence "{args.seq}"... ', end='')
    if re.search(motif, args.seq.upper()):
        print("Found")
    else:
        print("NOT FOUND")

=======
    # Convert the motif to uppercase for case-insensitive comparison
    args.motif = args.motif.upper()
    print(f'Motif search enabled: looking for motif "{args.motif}" in sequence "{args.seq}"... ', end='')

    # Check if the motif is present in the sequence
    if re.search(args.motif, args.seq):
        print("Found")
    else:
        print("NOT FOUND")
>>>>>>> 5ff4be76e4a8913780e017be0d3a08d93add67c1
