#!/usr/bin/env python


import sys
import re
from argparse import ArgumentParser

parser = ArgumentParser(description='Classify a sequence as DNA or RNA')
parser.add_argument("-s", "--seq", type=str, required=True, help="Input sequence")

args = parser.parse_args()

args.seq = args.seq.upper()

length = len(args.seq)
g_content = args.seq.count("G") / length * 100
c_content = args.seq.count("C") / length * 100
a_content = args.seq.count("A") / length * 100

# To identify if the given sequence is DNA or RNA:
if re.search("^[ACGTU]+$", args.seq):
    if re.search("^[^U]+$", args.seq) and re.search("T", args.seq):
        # If the sequence does not contain U, but has T, it is DNA for sure
        t_content = args.seq.count("T") / length * 100
        result = [
            ("G content", g_content),
            ("C content", c_content),
            ("A content", a_content),
            ("T content", t_content)
        ]
        print("Results for DNA:")
        for item in result:
            print(item[0] + ": " + str(item[1]))
    elif re.search("^[^T]+$", args.seq) and re.search("U", args.seq):
        # If the sequence does not contain T, but it contains U, for sure it is RNA
        u_content = args.seq.count("U") / length * 100
        result = [
            ("G content", g_content),
            ("C content", c_content),
            ("A content", a_content),
            ("U content", u_content)
        ]
        print("Results for RNA:")
        for item in result:
            print(item[0] + ": " + str(item[1]))
    else:
        # If the sequence contains only A, C, G, U, but neither T nor U, it's neither DNA nor RNA
        print("Not DNA or RNA")
else:
    # If the sequence contains characters other than A, C, G, T, U
    print("Not DNA or RNA")

