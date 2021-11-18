#!/usr/bin/env python3

import re  # Regular expressions to match locus_tag field
import argparse  # Handle command line arguments

# Set up command line arguments
parser = argparse.ArgumentParser(description="Create an output file with locus tags from a Galaxy RBH BLAST hits file and two fasta files.")
parser.add_argument("-f", metavar="seqs.fasta", type=str, nargs="+", help="fasta input files")
parser.add_argument("-r", type=str, metavar="results.txt", help="Galaxy RBH blast hits input file")
parser.add_argument("-o", type=str, metavar="out.txt", help="out file name")

# Parse arguments and assign to variables
args = parser.parse_args()
fastas = args.f
results = args.r
out_file = args.o

# A dictionary to hold the sequence IDs (key) and their locus tags (value)
id_dict = {}

# Parse the fasta files
for fn in fastas:
    with open(fn, "r") as f:
        for l in f.readlines():
            # Ignore any lines that don't have locus tags
            if "locus" in l:
                l = l.rstrip()
                # The ID is the first thing on the line
                id = l.split()[0][1:]
                # Horrid regular expression to match locus tag in description
                match = re.search(r"\[locus_tag\=(\S+)\]", l)
                # Get the locus tag from the match
                locus_tag = match.group(1)
                # Add the ID and locus tag to the dictionary
                id_dict[id] = locus_tag

# Open the BLAST results file
with open(results, "r") as f:
    # Open the out file to write the locus tags to
    with open(out_file, "w") as out:
        lines = f.readlines()
        # Ignore line 0 with the headers, and only loop through lines from 1
        for l in lines[1:]:
            # Split the line by spaces
            l = l.rstrip().split()
            # In theory, there can be multiple query IDs and hit IDs separated
            # by semicolons, so split them
            query_ids = l[0].split(";")
            hit_ids = l[1].split(";")
            # Get the locus tags from the hit IDs in the dictionary
            query_locus_tags = [id_dict[q] for q in query_ids]
            hit_locus_tags = [id_dict[h] for h in hit_ids]
            # Write the locus tags separated by a tab and joined by semicolons
            # if there are multiple of each
            out.write(";".join(query_locus_tags) + "\t" + ";".join(hit_locus_tags) + "\n")