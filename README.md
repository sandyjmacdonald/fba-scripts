# FBA Workshop Materials

Python scripts, files, and materials from the 3rd year computational modelling
and microbial metabolism FBA workshop.

## Filtering an existing FBA model by sequence similarity as the basis for a new one

These slides - in macOS Keynote, PDF, and Powerpoint format - walk through a 
basic example of how to take an existing flux balance analysis (FBA) metabolic
model (in .json format in the example) and filter it down based on BLAST hits 
to protein sequences from another genome. This filtered down model could be
the starting point for construction of a FBA model for the other genome.

In the example, it matches against each reaction's gene-reaction rule using 
gene locus tags, and assumes that the BLAST is a reciprocal best hit BLAST that 
produces a tab-separated table of query locus tags and hit locus tags.

Since many NCBI RefSeq genomes have the accession number as the fasta sequence 
ID, and Galaxy produces list only the fasta sequence IDs in their reciprocal 
BLAST results table, an `add_locus_tags.py` Python script is included, that can 
map the accession numbers to locus tags, given the Galaxy reciprocal BLAST 
results table, and the two fasta files. it assumes that the fasta sequence
headers have locus tags like so: `[locus_tag=b0001]`.

## The `add_locus_tags.py` Python script

This script takes a tab-separated table of reciprocal BLAST hits from Galaxy, 
with NCBI accession numbers as the query and hit IDS, and the two fasta files
used for the BLAST, and produces a tab-separated output file with two columns:
one with the query locus tags (one per line) and the other with the hits (one
or more per line, separated by semicolons). The output file can be used as the
locus tag ID file to filter an existing SBML model with locus tags from the 
query genome.

The script can be run as follows:

`python3 add_locus_tags.py -f a.fasta b.fasta -r blast_results.txt -o locus_tags_map.txt`

## The `filter_fba_model.py` Python script

This script takes an existing FBA model, in json format, and a tab-separated 
text file with locus tags from the existing model in the first column and 
corresponding locus tags from another genome in the second column. There can
be more than one locus tag in the second column, separated by semicolons.

The script considers each reaction's gene-reaction rule and keeps the reaction 
in the filtered model if sufficient genes match in the locus tag mapping file 
provided to the script. The gene-reaction rule's `and` and `or` conditions are
considered to ensure that reactions are kept/discarded correctly.

The end result is a stripped down model that will likely have gaps and will
require careful manual checking of those gaps, but it should provide a good
starting point for the new model.

You'll need to install the 
[COBRApy Python library](https://github.com/opencobra/cobrapy) before running 
the script.

The script can be run as follows:

`python3 filter_fba_model.py -l locus_tags_map.txt -m model.json`

A summary of the numbers and percentages of genes and reactions retained, and 
lists of the missing genes and reactions, is printed to the terminal also.

