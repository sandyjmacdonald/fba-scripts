#!/usr/bin/env python3

try:
    # Do this to avoid the annoying "Scaling..." message from glpk solver
    import optlang
    optlang.glpk_interface.Configuration()

    import cobra.io
    from cobra import Model, Reaction, Metabolite
except ImportError:
    print("You need to install the cobrapy library with: pip3 install cobra")

import argparse  # Handle command line arguments

# Set up command line arguments
parser = argparse.ArgumentParser(description="Create a filtered FBA model in json format from an existing model and a file of known locus tag matches.")
parser.add_argument("-l", metavar="locus_tags_map.txt", type=str, help="tab-separated mapping of known model locus tags to unknowns")
parser.add_argument("-m", type=str, metavar="model.json", help="FBA model in json format")

# Parse arguments and assign to variables
args = parser.parse_args()
tags_file = args.l  # locus tags map filename
model_file = args.m  # json model file

# A dictionary that maps locus tags in our known model to locus tags that will
# be in our new model
locus_tags = {}

# Read through the lines and construct the dictionary
with open(tags_file, "r") as f:
    for l in f.readlines():
        l = l.rstrip().split()
        a = l[0].split(";")
        b = l[1].split(";")
        for x in a:
            locus_tags[x] = b

# Make a list of the known genes' locus tags
known_genes = [str(k) for k in locus_tags.keys()]

# Read in the "known" model
model = cobra.io.load_json_model(model_file)

# Create a new model to add our matching reactions to
new_model = Model(model.id + "_filtered")

num_old_rxns = 0
num_new_rxns = 0

missing_rxns = []
matched_genes = []

for rxn in model.reactions:
    # Get the existing attributes
    id = rxn.id
    name = rxn.name
    ss = rxn.subsystem
    lb = rxn.lower_bound
    ub = rxn.upper_bound
    metabs = rxn.metabolites

    num_old_rxns += 1

    # Copy across all of the existing attributes to a new reaction, apart from
    # the GRR, which we'll check in detail
    new_rxn = Reaction(id)
    new_rxn.name = name
    new_rxn.subsystem = ss
    new_rxn.lower_bound = lb
    new_rxn.upper_bound = ub
    new_rxn.add_metabolites(metabs)

    # Set this to False by default and change to True if a match is found
    # (or no genes are associated with the reaction)
    add_reaction = False

    # Current GRR that we'll check through
    grr = rxn.gene_reaction_rule

    # Simple case where GRR is a single gene or blank
    if not " or " in grr and not " and " in grr:
        if grr == "":
            new_rxn.gene_reaction_rule = ""
            add_reaction = True
        else:
            g = grr.strip()
            # Check whether single gene in our known_genes list and, if it is,
            # then make a new GRR with the new genes
            if g in known_genes:
                new_grr = " or ".join(locus_tags[g])
                new_rxn.gene_reaction_rule = new_grr
                add_reaction = True
                matched_genes.append(str(g))

    # More than one gene, but if one or more match then we can add the reaction
    if " or " in grr and not " and " in grr:
        genes = grr.split(" or ")
        # If any one of the genes are in the known_genes list
        if any(g in known_genes for g in genes):
            # Empty list to add the new gene locus tags to
            new_genes = []
            # Go through genes and look for them in the known_genes list
            for g in genes:
                if g in known_genes:
                    new_gene = locus_tags[g]
                    new_genes += new_gene
                    matched_genes.append(str(g))
            # Construct the new GRR, joining with ors
            new_grr = " or ".join(new_genes)
            new_rxn.gene_reaction_rule = new_grr
            add_reaction = True

    # GRR is several genes with only ands, i.e. all must match
    if not " or " in grr and " and " in grr:
        genes = grr.split(" and ")
        # If all of the genes are in the known_genes list
        if all(g in known_genes for g in genes):
            new_genes = []
            # Go through genes and look for them in the known_genes list
            for g in genes:
                if g in known_genes:
                    new_gene = locus_tags[g]
                    new_genes += new_gene
                    matched_genes.append(str(g))
            # Construct the new GRR, joining with ors
            new_grr = " and ".join(new_genes)
            new_rxn.gene_reaction_rule = new_grr
            add_reaction = True

    # Complex cases with a mix of ands and ors
    if " or " in grr and " and " in grr:
        # Use this match variable to keep track of whether we have a valid match
        match = False
        # If there are a bunch of geneA and geneB ... joined by ors, split them
        alternatives = grr.split(" or ")
        # Remove brackets for each set of and genes
        alternatives = [a.replace("(", "").replace(")", "").split(" and ") for a in alternatives]
        # New list to add our new alternatives with the new locus tags
        new_alternatives = []

        for genes in alternatives:
            # We must match all of the genes
            if all(g in known_genes for g in genes):
                match = True
                new_genes = []
                # Go through genes and look for them in the known_genes list
                for g in genes:
                    if g in known_genes:
                        new_gene = locus_tags[g]
                        new_genes += new_gene
                        matched_genes.append(str(g))
                # If we have more than one gene, we need to join the new 
                # locus tags with ands and wrap in round brackets
                if len(new_genes) > 1:
                    new_genes = "(" + " and ".join(new_genes) + ")"
                # If there's only one, then just add it
                else:
                    new_genes = new_genes[0]
                new_alternatives.append(new_genes)
        # Join the sets of ands with ors
        new_grr = " or ".join(new_alternatives)
        new_rxn.gene_reaction_rule = new_grr
        if match:
            add_reaction = True

    if add_reaction:
        new_model.add_reactions([new_rxn])
        num_new_rxns += 1
    else:
        missing_rxns.append(id)

# Summary stats of our model mapping
original_genes = set([gene.id for gene in model.genes])
num_original_genes = len(original_genes)
num_matched_genes = len(set(matched_genes))

matched_genes = set(matched_genes)
missing_genes = original_genes - matched_genes

# Transfer the model objective over to the new model
obj = model.objective
new_model.objective = obj

# Write the new model a json file
new_model_file = model_file.split("/")
if len(new_model_file) > 1:
    new_model_file = "/".join(new_model_file[:-1]) + "/filtered_" + new_model_file[-1]
else:
    new_model_file = "filtered_" + new_model_file[0]
cobra.io.save_json_model(new_model, new_model_file)

# Print a summary, along with lists of missing genes and reactions
print(f"""
Model summary: 

{num_matched_genes} genes matched from {num_original_genes} genes in original model ({int(num_matched_genes/num_original_genes * 100)}%)

Missing genes:
{",".join(missing_genes)}

{num_new_rxns} reactions added from {num_old_rxns} reactions in original model ({int(num_new_rxns/num_old_rxns * 100)}%)

Missing reactions:
{",".join(missing_rxns)}

Filtered model written to {new_model_file}
""")
