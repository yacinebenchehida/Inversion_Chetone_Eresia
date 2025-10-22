#!/usr/bin/env python3

# Import required libraries
import argparse  # for parsing command-line arguments
from Bio import SeqIO  # for reading and writing FASTA files
from collections import defaultdict  # for storing sequences per gene
import random  # for random selection

def parse_gff(gff_file):
    """
    Parse the GFF file and extract ordered gene IDs per scaffold.
    Return a dictionary mapping scaffold -> ordered list of gene IDs.
    """
    scaffold_genes = defaultdict(list)  # initialize dictionary of lists
    with open(gff_file, 'r') as f:  # open GFF file
        for line in f:  # iterate through each line
            if line.startswith("#"):  # skip comment lines
                continue
            parts = line.strip().split("\t")  # split line into columns
            if len(parts) < 9:  # skip malformed lines
                continue
            feature_type = parts[2]  # third column = feature type
            if feature_type != "gene":  # process only "gene" entries
                continue
            scaffold = parts[0]  # first column = scaffold ID
            attributes = parts[8]  # ninth column = attributes
            gene_id = None  # initialize variable for gene ID
            for attr in attributes.split(";"):  # iterate through attributes
                if attr.strip().startswith("ID="):  # look for "ID="
                    gene_id = attr.strip().replace("ID=", "")  # extract gene ID
                    break  # stop after finding ID
            if gene_id:  # if a gene ID was found
                scaffold_genes[scaffold].append(gene_id)  # store it in order
    return scaffold_genes  # return the dictionary

def select_contiguous_sets(scaffold_genes, n_sets, n_genes):
    """
    Randomly select contiguous sets of genes.
    Return a list of lists of gene IDs (one per set).
    """
    remaining = {scaffold: genes[:] for scaffold, genes in scaffold_genes.items()}  # copy data
    selected_sets = []  # list to store gene sets
    eligible_scaffolds = [s for s, genes in remaining.items() if len(genes) >= n_genes]  # scaffolds with enough genes

    while len(selected_sets) < n_sets:  # repeat until enough sets are built
        if not eligible_scaffolds:  # if no scaffold remains eligible
            raise ValueError("Not enough genes to select requested sets.")
        scaffold = random.choice(eligible_scaffolds)  # randomly select a scaffold
        genes = remaining[scaffold]  # retrieve list of genes for that scaffold
        max_start = len(genes) - n_genes  # last possible starting index
        start_idx = random.randint(0, max_start)  # randomly pick a start index
        gene_set = genes[start_idx:start_idx + n_genes]  # slice contiguous genes
        selected_sets.append(gene_set)  # store this set
        del genes[start_idx:start_idx + n_genes]  # remove used genes
        if len(genes) < n_genes:  # if scaffold now too short, remove it
            eligible_scaffolds.remove(scaffold)
    return selected_sets  # return all selected sets

def load_sequences(fasta_path):
    """
    Load all sequences from the reference FASTA into a dictionary keyed by gene ID substring.
    Return a list of SeqRecord objects for later searching.
    """
    seq_records = list(SeqIO.parse(fasta_path, "fasta"))  # read all records
    return seq_records  # return list for substring search

def extract_sequences_for_sets(selected_sets, seq_records):
    """
    For each set, find sequences in seq_records whose header contains the gene ID.
    If multiple matches exist, keep only the longest one.
    Return a list of lists of SeqRecord objects, one sublist per set.
    """
    all_set_records = []  # store sequences for each set
    for gene_set in selected_sets:  # iterate through all sets
        set_records = []  # store sequences for this set
        for gene_id in gene_set:  # iterate through genes in the set
            matches = [rec for rec in seq_records if gene_id in rec.id]  # partial match search
            if not matches:  # if no match found
                continue  # skip this gene
            longest = max(matches, key=lambda x: len(x.seq))  # keep the longest sequence
            set_records.append(longest)  # store it
        all_set_records.append(set_records)  # store this set’s sequences
    return all_set_records  # return list of sets

def write_sets_to_fasta(all_set_records):
    """
    Write each gene set to its own FASTA file (set_1.fa, set_2.fa, ...).
    """
    for idx, records in enumerate(all_set_records, 1):  # iterate through sets
        output_file = f"set_{idx}.fa"  # define output filename
        SeqIO.write(records, output_file, "fasta")  # write sequences in FASTA format

def main():
    """Main function to parse arguments and run all steps."""
    parser = argparse.ArgumentParser(description="Select random contiguous gene sets and extract their sequences.")
    parser.add_argument("--gff", required=True, help="Input GFF file")
    parser.add_argument("--n_sets", type=int, required=True, help="Number of sets to select")
    parser.add_argument("--n_genes", type=int, required=True, help="Number of genes per set")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file containing CDS or protein sequences")
    args = parser.parse_args()  # parse input arguments

    scaffold_genes = parse_gff(args.gff)  # parse the GFF to get gene order
    selected_sets = select_contiguous_sets(scaffold_genes, args.n_sets, args.n_genes)  # build random sets
    seq_records = load_sequences(args.fasta)  # load reference FASTA
    all_set_records = extract_sequences_for_sets(selected_sets, seq_records)  # find sequences for each set
    write_sets_to_fasta(all_set_records)  # write each set to its own FASTA

if __name__ == "__main__":
    main()  # execute script
