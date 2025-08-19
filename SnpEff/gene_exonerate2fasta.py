#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO
import sys

def parse_gff_for_genes(gff_file):
    """
    Parses a GFF from exonerate to extract gene intervals.
    Keeps only 'gene' features with score >= 100.
    Returns:
        genes_by_id: dict of gene_id -> (seqid, strand, start, end)
    """
    genes_by_id = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            if feature.lower() != "gene":
                continue
            try:
                score = float(score)
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            if score < 100:
                continue
            gene_id = f"{seqid}_{strand}_{start}_{end}"
            genes_by_id[gene_id] = (seqid, strand, start, end)
    return genes_by_id

def extract_gene_sequences(genes_by_id, genome_fasta, output_fasta):
    """
    Extracts full gene sequences from the genome.
    """
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    n_written = 0

    with open(output_fasta, "w") as out:
        for gene_id, (seqid, strand, start, end) in genes_by_id.items():
            if seqid not in genome:
                print(f"[WARNING] seqid {seqid} not found in genome. Skipping.", file=sys.stderr)
                continue
            seq = genome[seqid].seq[start - 1:end]
            if strand == "-":
                seq = seq.reverse_complement()
            out.write(f">{gene_id}\n{seq}\n")
            n_written += 1

    print(f"[INFO] Genes written: {n_written}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-f", "--gff", required=True, help="Exonerate GFF file")
    parser.add_argument("-o", "--output", default="genes.fasta", help="Output FASTA file")
    args = parser.parse_args()

    genes = parse_gff_for_genes(args.gff)
    if not genes:
        print("[INFO] No genes passed the score threshold. Exiting.")
        sys.exit(0)

    extract_gene_sequences(genes, args.genome, args.output)

if __name__ == "__main__":
    main()


