#!/usr/bin/env python3  # Use Python 3 interpreter

import argparse  # Parse command-line arguments
from collections import defaultdict  # Kept to mirror prior script structure
from pathlib import Path  # Filesystem path utilities
from Bio import SeqIO  # FASTA I/O
from Bio.Seq import Seq  # Sequence operations (reverse_complement)

def parse_attrs(attr_str):  # Parse GFF attribute column into a dict
    d = {}  # Init dict for attributes
    for part in attr_str.strip().split(';'):  # Split attributes by semicolon
        if not part:  # Skip empty parts
            continue  # Continue loop
        if '=' in part:  # Only handle key=value pairs
            k, v = part.split('=', 1)  # Split into key and value
            d[k] = v  # Store in dict
    return d  # Return parsed attributes

def main():  # Entry point
    ap = argparse.ArgumentParser(  # Build CLI parser
        description="Extract gene DNA sequences from a GFF region; write one FASTA per gene."  # Help text
    )  # End parser build
    ap.add_argument("--gff", required=True, help="GFF3 file")  # GFF input path
    ap.add_argument("--fasta", required=True, help="Reference genome (FASTA)")  # FASTA input path
    ap.add_argument("--contig", required=True, help="Contig/seqid (e.g., ctg001860)")  # Target contig
    ap.add_argument("--start", type=int, required=True, help="Region start (1-based inclusive)")  # Region start
    ap.add_argument("--end", type=int, required=True, help="Region end (1-based inclusive)")  # Region end
    ap.add_argument("--outdir", required=True, help="Directory to write per-gene FASTA files")  # Output dir
    args = ap.parse_args()  # Parse CLI args

    outdir = Path(args.outdir)  # Convert outdir to Path
    outdir.mkdir(parents=True, exist_ok=True)  # Create outdir if needed

    genome = {rec.id: rec.seq for rec in SeqIO.parse(args.fasta, "fasta")}  # Load genome sequences into dict
    if args.contig not in genome:  # Ensure contig exists
        raise SystemExit(f"Contig {args.contig} not found in FASTA")  # Abort if missing

    region_start, region_end = args.start, args.end  # Assign region bounds

    with open(args.gff) as fh:  # Open GFF file
        for line in fh:  # Iterate lines
            if not line or line.startswith("#"):  # Skip blanks/comments
                continue  # Next line
            parts = line.rstrip("\n").split("\t")  # Split GFF columns
            if len(parts) < 9:  # Require 9 columns
                continue  # Skip malformed lines
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts  # Unpack columns
            if ftype != "gene":  # Keep only gene features
                continue  # Skip non-gene
            if seqid != args.contig:  # Enforce contig filter
                continue  # Skip other contigs
            start = int(start)  # Convert start to int
            end = int(end)  # Convert end to int
            if start < region_start or end > region_end:  # Enforce region bounds (fully inside)
                continue  # Skip out-of-region genes

            _ = parse_attrs(attrs)  # Parse attributes (kept for symmetry; not required below)

            dna = genome[args.contig][start - 1:end]  # Extract genomic slice (1-based to 0-based)
            if strand == "-":  # If on negative strand
                dna = dna.reverse_complement()  # Reverse-complement to 5'->3' orientation

            header = f"{args.contig}:{start}-{end}"  # Build header as coordinates
            out_path = outdir / f"Gene_{args.contig}_{start}-{end}.fa"  # Per-gene FASTA path
            with open(out_path, "w") as outfh:  # Open output FASTA
                outfh.write(f">{header}\n")  # Write FASTA header
                seq = str(dna)  # Convert sequence to string
                for i in range(0, len(seq), 60):  # Wrap sequence at 60 chars
                    outfh.write(seq[i:i + 60] + "\n")  # Write wrapped lines

if __name__ == "__main__":  # Standard Python entry guard
    main()  # Invoke main

