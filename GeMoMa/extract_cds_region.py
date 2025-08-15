#!/usr/bin/env python3
import argparse                     # To parse command-line arguments
from pathlib import Path             # For handling filesystem paths
from Bio import SeqIO                # For reading FASTA sequences
from Bio.Seq import Seq              # For sequence manipulation

# Function to parse the attributes column of a GFF3 line
def parse_attrs(attr_str):
    d = {}
    for part in attr_str.strip().split(';'):    # Split attributes by ';'
        if not part:
            continue
        if '=' in part:                         # Split key=value
            k, v = part.split('=', 1)
            d[k] = v
    return d

def main():
    # ---------------------------
    # Command-line arguments
    # ---------------------------
    ap = argparse.ArgumentParser(
        description="Extract each CDS in a GFF region as its own AA FASTA."
    )
    ap.add_argument("--gff", required=True, help="GFF3 file")
    ap.add_argument("--fasta", required=True, help="Reference genome (FASTA)")
    ap.add_argument("--contig", required=True, help="Contig/seqid (e.g., ctg001860)")
    ap.add_argument("--start", type=int, required=True, help="Region start (1-based inclusive)")
    ap.add_argument("--end", type=int, required=True, help="Region end (1-based inclusive)")
    ap.add_argument("--outdir", required=True, help="Directory to write per-CDS FASTA files")
    args = ap.parse_args()

    # ---------------------------
    # Create output directory
    # ---------------------------
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # ---------------------------
    # Load reference genome
    # ---------------------------
    genome = {rec.id: rec.seq for rec in SeqIO.parse(args.fasta, "fasta")}
    if args.contig not in genome:
        raise SystemExit(f"Contig {args.contig} not found in FASTA")

    region_start, region_end = args.start, args.end

    # ---------------------------
    # Process GFF3
    # ---------------------------
    with open(args.gff) as fh:
        for line in fh:
            if not line or line.startswith("#"):   # Skip comments and empty lines
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:                     # Malformed GFF line
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype != "CDS":                     # We only want CDS features
                continue
            if seqid != args.contig:               # Wrong contig
                continue
            start = int(start)
            end = int(end)
            if start < region_start or end > region_end:   # Outside our region
                continue
            attrd = parse_attrs(attrs)

            # Phase adjustment — remove extra nucleotides to preserve frame
            if phase == ".":
                ph = 0
            else:
                try:
                    ph = int(phase)
                except ValueError:
                    ph = 0

            # Extract nucleotide sequence for this CDS
            dna = genome[args.contig][start-1:end]  # Convert 1-based to 0-based indexing
            if strand == "+":
                if ph > 0:
                    dna = dna[ph:]                  # Trim from 5' end for + strand
            else:
                if ph > 0:
                    dna = dna[:-ph] if ph < len(dna) else dna[:0]  # Trim from 3' end
                dna = dna.reverse_complement()     # Reverse complement for - strand

            # Ensure length is a multiple of 3 before translation
            rem = len(dna) % 3
            if rem != 0:
                dna = dna[:len(dna) - rem]

            if len(dna) == 0:                       # Skip empty sequences
                continue

            # Translate to amino acids (no stop trimming — keep *)
            aa = dna.translate(table=1, to_stop=False)

            # Build header as contig:start-end
            header = f"{args.contig}:{start}-{end}"

            # Write to individual FASTA file
            out_path = outdir / f"CDS_{args.contig}_{start}-{end}.fa"
            with open(out_path, "w") as outfh:
                outfh.write(f">{header}\n")
                for i in range(0, len(aa), 60):     # Wrap sequence to 60 chars
                    outfh.write(str(aa)[i:i+60] + "\n")

if __name__ == "__main__":
    main()
