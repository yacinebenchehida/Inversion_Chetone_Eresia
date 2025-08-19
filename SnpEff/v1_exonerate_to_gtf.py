#!/usr/bin/env python3

import sys

def parse_exonerate_to_gtf(input_path, output_path):
    # Open input and output files
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        inside_gff_block = False  # Flag to track if inside GFF dump block
        gene_count = 0            # Counter to assign unique gene IDs
        transcript_count = 0      # Counter for transcripts
        current_gene_id = None    # Stores current gene ID
        current_transcript_id = None  # Stores current transcript ID

        # Write GTF header (optional but recommended)
        outfile.write('##gtf-version 2\n')

        # Iterate through each line of input
        for line in infile:
            line = line.strip()  # Remove trailing whitespace

            # Detect start and end of GFF dump block from Exonerate output
            if line.startswith("# --- START OF GFF DUMP ---"):
                inside_gff_block = True
                continue
            elif line.startswith("# --- END OF GFF DUMP ---"):
                inside_gff_block = False
                continue

            # Skip lines outside the GFF dump block or comment lines
            if not inside_gff_block or line.startswith("#"):
                continue

            fields = line.split('\t')
            # Skip malformed lines that don't have 9 fields
            if len(fields) != 9:
                continue

            # Unpack fields according to GFF format
            seqname, source, feature, start, end, score, strand, frame, attributes = fields

            # When encountering a gene feature, increment gene counter and create new gene and transcript IDs
            if feature == "gene":
                gene_count += 1
                transcript_count = 1  # Reset transcript count for new gene
                current_gene_id = f"gene{gene_count}"
                current_transcript_id = f"{current_gene_id}.t{transcript_count}"

                # Write gene feature line in GTF format with required attributes
                attr_str = f'gene_id "{current_gene_id}";'
                gtf_line = '\t'.join([
                    seqname, source, feature, start, end, score, strand, frame, attr_str
                ])
                outfile.write(gtf_line + '\n')

                # Also write transcript feature line (recommended for snpEff)
                attr_str = f'gene_id "{current_gene_id}"; transcript_id "{current_transcript_id}";'
                gtf_line = '\t'.join([
                    seqname, source, "transcript", start, end, score, strand, frame, attr_str
                ])
                outfile.write(gtf_line + '\n')

                continue  # move to next line

            # For exon or CDS features, write with current gene and transcript IDs
            if feature.lower() in ("exon", "cds"):
                # If no current gene, skip (malformed input)
                if current_gene_id is None:
                    continue

                # For snpEff, CDS feature should be uppercase "CDS"
                feature_type = "CDS" if feature.lower() == "cds" else "exon"
                attr_str = f'gene_id "{current_gene_id}"; transcript_id "{current_transcript_id}";'

                # Write feature line in GTF format
                gtf_line = '\t'.join([
                    seqname, source, feature_type, start, end, score, strand, frame, attr_str
                ])
                outfile.write(gtf_line + '\n')

if __name__ == "__main__":
    # Check for correct number of command line arguments
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: python exonerate_to_gtf.py <input_exonerate_gff> <output_gtf>\n")
        sys.exit(1)

    # Call parser function with input and output paths from command line
    parse_exonerate_to_gtf(sys.argv[1], sys.argv[2])
