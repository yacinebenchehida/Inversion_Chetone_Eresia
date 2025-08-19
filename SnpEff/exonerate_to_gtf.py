#!/usr/bin/env python3

import argparse  # For argument parsing
from collections import defaultdict  # For default dictionary
from Bio import SeqIO  # For reading FASTA files
from Bio.Seq import Seq  # For sequence manipulation

def parse_gff(gff_file):
    """
    Parse GFF file to extract CDS coordinates grouped by gene.
    Args:
        gff_file: Path to Exonerate GFF output
    Returns:
        cds_dict: dict with gene_key -> list of (start, end, strand)
        gene_ids: dict gene_key -> gene_id string
    """
    cds_dict = defaultdict(list)  # key: (seqid, gene_start, gene_end), value: list of CDS tuples
    gene_ids = {}
    current_gene = None
    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            if len(parts) == 8:
                parts.append(".")
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            
            if feature == "gene":
                gene_key = (seqid, int(start), int(end))
                # Parse gene id from attributes reliably
                gene_id = None
                for attr in attributes.split(";"):
                    if "gene_id" in attr:
                        gene_id = attr.split("gene_id")[-1].replace("\"", "").strip()
                        break
                if not gene_id:
                    gene_id = attributes.split("sequence")[-1].strip().replace(";", "").replace(" ", "_")
                gene_ids[gene_key] = gene_id
                current_gene = gene_key
            
            elif feature == "cds" and current_gene is not None:
                cds_dict[current_gene].append((int(start), int(end), strand))
    return cds_dict, gene_ids

def trim_cds_to_protein(cds_coords, strand, genome_seq, protein_seq):
    """
    Adjust CDS boundaries so that concatenated CDS translate exactly to protein_seq.
    Trim or extend CDS ends as needed, working exon by exon.
    Args:
        cds_coords: list of (start, end, strand) tuples in genomic order
        strand: '+' or '-'
        genome_seq: reference Seq object for the entire contig/scaffold
        protein_seq: expected protein sequence as a Seq object
    Returns:
        adjusted_cds: list of (start, end) tuples adjusted for perfect translation
    """
    # Sort CDS coords by strand
    exons = sorted(cds_coords, key=lambda x: x[0], reverse=(strand == "-"))
    concatenated_seq = Seq("")
    adjusted_cds = []

    for i, (start, end, _) in enumerate(exons):
        exon_seq = genome_seq[start-1:end]  # extract exon genomic sequence (1-based coords)
        concatenated_seq += exon_seq

    # Translate concatenated CDS (as is)
    translated = concatenated_seq.translate(to_stop=False)
    
    # If translation matches protein exactly, no trimming needed
    if str(translated) == str(protein_seq):
        return [(start, end) for (start, end, _) in exons]

    # Otherwise, try trimming nucleotides at exon boundaries from the last exon progressively
    # to fix frameshift or stop codons in translation
    # We only trim at the 3' end of the concatenated CDS in genomic coordinates
    # (For negative strand, trimming happens at the start of the first exon because of reverse orientation)
    
    # Convert protein to string for comparisons
    prot_str = str(protein_seq)
    
    # Attempt trimming 0 to 2 nucleotides from the last exon end to fix frame
    for trim in range(3):
        # Build adjusted concatenated seq with trimming applied at 3' end
        if strand == "+":
            # Trim last exon end coordinate by trim nt
            new_exons = exons[:-1] + [(exons[-1][0], exons[-1][1] - trim)]
        else:
            # For '-' strand, trim at start coordinate of first exon (lowest genomic coordinate)
            new_exons = [(exons[0][0] + trim, exons[0][1])] + exons[1:]

        # Concatenate adjusted exons
        adj_seq = Seq("")
        for s, e in new_exons:
            adj_seq += genome_seq[s-1:e]
        
        # Reverse complement if negative strand
        if strand == "-":
            adj_seq = adj_seq.reverse_complement()

        # Translate and compare to protein sequence
        adj_trans = adj_seq.translate(to_stop=False)
        if str(adj_trans) == prot_str:
            # Found correct trimming
            return [(s, e) for s, e in new_exons]

    # If no trimming matches, raise error (could implement more complex handling here)
    raise ValueError("Cannot adjust CDS boundaries to match protein sequence exactly.")

def calculate_phases(cds_coords, strand):
    """
    Calculate phases for each CDS based on concatenated CDS length.
    Args:
        cds_coords: list of (start, end) tuples sorted in transcription order
        strand: '+' or '-'
    Returns:
        List of phases (0,1,2) corresponding to each CDS in order
    """
    phases = []  # list to store phases
    length_sum = 0  # cumulative CDS length counter
    for start, end in cds_coords:  # iterate over CDS coords
        phases.append(length_sum % 3)  # current phase is cumulative length mod 3
        exon_len = end - start + 1  # calculate exon length
        length_sum += exon_len  # update cumulative length
    return phases  # return list of phases

def write_snpEff_gtf(output_gtf, cds_dict, gene_ids, genome_fasta, protein_fasta):
    """
    Write a snpEff-compatible GTF file using adjusted CDS with correct phase.
    Args:
        output_gtf: output GTF filename
        cds_dict: dict gene_key -> list of (start, end, strand)
        gene_ids: dict gene_key -> gene_id string
        genome_fasta: genome fasta filename
        protein_fasta: protein fasta filename with expected protein sequences
    """
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    proteins = SeqIO.to_dict(SeqIO.parse(protein_fasta, "fasta"))
    
    with open(output_gtf, "w") as gtf_out:
        gtf_out.write("##gtf-version 2\n")  # header

        for gene_key, cds_list in cds_dict.items():
            seqid = gene_key[0]
            if seqid not in genome:
                continue
            strand = cds_list[0][2]
            
            # Get expected protein sequence from protein fasta, match by gene_id
            gene_id = gene_ids.get(gene_key, None)
            if gene_id is None or gene_id not in proteins:
                raise ValueError(f"Protein for gene {gene_id} not found in protein FASTA.")
            protein_seq = proteins[gene_id].seq
            
            # Adjust CDS coords to fit protein exactly
            adjusted_cds = trim_cds_to_protein(cds_list, strand, genome[seqid].seq, protein_seq)

            # Calculate phases for adjusted CDS
            # Sort CDS in transcription order for phase calculation
            cds_sorted = sorted(adjusted_cds, reverse=(strand == "-"))
            phases = calculate_phases(cds_sorted, strand)
            
            # Gene line (start = min CDS start, end = max CDS end)
            gene_start = min([start for start, end in cds_sorted])
            gene_end = max([end for start, end in cds_sorted])
            
            # IDs
            transcript_id = f"{gene_id}.t1"
            
            # Write gene feature
            attr_gene = f'gene_id "{gene_id}";'
            gtf_out.write(f"{seqid}\tsnpEff\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attr_gene}\n")
            
            # Write transcript feature
            attr_tr = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
            gtf_out.write(f"{seqid}\tsnpEff\ttranscript\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t{attr_tr}\n")
            
            # Write CDS features with phases
            for (start, end), phase in zip(cds_sorted, phases):
                attr_cds = f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'
                gtf_out.write(f"{seqid}\tsnpEff\tCDS\t{start}\t{end}\t.\t{strand}\t{phase}\t{attr_cds}\n")


def write_translated_proteins(cds_dict, gene_ids, genome_fasta, output_fasta):
    """
    Extract CDS sequences from genome, concatenate and translate them,
    then write translated protein sequences to a FASTA file.
    Args:
        cds_dict: dict of CDS coords keyed by gene (list of (start, end, strand))
        gene_ids: dict of gene ids keyed by gene
        genome_fasta: path to genome FASTA file
        output_fasta: path to output protein FASTA file
    """
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))  # Load genome sequences into dictionary keyed by sequence ID

    with open(output_fasta, "w") as out_fa:  # Open output FASTA file for writing translated proteins
        for gene_key in cds_dict:  # Iterate over each gene key in CDS dictionary
            seqid = gene_key[0]  # Extract scaffold/contig ID from gene key tuple

            if seqid not in genome:  # Check if scaffold exists in genome sequences; skip if not found
                continue

            strand = cds_dict[gene_key][0][2]  # Get strand information from first CDS of gene

            # Sort CDS entries by start coordinate in transcription order (reverse if negative strand)
            sorted_cds = sorted(cds_dict[gene_key], key=lambda x: x[0], reverse=(strand == "-"))

            full_cds_seq = ""  # Initialize empty string to accumulate full CDS nucleotide sequence

            for start, end, _ in sorted_cds:  # Iterate over sorted CDS coordinates
                exon_seq = genome[seqid].seq[start - 1:end]  # Extract exon nucleotide sequence (1-based inclusive coordinates)
                full_cds_seq += str(exon_seq)  # Append exon sequence as string to full CDS sequence

            if strand == "-":  # If gene is on negative strand
                full_cds_seq = str(Seq(full_cds_seq).reverse_complement())  # Reverse complement the concatenated CDS sequence

            trim_len = len(full_cds_seq) % 3  # Calculate number of nucleotides to trim to reach multiple of 3 length
            if trim_len != 0:  # If length is not a multiple of 3
                full_cds_seq = full_cds_seq[:-trim_len]  # Trim nucleotides from the end to make length divisible by 3

            protein_seq = str(Seq(full_cds_seq).translate(to_stop=False))  # Translate nucleotide CDS to protein sequence (include stop codons as '*')

            gene_id = gene_ids.get(gene_key, f"{seqid}_{gene_key[1]}_{gene_key[2]}")  # Get gene ID or fallback to formatted scaffold positions

            out_fa.write(f">{gene_id}\n{protein_seq}\n")  # Write FASTA header and protein sequence to output file



def main():
    parser = argparse.ArgumentParser(description="Generate snpEff-ready GTF and protein FASTA from Exonerate GFF and proteins")
    parser.add_argument("-gff", required=True, help="Exonerate GFF file")
    parser.add_argument("-genome", required=True, help="Genome FASTA file")
    parser.add_argument("-prot", required=True, help="Protein FASTA file matching genes in GFF")
    parser.add_argument("-o", "--output", required=True, help="Output prefix (without extension)")
    args = parser.parse_args()

    cds_dict, gene_ids = parse_gff(args.gff)  # Parse GFF to get CDS and gene info

    # Write snpEff-compatible GTF
    write_snpEff_gtf(args.output + ".gtf", cds_dict, gene_ids, args.genome, args.prot)

    # Write translated protein sequences from genome CDS for verification
    write_translated_proteins(cds_dict, gene_ids, args.genome, args.output + ".proteins.fasta")

if __name__ == "__main__":
    main()  # run main function