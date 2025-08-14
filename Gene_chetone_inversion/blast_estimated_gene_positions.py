import argparse
import subprocess
from Bio import SeqIO

# Argument parser to get file paths from the command line
parser = argparse.ArgumentParser(description="Extract sequences and perform tblastn.")
parser.add_argument("--gene_size", required=True, help="Path to the gene_size.txt file")
parser.add_argument("--fasta", required=True, help="Path to the genes FASTA file")
parser.add_argument("--db", required=True, help="BLAST database name")
parser.add_argument("--output_file", required=True, help="Path to output the results to")

args = parser.parse_args()

# Read gene names from gene_size.txt
with open(args.gene_size) as f:
    gene_names = set(line.split("\t")[0].strip() for line in f)
    gene_sizes = {}
    f.seek(0)  # Revenir au début du fichier pour lire à nouveau
    for line in f:
        fields = line.split("\t")
        gene_name = fields[0].strip()  # Nom du gène
        gene_size = int(fields[1].strip())  # Taille du gène
        gene_sizes[gene_name] = gene_size  # Ajouter au dictionnaire

# Open the output file in write mode
with open(args.output_file, 'w') as output_file:
    # Run tblastn for each gene sequence and extract the best hit
    for record in SeqIO.parse(args.fasta, "fasta"):
        for gene in gene_names:
            if gene in record.id: 
                print(f"Processing {record.id}...")

                # Run tblastn directly on the sequence
                blast_command = [
                "tblastn",
                "-db", args.db,  # BLAST database
                "-query", "-",  # Use stdin to pass sequence
                "-outfmt", "7",  # Tabular output for easy parsing
                "-word_size", "5",
                "-num_threads", "4",
            ]

                # Run BLAST and capture output
                process = subprocess.run(
                blast_command,
                input=f">{record.id}\n{record.seq}\n",  # Pass the sequence directly to stdin
                text=True,
                capture_output=True,
                )

                # Process BLAST results
                best_hit = None
                for line in process.stdout.splitlines():
                    if line.startswith("#"):
                        continue
                    fields = line.split("\t")
                    subject = fields[1]  # The second field is the subject hit (the database name)

                    # Only keep hits with the relevant database name (e.g., "ctg001860_1_np1212")
                    if "ctg001860_1_np1212" in subject:
                        if not best_hit or float(fields[11]) > float(best_hit[11]):  # Compare by bit score (field 11)
                            best_hit = fields
                            start_position = int(fields[8])  # Colonne 9 (start)
                            end_position = int(fields[9])    # Colonne 10 (end)
                            mid_position = (start_position + end_position) / 2

                            gene_size = gene_sizes[gene]
                            half_gene_size = gene_size / 2

                            # Calculer le début et la fin du gène fictif
                            gene_start = mid_position - half_gene_size
                            gene_end = mid_position + half_gene_size

                            # Affichage des résultats
                            print(f"Gene: {gene}")
                            print(f"Gene size: {gene_size}")
                            print(f"Middle position: {mid_position}")
                            print(f"Adjusted start: {gene_start}")
                            print(f"Adjusted end: {gene_end}")
                            
                            # Write to the output file
                            output_file.write(f"{gene}\t{gene_start}\t{gene_end}\n")
