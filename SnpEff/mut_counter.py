import vcfpy  # Import vcfpy library to parse VCF files
import gffutils  # Import gffutils to parse GFF files and compute CDS lengths
import os         # File operations
import argparse   # To parse command-line arguments


def count_mutation_type(vcf_file, impact_type):
    """
    Count variants in a VCF with a specific impact type from SnpEff annotation.

    Parameters:
        vcf_file (str): path to the VCF file
        impact_type (str): 'HIGH', 'MODERATE', 'LOW', 'MODIFIER'

    Returns:
        int: number of variants matching the impact type
    """
    count = 0  # Initialize counter for variants of specified impact type
    reader = vcfpy.Reader.from_path(vcf_file)  # Open VCF file for reading
    for record in reader:  # Iterate over each variant record in VCF
        if 'ANN' in record.INFO:  # Check if SnpEff annotation exists in INFO field
            ann_list = record.INFO['ANN']  # Extract list of annotations
            for ann in ann_list:  # Iterate over each annotation
                fields = ann.split('|')  # SnpEff annotation format: split by '|'
                if fields[2] == impact_type:  # Compare impact type field
                    count += 1  # Increment counter if impact type matches
                    break  # Only count variant once even if multiple transcripts
    return count  # Return total count of variants with given impact type


def compute_cds_length_in_window(gff_file, scaffold, window_start, window_end):
    """
    Compute total CDS length in a specific window from a GFF file.

    Parameters:
        gff_file (str): path to GFF annotation
        scaffold (str): scaffold/chromosome name
        start (int): window start position (1-based)
        end (int): window end position (inclusive)

    Returns:
        int: total length of CDS features within the window
    """
    db = gffutils.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    cds_length = 0  # Initialize counter for CDS length
    for feature in db.region(region=(scaffold, window_start, window_end), featuretype='CDS'):  # Iterate over CDS features in window
        # Determine the portion of CDS overlapping the window
        left_edge_overlap = max(feature.start, window_start)  # If feature starts before window, crop left side; else use feature.start
        right_edge_overlap = min(feature.end, window_end)     # If feature ends after window, crop right side; else use feature.end
        if left_edge_overlap <= right_edge_overlap:  # Only add if there is actual overlap
            cds_length += right_edge_overlap - left_edge_overlap + 1  # Add only overlapping portion of CDS
    return cds_length  # Return total CDS length within window

def count_corrected_variants(vcf_file, gff_file, scaffold, start, end, impact_type):
    """
    Compute corrected number of variants per window, normalizing by CDS length.

    Parameters:
        vcf_file (str): path to VCF file
        gff_file (str): path to GFF annotation
        scaffold (str): scaffold/chromosome name
        start (int): window start position
        end (int): window end position
        impact_type (str): 'HIGH', 'MODERATE', 'LOW', 'MODIFIER'

    Returns:
        float: corrected variant count (variants per base of CDS)
    """
    raw_count = count_mutation_type(vcf_file, impact_type)
    cds_length = compute_cds_length_in_window(gff_file, scaffold, start, end)
    corrected_count = raw_count / cds_length if cds_length > 0 else 0.0
    return raw_count, corrected_count

def write_window_result(outdir, scaffold, start, end, raw_count, corrected_count, impact_type):
    """
    Append window result to a single tab-delimited file without erasing it.
    """
    results_file = os.path.join(outdir, "all_windows_results.tsv")  # Master file path
    header = ["scaffold","window_start","window_end","raw_count","corrected_count","impact_type"]

    os.makedirs(outdir, exist_ok=True)  # Ensure output directory exists

    file_exists = os.path.isfile(results_file)  # Check if the file exists

    with open(results_file, "a", newline="") as f:  # Open in append mode
        if not file_exists:
            f.write("\t".join(header) + "\n")  # Write header if missing
        line = [str(scaffold), str(start), str(end),
                str(raw_count), str(corrected_count), str(impact_type)]
        f.write("\t".join(line) + "\n")  # Append result line

if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # Argument parser
    parser.add_argument('--vcf', required=True, help='VCF file path')  # VCF
    parser.add_argument('--gff', required=True, help='GFF file path')  # GFF
    parser.add_argument('--scaffold', required=True, help='Scaffold name')  # Scaffold
    parser.add_argument('--start', type=int, required=True, help='Window start')  # Window start
    parser.add_argument('--end', type=int, required=True, help='Window end')  # Window end
    parser.add_argument('--impact', required=True, help='Mutation impact to count')  # Impact type
    parser.add_argument('--outdir', required=True, help='Output directory')  # Output dir
    args = parser.parse_args()  # Parse args

    raw_count, corrected_count = count_corrected_variants(
        args.vcf, args.gff, args.scaffold, args.start, args.end, args.impact
    )  # Compute counts

    os.makedirs(args.outdir, exist_ok=True)  # Ensure output dir exists

    write_window_result(
        outdir=args.outdir,
        scaffold=args.scaffold,
        start=args.start,
        end=args.end,
        raw_count=raw_count,
        corrected_count=corrected_count,
        impact_type=args.impact
    )  # Write to master tab-delimited file