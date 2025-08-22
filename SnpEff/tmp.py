# Step 4-6: Python script to count variants and normalize by CDS length
# pseudocode / blueprint

import vcf  # pyVCF or equivalent
import pandas as pd

# Load CDS coordinates per scaffold (could come from GFF)
cds_df = pd.read_csv("cds_coordinates.txt")  # scaffold, start, end

windows = pd.read_csv("windows.txt", sep="\t", names=["scaffold","start","end"])
results = []

for index, row in windows.iterrows():
    scaffold, w_start, w_end = row.scaffold, row.start, row.end
    
    # Determine CDS length overlapping this window
    overlapping_cds = cds_df[(cds_df.scaffold == scaffold) &
                              (cds_df.end >= w_start) &
                              (cds_df.start <= w_end)]
    
    cds_length = 0
    for _, cds in overlapping_cds.iterrows():
        cds_start = max(cds.start, w_start)
        cds_end = min(cds.end, w_end)
        cds_length += cds_end - cds_start + 1
    
    # Parse VCF and count variants by IMPACT
    vcf_reader = vcf.Reader(open(f"sub_window_{scaffold}_{w_start}_{w_end}.vcf"))
    counts = {"HIGH":0, "MODERATE":0, "LOW":0}
    
    for record in vcf_reader:
        impact = record.INFO['ANN'][0].split('|')[2]  # SnpEff IMPACT field
        if impact in counts:
            counts[impact] += 1
    
    # Normalize by CDS length
    for key in counts:
        counts[key] = counts[key] / cds_length
    
    results.append({"scaffold":scaffold, "start":w_start, "end":w_end, **counts})

# Step 7: Aggregate results
results_df = pd.DataFrame(results)
results_df.to_csv("variant_density_per_window.csv", index=False)
