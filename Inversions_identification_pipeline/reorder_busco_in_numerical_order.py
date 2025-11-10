import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep="\t", header=None)

df["start"] = df[0].str.extract(r":(\d+)-").astype(int)

df = df.sort_values(by="start")

df.to_csv(output_file, sep="\t", header=False, index=False)
