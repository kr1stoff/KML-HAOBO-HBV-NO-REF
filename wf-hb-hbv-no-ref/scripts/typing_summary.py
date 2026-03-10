import sys
import pandas as pd
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

outitems = []

for type_file in snakemake.input:
    sample = Path(type_file).stem.split('.')[0]
    hbv_type = open(type_file, "r").read().strip()
    outitems.append([sample, hbv_type])

df = pd.DataFrame(outitems, columns=["Sample", "HBV_Genoype"])
df.to_csv(snakemake.output[0], index=False, sep="\t")
