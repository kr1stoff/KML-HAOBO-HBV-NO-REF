import pandas as pd
from pathlib import Path
import sys
sys.stderr = open(snakemake.log[0], "w")

input_tabs = snakemake.input
output_excel = snakemake.output[0]


with pd.ExcelWriter(output_excel) as writer:
    for it in input_tabs:
        sample_name = Path(it).stem.split('.')[0]
        df = pd.read_csv(it, sep='\t')
        df.to_excel(writer, sheet_name=sample_name, index=False)
