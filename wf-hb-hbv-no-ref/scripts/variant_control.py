import pandas as pd
import sys
sys.stderr = open(snakemake.log[0], "w")


input_tab = snakemake.input[0]
batch_control_tab = snakemake.input[1]
output_tab = snakemake.output[0]

df = pd.read_csv(input_tab, sep="\t")
batch_df = pd.read_csv(batch_control_tab)
outdf = pd.merge(df, batch_df, on=["Chrom", "Pos", "Ref", "Alt"], how="left")
outdf.to_csv(output_tab, sep="\t", index=False)
