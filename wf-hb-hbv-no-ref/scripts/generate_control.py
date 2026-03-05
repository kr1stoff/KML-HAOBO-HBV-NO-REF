# 合并所有样本变异结果，然后使用 Chrom, Pos, Ref, Alt 作为 key 合并 AltFreq。然后计算每个变异在批次内的出现频率

import pandas as pd
from functools import reduce
from pathlib import Path

import sys
sys.stderr = open(snakemake.log[0], "w")


input_tabs = snakemake.input
output_tab = snakemake.output[0]


dfs = []
for tab in input_tabs:
    df = pd.read_csv(tab, sep='\t', usecols=[
                     "Chrom", "Pos", "Ref", "Alt", "AltFreq"])
    df.rename(columns={"AltFreq": Path(tab).stem.split('.')[0]}, inplace=True)
    dfs.append(df)
df = reduce(lambda x, y: pd.merge(
    x, y, on=["Chrom", "Pos", "Ref", "Alt"], how="outer"), dfs)
df.set_index(["Chrom", "Pos", "Ref", "Alt"], inplace=True)

# 频率非零的样本数 / 总样本数
batch_ser = round((1 - df.isna().sum(axis=1) / df.shape[1]), 4)
batch_ser.name = 'IntraBatch'
batch_ser.reset_index().to_csv(output_tab, index=False)
