import pandas as pd
import sys

sys.stderr = open(snakemake.log[0], "w")

in_tab = snakemake.input[0]
lcr_bed = snakemake.input[1]
out_tab = snakemake.output[0]

df = pd.read_table(in_tab)
# 获取低复杂度区域 BED 区域
low_cpx_regions = {}
with open(lcr_bed, "r") as f:
    for line in f:
        fields = line.strip().split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        low_cpx_regions[chrom] = low_cpx_regions.get(chrom, [])
        low_cpx_regions[chrom].append((start, end))
with open(out_tab, "w") as f:
    f.write('\t'.join(df.columns.values.tolist() + ['Filter']) + '\n')
    for item in df.itertuples():
        filters = []
        # ! 过滤标签
        if (item.TotalDepth < 1000) or (item.AltDepth < 10):
            filters.append("LowQuality")
        if item.AltFreq < 0.01:
            filters.append("BelowLOD")
        if (item.SAF < 2) or (item.SAR < 2) or (item.ForwardStrandRate < 0.1) or (item.ReverseStrandRate < 0.1):
            filters.append("StrandBias")
        if (item.RPL < 2) or (item.RPR < 2) or (item.PlacedLeftRate < 0.1) or (item.PlacedRightRate < 0.1):
            filters.append("PosBias")
        if item.Chrom in low_cpx_regions.keys():
            for region in low_cpx_regions[item.Chrom]:
                if region[0] <= (item.Pos - 1) <= region[1]:
                    filters.append("LowComplexity")
        filter_tag = ";".join(filters) if filters else "PASS"
        f.write('\t'.join([str(item)
                for item in item[1:]] + [filter_tag]) + '\n')
