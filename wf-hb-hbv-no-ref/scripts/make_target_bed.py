
import sys

sys.stderr = open(snakemake.log[0], "w")

in_bed = snakemake.input[0]
out_bed = snakemake.output[0]

# 初始化染色体、左边界、右边界
chrom, left_poss, right_poss = '', [], []

with open(in_bed, "r") as f:
    for line in f:
        line = line.strip()
        fields = line.split("\t")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        orient = fields[5]
        if orient == "+":
            left_poss.append(end)
        elif orient == "-":
            right_poss.append(start)

with open(out_bed, "w") as f:
    f.write(f"{chrom}\t{max(left_poss)}\t{min(right_poss)}\n")
