import sys
import pandas as pd
from pathlib import Path

sys.stderr = open(snakemake.log[0], "w")

# IO
blastn_res = snakemake.input[0]
hbv_type = snakemake.output['hbv_type']
ref = snakemake.output['ref']

# 登录号-分型 映射表
acc_type_dict = snakemake.params.acc_type_dict
# 序列目录路径
seq_dir = snakemake.params.seq_dir

# 判断样本分型
df = pd.read_table(blastn_res, sep="\t", header=None,
                   names=["qseqid", "sseqid", "pident", "length", "bitscore"])
# 使用 pident * length 作为排序的指标，选择最高的作为样本分型
df['pident_length_product'] = df['pident'] * df['length']
df.sort_values(by='pident_length_product', ascending=False, inplace=True)
top1_acc = df.iloc[0]['sseqid']
top1_type = acc_type_dict[top1_acc]

# 输出样本分型
with open(hbv_type, "w") as f:
    f.write(f"{top1_type}\n")

# 复制参考序列
ref_origin_path = Path(seq_dir) / f"{top1_acc}.fasta"
with open(ref_origin_path, "r") as f, open(ref, "w") as out:
    out.write(f.read())
