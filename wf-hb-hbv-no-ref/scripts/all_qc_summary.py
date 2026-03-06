import pandas as pd
import sys
import numpy as np

sys.stderr = open(snakemake.log[0], "w")

dfqc = pd.read_csv(str(snakemake.input.fq), sep="\t")
dfbam = pd.read_csv(str(snakemake.input.bam), sep="\t")
dftyping = pd.read_csv(str(snakemake.input.typing), sep="\t")
dfmerged = pd.merge(dftyping, dfqc, on='Sample', how='left')
dfmerged = pd.merge(dfmerged, dfbam, on='Sample', how='left')

# [20260112 HWY SOP QC 参数] 1.≥0.35Gb, 2. Q30≥85%
# 使用numpy的where函数生成QC结果，大型数据集下性能最优
dfmerged.insert(
    1,
    'QCJudgment',
    np.where((dfmerged['RawBases'] >= 0.35*1e9) &
             (dfmerged['CleanQ30'] >= 0.85), 'Pass', 'Fail')
)

# 输出结果
dfmerged.to_csv(str(snakemake.output.tsv), index=False, sep='\t')
dfmerged.to_excel(str(snakemake.output.excel), index=False)
