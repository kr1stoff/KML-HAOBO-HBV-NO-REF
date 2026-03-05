# 变异位点检验. 二项分布, 泊松分布
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.stats.multitest as multi
import sys

sys.stderr = open(snakemake.log[0], "w")

var_file = snakemake.input[0]
out_file = snakemake.output[0]


def multiple_testing_correction(pvalues, alpha=0.05, method='fdr_bh'):
    """
    对变异检测结果进行多重检验校正
    :param pvalues: 所有变异位点的p值列表
    :param alpha: 显著性水平(默认0.05)
    :param method: 多重检验校正方法(默认FDR-BH)
    :return corrected_pvals: 校正后的p值数组
    """
    try:
        # FDR校正
        _, corrected_pvals, _, _ = multi.multipletests(
            pvalues, alpha=alpha, method=method)
    except ZeroDivisionError:
        sys.stderr.write(
            "Warning: ZeroDivisionError in multiple_testing_correction. Setting all p-values to 1.0.\n")
        corrected_pvals = np.ones_like(pvalues)
    return corrected_pvals


# main
df = pd.read_csv(var_file, sep='\t')
total = df['TotalDepth'].to_numpy(dtype=int)
alt = df['AltDepth'].to_numpy(dtype=int)
# 错误率, 多重检验校正参数
error_rate = 0.002
mtc_alpha = 0.05
mtc_method = 'fdr_bh'

# 二项分布
with np.errstate(invalid='ignore', divide='ignore'):
    binom_pvals = stats.binom.sf(alt - 1, total, error_rate)
binom_pvals = np.where(total <= 0, 1.0, binom_pvals)

# 泊松分布
expected_errors = error_rate * total
with np.errstate(invalid='ignore', divide='ignore'):
    poisson_pvals = stats.poisson.sf(alt - 1, expected_errors)
poisson_pvals = np.where(total <= 0, 1.0, poisson_pvals)

# 多重检验校正 (FDR-BH)
binom_adj = multiple_testing_correction(binom_pvals)
poisson_adj = multiple_testing_correction(poisson_pvals)

# 写入 DataFrame
df['PvalBinomial'] = binom_pvals
df['PvalBinomialAdj'] = binom_adj
df['PvalPoisson'] = poisson_pvals
df['PvalPoissonAdj'] = poisson_adj

# 输出
df.to_csv(out_file, sep='\t', index=False)
