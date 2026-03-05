# KML-HAOBO-HBV-NO-REF

HAOBO HBV 靶向基因测序，无预设参考分析流程。先做A-H分型，再分析变异

## 命令行

### 运行

在 `wf-hb-hbv-no-ref` 目录下运行以下命令

```bash
snakemake --cores 32 --use-conda --rerun-incomplete --scheduler greedy \
  --config samples_tsv=$PWD/tests/example.tsv \
  --resources freebayes_jobs=10 \
  --directory /data/mengxf/Project/KML260302-HAOBOHBV-NO-REF/results/260303
```

### 自定义参数

```text
--config sample_tsv:          制表符分割的样本信息表, 包含样本名称, FQ1, FQ2路径
--resources freebayes_jobs:   freebayes 并行线程数, 默认 10
```

## 开发

### typing

分型部分分析流程.

1. 先抽取 R1+R2 共 100K reads
2. 比对到默认 B 型参考基因组
3. 屏蔽低深度区域, 100X 覆盖度以下的区域
4. `bcftools` 分析变异
5. 生成一致性序列
6. `blastn` 比对到 A-H 分型参考基因组
7. 按照 identity,length 最高的比对结果作为最终分型
