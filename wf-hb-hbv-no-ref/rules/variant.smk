rule freebayes:
    input:
        alns=rules.trim_samtools_sort_and_index2.output.bam,
        ref=rules.lowcov_mask.output.mask,
        fai=rules.lowcov_samtools_faidx.output,
        targets=rules.prepare_make_target_bed.output,
    output:
        vcf="variant/{sample}.freebayes.vcf",
    log:
        ".log/variant/{sample}.freebayes.log",
    benchmark:
        ".log/variant/{sample}.freebayes.bm"
    # resources 中 freebayes_jobs 控制 freebayes 的并行线程数，减小内存压力
    resources:
        freebayes_jobs=1
    # [20251014] 浩博建议的参数
    params:
        (
            "--pooled-continuous --min-repeat-size 10 --read-indel-limit 15 --use-best-n-alleles 4 "
            "--theta 0.005 --haplotype-length 0 --min-alternate-fraction 0.005 --min-base-quality 30 "
            "--min-coverage 1000 --min-alternate-count 10 --min-mapping-quality 30 --max-complex-gap 1 --trim-complex-tail"
        ),
    conda:
        config["conda"]["freebayes"]
    shell:
        "freebayes {params} --targets {input.targets} --fasta-reference {input.ref} {input.alns} > {output.vcf} 2> {log}"


# 多等位基因变异拆成多个单等位基因变异
rule variant_bcftools_norm:
    input:
        rules.freebayes.output,
    output:
        "variant/{sample}.norm.vcf",
    log:
        ".log/variant/{sample}.variant_bcftools_norm.log",
    benchmark:
        ".log/variant/{sample}.variant_bcftools_norm.bm"
    conda:
        config["conda"]["bcftools"]
    shell:
        "bcftools norm -m -both {input} > {output} 2> {log}"


# MNP 拆成 SNP
rule variant_vt_decompose_blocksub:
    input:
        rules.variant_bcftools_norm.output,
    output:
        "variant/{sample}.decomposed.vcf",
    log:
        ".log/variant/{sample}.variant_vt_decompose_blocksub.log",
    benchmark:
        ".log/variant/{sample}.variant_vt_decompose_blocksub.bm"
    conda:
        config["conda"]["vt"]
    shell:
        "vt decompose_blocksub {input} > {output} 2> {log}"


rule variant_vcf2tab:
    input:
        rules.variant_vt_decompose_blocksub.output,
    output:
        "variant/{sample}.vcf2tab.tsv",
    log:
        ".log/variant/{sample}.variant_vcf2tab.log",
    benchmark:
        ".log/variant/{sample}.variant_vcf2tab.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf2tab.py"


# Description: 变异过滤(错误模型)第一步
#               1. 低质量: 深度 < 1000X, 等位基因深度 < 10X
#               2. 低于检测限: 最小等位基因频率 < 1%
#               3. 链偏倚: 正链/负链支持的 reads 数量 < 2, 正链/负链支持的 reads < 该变异总 reads 的 10%
#               4. 位置偏倚: 变异位置在 read 左侧或右侧支持的测序深度 < 2, 变异位置在 read 左侧或右侧支持的测序深度 < 该变异总深度的 10%
#               5. 预先注释好基因组低复杂度 BED 区域, 注释 LowComplexity 标签
# Date: 20251114
rule variant_vcf_filter:
    input:
        rules.variant_vcf2tab.output,
        rules.prepare_extract_lcr_bed.output,
    output:
        "variant/{sample}.filter.tsv",
    log:
        ".log/variant/{sample}.variant_vcf_filter.log",
    benchmark:
        ".log/variant/{sample}.variant_vcf_filter.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf_filter.py"


# Description: 变异检验
#               1. 二项分布, scipy.stats.binom_test, 使用 测序深度, 最小等位基因频率, 测序错误率
#               2. 泊松分布, scipy.stats.poisson_test, 使用 测序深度, 最小等位基因频率, 测序错误率
# Date: 20251015
rule variant_test:
    input:
        rules.variant_vcf_filter.output,
    output:
        "variant/{sample}.test.tsv",
    log:
        ".log/variant/{sample}.variant_test.log",
    benchmark:
        ".log/variant/{sample}.variant_test.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/vcf_test.py"


# Description: 生成质控信息, 目前只有批内质控
# Date: 20251202
rule variant_generate_control:
    input:
        expand("variant/{sample}.test.tsv", sample=samples),
    output:
        "variant/control.csv",
    log:
        ".log/variant/variant_generate_control.log",
    benchmark:
        ".log/variant/variant_generate_control.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/generate_control.py"


rule variant_control:
    input:
        rules.variant_test.output,
        rules.variant_generate_control.output,
    output:
        "variant/{sample}.control.tsv",
    log:
        ".log/variant/{sample}.variant_control.log",
    benchmark:
        ".log/variant/{sample}.variant_control.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/variant_control.py"
