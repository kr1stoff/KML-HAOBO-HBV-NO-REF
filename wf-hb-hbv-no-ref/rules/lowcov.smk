use rule typing_bedtools_lowcov_mask as lowcov_mask with:
    input:
        bam=rules.trim_samtools_sort_and_index2.output.bam,
        ref=rules.typing_and_copy_reference.output.ref,
    output:
        genomecov="lowcov/{sample}.genomecov.bed",
        lowcov="lowcov/{sample}.lowcov.bed",
        mask="lowcov/{sample}.masked.fa",
    log:
        ".log/lowcov/{sample}.lowcov_mask.log",
    benchmark:
        ".log/lowcov/{sample}.lowcov_mask.bm"
    params:
        depth_threshold=1000,
        length_threshold=20,


use rule prepare_samtools_faidx as lowcov_samtools_faidx with:
    input:
        rules.lowcov_mask.output.mask,
    output:
        "lowcov/{sample}.masked.fa.fai",
    log:
        ".log/lowcov/{sample}.lowcov_samtools_faidx.log",
    benchmark:
        ".log/lowcov/{sample}.lowcov_samtools_faidx.bm"
