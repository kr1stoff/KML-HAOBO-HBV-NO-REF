rule map_bwa_mem:
    input:
        fq1=rules.fastp.output.fq1,
        fq2=rules.fastp.output.fq2,
        ref=rules.typing_and_copy_reference.output.ref,
        idx=rules.prepare_bwa_index.output,
    output:
        temp("align/{sample}.sam"),
    benchmark:
        ".log/align/{sample}.map_bwa_mem.bm"
    log:
        ".log/align/{sample}.map_bwa_mem.log",
    conda:
        config["conda"]["bwa"]
    params:
        extra=r"-M -Y -R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
    threads: config["threads"]["medium"]
    shell:
        "bwa mem -t {threads} {params.extra} {input.ref} {input.fq1} {input.fq2} -o {output} 2> {log}"


use rule typing_samtools_sort_and_index as map_samtools_sort_and_index with:
    input:
        rules.map_bwa_mem.output,
    output:
        bam="align/{sample}.bam",
        bai="align/{sample}.bam.bai",
    benchmark:
        ".log/align/{sample}.map_samtools_sort_and_index.bm"
    log:
        ".log/align/{sample}.map_samtools_sort_and_index.log",


rule map_samtools_stats_target:
    input:
        bam=rules.map_samtools_sort_and_index.output.bam,
        bed=rules.prepare_make_target_bed.output,
    output:
        "align/{sample}.bam.target.stat",
    benchmark:
        ".log/align/{sample}.map_samtools_stats_target.bm"
    log:
        ".log/align/{sample}.map_samtools_stats_target.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} --target-regions {input.bed} {input.bam} > {output} 2> {log}"


rule map_samtools_stats_all_position:
    input:
        rules.map_samtools_sort_and_index.output.bam,
    output:
        "align/{sample}.bam.stat",
    benchmark:
        ".log/align/{sample}.map_samtools_stats_all_position.bm"
    log:
        ".log/align/{sample}.map_samtools_stats_all_position.log",
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools stats --threads {threads} {input} > {output} 2> {log}"


rule map_samtools_depth_target:
    input:
        bam=rules.map_samtools_sort_and_index.output.bam,
        bed=rules.prepare_make_target_bed.output,
    output:
        "align/{sample}.bam.target.depth",
    benchmark:
        ".log/align/{sample}.map_samtools_depth_target.bm"
    log:
        ".log/align/{sample}.map_samtools_depth_target.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools depth -a -b {input.bed} {input.bam} > {output} 2> {log}"


rule map_bam_stats:
    input:
        rules.map_samtools_stats_all_position.output,
        rules.map_samtools_stats_target.output,
        rules.map_samtools_depth_target.output,
    output:
        "align/{sample}.stats.csv",
    benchmark:
        ".log/align/{sample}.map_bam_stats.bm"
    log:
        ".log/align/{sample}.map_bam_stats.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats.py"


rule bam_stats_summary:
    input:
        expand("align/{sample}.stats.csv", sample=samples),
    output:
        "align/bam_summary.tsv",
    benchmark:
        ".log/align/bam_stats_summary.bm"
    log:
        ".log/align/bam_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/bam_stats_summary.py"


rule map_samtools_bedcov_target:
    input:
        bam=rules.map_samtools_sort_and_index.output.bam,
        bed=rules.prepare_make_target_bed.output,
        index=rules.map_samtools_sort_and_index.output.bai,
    output:
        "align/{sample}.bam.target.bedcov",
    benchmark:
        ".log/align/{sample}.map_samtools_bedcov_target.bm"
    log:
        ".log/align/{sample}.map_samtools_bedcov_target.log",
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools bedcov -c {input.bed} {input.bam} > {output} 2> {log}"
