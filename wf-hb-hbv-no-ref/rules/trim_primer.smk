rule trim_primer_ivar:
    input:
        bam=rules.map_bwa_mem.output,
        bed=rules.prepare_make_primer_bed.output.bed,
    output:
        temp("trim-primer/{sample}.trimmed.unsorted.bam"),
    log:
        ".log/trim_primer/{sample}.trim_primer_ivar.log",
    benchmark:
        ".log/trim_primer/{sample}.trim_primer_ivar.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["ivar"]
    shell:
        """
        prefix=$(dirname {output})/$(basename {output} .bam)
        ivar trim -i {input.bam} -b {input.bed} -p $prefix -e 2> {log}
        """


rule trim_primer_ivar_sort_and_index:
    input:
        rules.trim_primer_ivar.output,
    output:
        bam="trim-primer/{sample}.trimmed.bam",
        bai="trim-primer/{sample}.trimmed.bam.bai",
    log:
        ".log/trim_primer/{sample}.trim_primer_ivar_sort_and_index.log",
    benchmark:
        ".log/trim_primer/{sample}.trim_primer_ivar_sort_and_index.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["samtools"]
    shell:
        """
        # ! ivar trim 如果输入没有条目会删掉信息头输出空文件
        # 检查unsorted文件是否有reads(除了信息头)
        if [[ $(samtools view -c {input} 2>/dev/null || echo 0) -gt 0 ]]; then
            # 有reads时正常处理
            samtools sort -@ {threads} -o {output.bam} {input} 2>> {log}
        else
            # 没有reads时创建空的已排序BAM文件
            echo "警告: {input} 没有reads，创建空输出" >> {log}
            samtools view -hbS {input} > {output.bam} 2>> {log}
        fi
        samtools index {output.bam} 2>> {log}
        """


rule trim_primer_bedtools_intersect:
    input:
        bam=rules.trim_primer_ivar_sort_and_index.output.bam,
        bed=rules.prepare_primer_bedtools_slop.output.primer,
    output:
        temp("trim-primer/{sample}.filtered.unsorted.bam"),
    log:
        ".log/trim_primer/{sample}.trim_primer_bedtools_intersect.log",
    benchmark:
        ".log/trim_primer/{sample}.trim_primer_bedtools_intersect.bm"
    threads: config["threads"]["low"]
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools intersect -v -a {input.bam} -b {input.bed} > {output} 2> {log}"


use rule typing_samtools_sort_and_index as trim_primer_sort_and_index with:
    input:
        rules.trim_primer_bedtools_intersect.output,
    output:
        bam="trim-primer/{sample}.filtered.bam",
        bai="trim-primer/{sample}.filtered.bam.bai",
    log:
        ".log/trim_primer/{sample}.trim_primer_sort_and_index.log",
    benchmark:
        ".log/trim_primer/{sample}.trim_primer_sort_and_index.bm"
