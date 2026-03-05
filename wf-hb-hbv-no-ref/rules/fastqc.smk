rule fastqc:
    input:
        rules.create_symlinks.output.fq1,
        rules.create_symlinks.output.fq2,
    output:
        directory("qc/fastqc/{sample}"),
    benchmark:
        ".log/qc/fastqc/{sample}.bm"
    log:
        ".log/qc/fastqc/{sample}.log",
    conda:
        config["conda"]["fastqc"]
    threads: config["threads"]["low"]
    shell:
        "mkdir {output} && fastqc {input} -o {output} -t {threads} --extract &> {log}"


rule multiqc:
    input:
        expand("qc/fastqc/{sample}", sample=samples),
    output:
        directory("qc/multiqc"),
    benchmark:
        ".log/qc/multiqc/multiqc.bm"
    log:
        ".log/qc/multiqc/multiqc.log",
    conda:
        config["conda"]["multiqc"]
    shell:
        "multiqc {input} --outdir {output} 2> {log}"


rule fastp:
    input:
        rules.create_symlinks.output.fq1,
        rules.create_symlinks.output.fq2,
    output:
        fq1="qc/fastp/{sample}.1.fastq.gz",
        fq2="qc/fastp/{sample}.2.fastq.gz",
        html="qc/fastp/{sample}.html",
        json="qc/fastp/{sample}.json",
    log:
        ".log/qc/fastp/{sample}.fastp.log",
    benchmark:
        ".log/qc/fastp/{sample}.fastp.bm"
    conda:
        config["conda"]["fastp"]
    threads: config["threads"]["low"]
    params:
        # * [20251013] 浩博 (ChatGPT) 建议的参数
        extra=(
            "--cut_front --cut_tail --cut_mean_quality 20 --qualified_quality_phred 20 "
            "--length_required 75 --detect_adapter_for_pe --trim_poly_g"
        ),
    shell:
        """
        fastp -w {threads} {params.extra} \
            -i {input[0]} -I {input[1]} \
            -o {output.fq1} -O {output.fq2} \
            -h {output.html} -j {output.json} 2> {log}
        """


# Description:  对 cleaned fastq 进行 fastp 质控, 查看 adapter 残留比例.
#               若 cleaned fastq 为空, 则创建空的 trimmed fastq 文件, html 文件和 json 文件.
rule cleaned_fastp:
    input:
        rules.fastp.output.fq1,
        rules.fastp.output.fq2,
    output:
        fq1=temp("qc/fastp/{sample}.trimmed.1.fastq.gz"),
        fq2=temp("qc/fastp/{sample}.trimmed.2.fastq.gz"),
        html=temp("qc/fastp/{sample}.trimmed.html"),
        json="qc/fastp/{sample}.trimmed.json",
    log:
        ".log/qc/fastp/{sample}.cleaned_fastp.log",
    benchmark:
        ".log/qc/fastp/{sample}.cleaned_fastp.bm"
    conda:
        config["conda"]["fastp"]
    threads: config["threads"]["low"]
    params:
        extra="",
    run:
        import gzip
        import os

        if (os.stat(input[0]).st_size <= 100) or (os.stat(input[1]).st_size <= 100):
            # 创建空的 trimmed fastq 文件
            for out_file in [output.fq1, output.fq2]:
                with gzip.open(out_file, "wb") as f:
                    f.write(b"")
                    # 创建空的 html 文件
            with open(output.html, "w") as f:
                f.write("<html><body>No reads after first filtering</body></html>")
                # 创建空的 json 文件
            with open(output.json, "w") as f:
                f.write(
                    '{"adapter_cutting": {"adapter_trimmed_reads": 0}, "summary": {"before_filtering": {"total_reads": 0}}}'
                )
        else:
            # 正常运行 fastp
            shell(
                """
            fastp -w {threads} {params.extra} \
                -i {input[0]} -I {input[1]} \
                -o {output.fq1} -O {output.fq2} \
                -h {output.html} -j {output.json} 2> {log}
            """
            )


rule fq_stats_summary:
    input:
        fastp_json=expand("qc/fastp/{sample}.json", sample=samples),
        cleaned_json=expand("qc/fastp/{sample}.trimmed.json", sample=samples),
    output:
        "qc/fastp/fq_summary.tsv",
    benchmark:
        ".log/qc/fastp/fq_stats_summary.bm"
    log:
        ".log/qc/fastp/fq_stats_summary.log",
    conda:
        config["conda"]["python"]
    script:
        "../scripts/fq_stats_summary.py"
