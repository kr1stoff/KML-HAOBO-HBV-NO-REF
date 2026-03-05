rule seqtk_sample_50k:
    input:
        rules.fastp.output.fq1,
    output:
        temp("{sample}.50k.1.fq"),
    log:
        ".log/typing/{sample}.seqtk_sample_50k.log",
    benchmark:
        ".log/typing/{sample}.seqtk_sample_50k.bm"
    conda:
        config["conda"]["seqtk"]
    shell:
        "seqtk sample {input} 50000 > {output} 2> {log}"


use rule seqtk_sample_50k as seqtk_sample_50k_r2 with:
    input:
        rules.fastp.output.fq2,
    output:
        temp("{sample}.50k.2.fq"),
    log:
        ".log/typing/{sample}.seqtk_sample_50k_r2.log",
    benchmark:
        ".log/typing/{sample}.seqtk_sample_50k_r2.bm"


rule typing_bwa_mem:
    input:
        fq1=rules.seqtk_sample_50k.output,
        fq2=rules.seqtk_sample_50k_r2.output,
        ref=config["database"]["ref"],
    output:
        temp("typing/{sample}.sam"),
    benchmark:
        ".log/typing/{sample}.typing_bwa_mem.bm"
    log:
        ".log/typing/{sample}.typing_bwa_mem.log",
    conda:
        config["conda"]["bwa"]
    threads: config["threads"]["medium"]
    shell:
        "bwa mem -t {threads} {input.ref} {input.fq1} {input.fq2} -o {output} 2> {log}"


rule typing_samtools_sort_and_index:
    input:
        rules.typing_bwa_mem.output,
    output:
        bam=temp("typing/{sample}.sorted.bam"),
        bai=temp("typing/{sample}.sorted.bam.bai"),
    log:
        ".log/typing/{sample}.typing_samtools_sort_and_index.log",
    benchmark:
        ".log/typing/{sample}.typing_samtools_sort_and_index.bm"
    conda:
        config["conda"]["samtools"]
    shell:
        """
        samtools sort -o {output.bam} {input} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule typing_bedtools_lowcov_mask:
    input:
        bam=rules.typing_bwa_mem.output,
        ref=config["database"]["ref"],
    output:
        genomecov="typing/{sample}.genomecov.bed",
        lowcov="typing/{sample}.lowcov.bed",
        mask="typing/{sample}.masked.fa",
    log:
        ".log/typing/{sample}.typing_bedtools_lowcov_mask.log",
    benchmark:
        ".log/typing/{sample}.typing_bedtools_lowcov_mask.bm"
    conda:
        config["conda"]["bedtools"]
    params:
        depth_threshold=100,
        length_threshold=20,
    shell:
        """
        bedtools genomecov -bga -ibam {input.bam} > {output.genomecov} 2> {log}
        awk '$4<{params.depth_threshold}' {output.genomecov} | bedtools merge | awk '$3-$2>{params.length_threshold}' > {output.lowcov} 2>> {log}
        bedtools maskfasta -fi {input.ref} -bed {output.lowcov} -fo {output.mask} 2>> {log}
        """


rule typing_bcftools_concensus:
    input:
        bam=rules.typing_samtools_sort_and_index.output.bam,
        bai=rules.typing_samtools_sort_and_index.output.bai,
        ref=config["database"]["ref"],
        mask=rules.typing_bedtools_lowcov_mask.output.mask,
        lowcov=rules.typing_bedtools_lowcov_mask.output.lowcov,
    output:
        consensus="typing/{sample}.consensus.fasta",
        vcf="typing/{sample}.vcf.gz",
    log:
        ".log/typing/{sample}.typing_bcftools_concensus.log",
    benchmark:
        ".log/typing/{sample}.typing_bcftools_concensus.bm"
    conda:
        config["conda"]["bcftools"]
    shell:
        """
        bcftools mpileup -f {input.ref} -a FORMAT/DP -Ou {input.bam} | bcftools call -mv -Oz -o {output.vcf} 2> {log}
        bcftools index {output.vcf} 2>> {log}
        bcftools consensus -f {input.mask} -m {input.lowcov} {output.vcf} > {output.consensus} 2>> {log}
        """


rule typing_blastn:
    input:
        consensus=rules.typing_bcftools_concensus.output.consensus,
        blastn=config["database"]["blastn"],
    output:
        "typing/{sample}.blastn.txt",
    log:
        ".log/typing/{sample}.typing_blastn.log",
    benchmark:
        ".log/typing/{sample}.typing_blastn.bm"
    conda:
        config["conda"]["blastn"]
    shell:
        """
        blastn -max_target_seqs 8 -num_threads 16 \
            -query {input.consensus} \
            -db {input.blastn} \
            -outfmt "6 qseqid sseqid pident length bitscore" \
            -out {output} 2> {log}
        """


rule typing_and_copy_reference:
    input:
        rules.typing_blastn.output[0],
    output:
        hbv_type="typing/{sample}.typing.txt",
        ref="prepare/{sample}.fasta",
    log:
        ".log/typing/{sample}.typing_and_copy_reference.log",
    benchmark:
        ".log/typing/{sample}.typing_and_copy_reference.bm"
    conda:
        config["conda"]["python"]
    params:
        acc_type_dict=config["custom"]["accession_type"],
        seq_dir=config["custom"]["sequences_dir"],
    script:
        "../scripts/type_and_copy_ref.py"


rule typing_summary:
    input:
        expand("typing/{sample}.typing.txt", sample=samples),
    output:
        "typing/typing_summary.tsv",
    log:
        ".log/typing/typing_summary.log",
    benchmark:
        ".log/typing/typing_summary.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/typing_summary.py"
