rule prepare_bwa_index:
    input:
        rules.typing_and_copy_reference.output.ref,
    output:
        multiext("prepare/{sample}.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        ".log/prepare/{sample}.bwa_index.log",
    benchmark:
        ".log/prepare/{sample}.bwa_index.bm"
    conda:
        config["conda"]["bwa"]
    shell:
        "bwa index {input} 2> {log}"


rule prepare_samtools_faidx:
    input:
        rules.typing_and_copy_reference.output.ref,
    output:
        "prepare/{sample}.fasta.fai",
    log:
        ".log/prepare/{sample}.samtools_faidx.log",
    benchmark:
        ".log/prepare/{sample}.samtools_faidx.bm"
    conda:
        config["conda"]["samtools"]
    shell:
        "samtools faidx {input} 2> {log}"


rule prepare_tantan:
    input:
        rules.typing_and_copy_reference.output.ref,
    output:
        "prepare/{sample}.tantan.fasta",
    log:
        ".log/prepare/{sample}.prepare_tantan.log",
    benchmark:
        ".log/prepare/{sample}.prepare_tantan.bm"
    conda:
        config["conda"]["tantan"]
    shell:
        "tantan {input} > {output} 2> {log}"


rule prepare_extract_lcr_bed:
    input:
        rules.prepare_tantan.output,
    output:
        "prepare/{sample}.lcr_regions.bed",
    log:
        ".log/prepare/{sample}.prepare_extract_lcr_bed.log",
    benchmark:
        ".log/prepare/{sample}.prepare_extract_lcr_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/tantan_lcr_extractor.py"


rule prepare_bwa_mem_primer:
    input:
        primer=f"{workflow.basedir}/assets/primer.fasta",
        ref=rules.typing_and_copy_reference.output.ref,
        idx=rules.prepare_bwa_index.output,
    output:
        temp("prepare/{sample}.primer.sam"),
    log:
        ".log/prepare/{sample}.prepare_bwa_mem_primer.log",
    benchmark:
        ".log/prepare/{sample}.prepare_bwa_mem_primer.bm"
    params:
        "-k 8 -T 8",
    threads: config["threads"]["low"]
    conda:
        config["conda"]["bwa"]
    shell:
        "bwa mem -t {threads} {params} {input.ref} {input.primer} > {output} 2> {log}"


rule prepare_make_primer_bed:
    input:
        sam=rules.prepare_bwa_mem_primer.output[0],
        fai=rules.prepare_samtools_faidx.output[0],
    output:
        mapped=temp("prepare/{sample}.mapped.sam"),
        bed="prepare/{sample}.primer.bed",
    log:
        ".log/prepare/{sample}.prepare_make_primer_bed.log",
    benchmark:
        ".log/prepare/{sample}.prepare_make_primer_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_primer_bed.py"


rule prepare_primer_bedtools_slop:
    input:
        bed=rules.prepare_make_primer_bed.output.bed,
        fai=rules.prepare_samtools_faidx.output,
    output:
        primer="prepare/{sample}.primer.mask.bed",
    log:
        ".log/prepare/{sample}.prepare_primer_bedtools_slop.log",
    benchmark:
        ".log/prepare/{sample}.prepare_primer_bedtools_slop.bm"
    params:
        flank=50,
    conda:
        config["conda"]["bedtools"]
    shell:
        "bedtools slop -b {params.flank} -i {input.bed} -g {input.fai} > {output.primer} 2> {log}"


rule prepare_make_target_bed:
    input:
        rules.prepare_primer_bedtools_slop.output.primer,
    output:
        "prepare/{sample}.target.bed",
    log:
        ".log/prepare/{sample}.prepare_make_target_bed.log",
    benchmark:
        ".log/prepare/{sample}.prepare_make_target_bed.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/make_target_bed.py"
