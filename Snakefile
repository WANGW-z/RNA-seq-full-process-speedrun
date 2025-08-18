import pandas as pd

# 加载配置和样本信息
configfile: "config/config.json"
samples = pd.read_csv(config["samples"], index_col="sample")


# 定义最终输出目标
rule all:
    input:
        "results/multiqc_report.html",
        "results/counts/gene_counts.tsv",
        "results/deseq2/DE_results.csv",
        expand("results/fastqc/trimmed/{sample}_R{read}_fastqc.html", sample=samples.index, read=[1, 2])  # 确保post-trim QC执行

# ----------- 质控1 -----------
rule fastqc_raw:
    input:
        fq1 = lambda wildcards: samples.loc[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.loc[wildcards.sample, "fq2"]
    output:
        html1 = "results/fastqc/raw/{sample}_R1_fastqc.html",
        zip1 = "results/fastqc/raw/{sample}_R1_fastqc.zip",
        html2 = "results/fastqc/raw/{sample}_R2_fastqc.html", 
        zip2 = "results/fastqc/raw/{sample}_R2_fastqc.zip"
    log:
        "logs/fastqc/raw/{sample}.log"
    params:
        prefix = "E250715005_L01_"
    threads: 2
    shell:
        """
        # 创建目录并运行FastQC（所有命令写在一行，避免换行问题）
        mkdir -p results/fastqc/raw && \
        fastqc {input.fq1} {input.fq2} -o results/fastqc/raw/ 2> {log} && \
        mv "results/fastqc/raw/{params.prefix}{wildcards.sample}_1_fastqc.html" {output.html1} && \
        mv "results/fastqc/raw/{params.prefix}{wildcards.sample}_1_fastqc.zip" {output.zip1} && \
        mv "results/fastqc/raw/{params.prefix}{wildcards.sample}_2_fastqc.html" {output.html2} && \
        mv "results/fastqc/raw/{params.prefix}{wildcards.sample}_2_fastqc.zip" {output.zip2}
        """


# ----------- 预处理 -----------
rule trim:
    input:
        fq1 = lambda wildcards: samples.loc[wildcards.sample, "fq1"],
        fq2 = lambda wildcards: samples.loc[wildcards.sample, "fq2"]
    output:
        r1 = "results/trimmed/{sample}_R1_clean.fastq.gz",
        r1_unpaired = "results/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2 = "results/trimmed/{sample}_R2_clean.fastq.gz",
        r2_unpaired = "results/trimmed/{sample}_R2_unpaired.fastq.gz"
    log:
        "logs/trim/{sample}.log"
    threads: 4
    conda:
        "envs/rnaseq.yml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 "
        "{input.fq1} {input.fq2} "
        "{output.r1} {output.r1_unpaired} "
        "{output.r2} {output.r2_unpaired} "
        "{config[params][trimmomatic]} 2> {log}"

# ----------- 质控2 -----------

rule fastqc_posttrim:
    input:
        r1 = rules.trim.output.r1,  # 来自trim规则的clean输出
        r2 = rules.trim.output.r2
    output:
        html1 = "results/fastqc/trimmed/{sample}_R1_fastqc.html",  # R1报告
        zip1 = "results/fastqc/trimmed/{sample}_R1_fastqc.zip",    # R1数据
        html2 = "results/fastqc/trimmed/{sample}_R2_fastqc.html",  # R2报告  
        zip2 = "results/fastqc/trimmed/{sample}_R2_fastqc.zip"     # R2数据
    log:
        "logs/fastqc/trimmed/{sample}.log"  # 合并日志
    threads: 2
    conda:
        "envs/rnaseq.yml"
    shell:
        """
        # 运行FastQC（自动生成*_clean_fastqc.*）
        fastqc {input.r1} {input.r2} -o results/fastqc/trimmed/ 2> {log} && \
        # 重命名R1文件
        mv results/fastqc/trimmed/{wildcards.sample}_R1_clean_fastqc.html {output.html1} && \
        mv results/fastqc/trimmed/{wildcards.sample}_R1_clean_fastqc.zip {output.zip1} && \
        # 重命名R2文件
        mv results/fastqc/trimmed/{wildcards.sample}_R2_clean_fastqc.html {output.html2} && \
        mv results/fastqc/trimmed/{wildcards.sample}_R2_clean_fastqc.zip {output.zip2}
        """

# ----------- 质控3 -----------
rule multiqc:
    input:
        raw = expand("results/fastqc/raw/{sample}_R{read}_fastqc.zip", sample=samples.index, read=[1, 2]),
        trimmed = expand("results/fastqc/trimmed/{sample}_R{read}_fastqc.zip", sample=samples.index, read=[1, 2])
    output:
        "results/multiqc_report.html"
    conda:
        "envs/rnaseq.yml"
    shell:
        "multiqc results/fastqc/ -o results/"  # 自动合并raw和trimmed目录
# ----------- 比对与定量 -----------
rule hisat2_align:
    input:
        r1 = rules.trim.output.r1,
        r2 = rules.trim.output.r2
    output:
        bam = "results/align/{sample}.bam",
        sam = temp("results/align/{sample}.sam")
    params:
        index = config["genome"]["index"],  # 必须完整路径
        tmpdir = "tmp/{sample}"
    log:
        hisat2 = "logs/hisat2/{sample}.log",  # 定义 HISAT2 的日志文件
        samtools = "logs/samtools/{sample}.log"  # 定义 Samtools 的日志文件		
    threads: 4
    resources:
        mem_mb=16000
    shell:
        """
        # 创建临时目录
        mkdir -p {params.tmpdir}

        # HISAT2比对（日志重定向到单独文件）
        hisat2 -x {params.index} \
               -1 {input.r1} \
               -2 {input.r2} \
               --dta --threads {threads} \
               -S {output.sam} 2> {log.hisat2}

        # SAM转BAM并排序
        samtools view -bS {output.sam} | samtools sort -m 2G -T {params.tmpdir} -o {output.bam} 2> {log.samtools}


        # 清理
        rm -rf {params.tmpdir}
        """

rule featurecounts:
    input:
        bams = expand("results/align/{sample}.bam", sample=samples.index),
        gtf = config["genome"]["gtf"]
    output:
        "results/counts/gene_counts.tsv"
    log:
        "logs/featurecounts.log"
    conda:
        "envs/rnaseq.yml"
    shell:
        "featureCounts -a {input.gtf} -o {output} "
        "-T {threads} -p -t exon -g gene_id "
        "results/align/*.bam 2> {log}"

# ----------- 差异分析 -----------
rule deseq2:
    input:
        counts = "results/counts/gene_counts.tsv",
        samples = config["samples"],
        config = "config/config.json"  # 将配置文件作为输入
    output:
        "results/deseq2/DE_results.csv"
    conda:  
        "/home/wangwen/Biosoft/ENTER/envs/r-env"
    script:
        "scripts/deseq2.R"