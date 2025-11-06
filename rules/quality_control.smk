rule raw_fastqc:
    input: raw_fastq = "raw_fastq/{sample}_{read_pair_tag}.fastq.gz"
    output: html="qc_reports/{sample}/raw_fastqc/{sample}_{read_pair_tag}_fastqc.html"
    log: "logs/{sample}/raw_fastqc_{read_pair_tag}.log"
    threads: 2
    conda: "../wrappers/raw_fastqc/env.yaml"
    script: "../wrappers/raw_fastqc/script.py"

def processed_fastq_input(wildcards):
    preprocessed = "processed_fastq"
    if read_pair_tags == ["SE"]:
        return os.path.join(preprocessed, "{sample}_R1.fastq.gz")
    else:
        return os.path.join(preprocessed, "{sample}_{read_pair_tag}.fastq.gz")

rule processed_fastqc:
    input: fastq = processed_fastq_input,
    output: html = "qc_reports/{sample}/processed_fastqc/{sample}_{read_pair_tag}_fastqc.html",
    log: "logs/{sample}/processed_fastqc_{read_pair_tag}.log"
    threads: 2
    conda: "../wrappers/processed_fastqc/env.yaml"
    script: "../wrappers/processed_fastqc/script.py"

def merge_raw_qc_input(wildcards):
    if config["long_reads"]:
        return directory(expand("qc_reports/{sample}/raw_nanoplot", sample = sample_name))
    else:
        return expand("qc_reports/{sample}/raw_fastqc/{sample}_{read_pair_tag}_fastqc.html", sample=sample_name, read_pair_tag=read_pair_tags)

def merge_processed_qc_input(wildcards):
    if config["long_reads"]:
        return directory(expand("qc_reports/{sample}/processed_nanoplot", sample = sample_name))
    else:
        return expand("qc_reports/{sample}/processed_fastqc/{sample}_{read_pair_tag}_fastqc.html", sample=sample_name, read_pair_tag=read_pair_tags)

rule merge_raw_fastqc_qc:
    input: html = merge_raw_qc_input
    output: report = "qc_reports/raw_multiqc_report.html"
    params: long_reads = config["long_reads"]
    log: "logs/all_samples/merge_raw_fastq_qc.log"
    conda: "../wrappers/merge_raw_qc/env.yaml"
    script: "../wrappers/merge_raw_qc/script.py"

rule processed_nanoplot:
    input: fastq= "processed_fastq/{sample}_R1.fastq.gz",
    output: out = directory("qc_reports/{sample}/processed_nanoplot")
    log: "logs/{sample}/processed_nanoplot.log"
    params: sample = "{sample}"
    threads: 2
    conda: "../wrappers/nanoplot/env.yaml"
    script: "../wrappers/nanoplot/script.py"

rule raw_nanoplot:
    input: fastq = "raw_fastq/{sample}_R1.fastq.gz"
    output: out =  directory("qc_reports/{sample}/raw_nanoplot")
    log: "logs/{sample}/raw_nanoplot.log"
    threads: 2
    params: sample = "{sample}"
    conda: "../wrappers/nanoplot/env.yaml"
    script: "../wrappers/nanoplot/script.py"

rule merge_processed_fastqc_qc:
    input: merge_processed_qc_input
    output: report = "qc_reports/processed_multiqc_report.html"
    params: long_reads = config["long_reads"]
    log: "logs/all_samples/merge_processed_fastq_qc.log"
    conda: "../wrappers/merge_processed_qc/env.yaml"
    script: "../wrappers/merge_processed_qc/script.py"

