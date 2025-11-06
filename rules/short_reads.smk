rule kraken2_short:
    input:  paired = expand("processed_fastq/{{sample}}_{read_pair_tag}.fastq.gz", read_pair_tag = read_pair_tags),
            database = config["kraken2_db"]
    output: kraken_report = "results/kraken/{sample}_kraken2_report.txt",
            classified = "results/kraken/{sample}_kraken2_output.txt"
    log: "logs/{sample}/kraken2_short.log"
    params: is_paired = config["is_paired"],
            conf_interval = config["kraken2_conf_interval"],
            min_base_qual = config["kraken2_min_base_qual"],
    threads: 2
    conda: "../wrappers/kraken2/env.yaml"
    script: "../wrappers/kraken2/script.py"

rule dada2:
    input: fastq = expand("processed_fastq/{sample}_{read_pair_tag}.fastq.gz", sample = sample_name, read_pair_tag = read_pair_tags),
           train_set = config["dada2_train_set"],
           sp_assign = config["dada2_sp_assign"],
           experiment_design = "experiment_design/metadata.tsv"
    output: summary = "results/dada2/summary.tsv"
    log: "logs/all_samples/dada2/dada2.log"
    params: for_trunc = config["dada2_for_trunc"],
            rev_trunc = config["dada2_rev_trunc"],
            multithread = config["dada2_multithread"],
            for_error = config["dada2_for_error"],
            rev_error = config["dada2_rev_error"]
    threads: 10
    conda: "../wrappers/dada2/env.yaml"
    script: "../wrappers/dada2/script.py"


rule centrifuge:
    input: processed = expand("processed_fastq/{{sample}}_{read_pair_tag}.fastq.gz", read_pair_tag = read_pair_tags)
    output: classified = "results/centrifuger/{sample}_centrifuger.tsv",
            kreport = "results/centrifuger/{sample}_kraken2_report.txt"
    log: "logs/{sample}/centrifuger.log"
    threads: 30
    params:
        is_paired = config["is_paired"],
        centrifuge_db = config["centrifuger_db_path"],
        min_hitlen = config["centrifuger_min_hitlen"],
        max_hits = config["centrifuger_max_hits"]
    conda:  "../wrappers/centrifuger/env.yaml"
    script: "../wrappers/centrifuger/script.py"