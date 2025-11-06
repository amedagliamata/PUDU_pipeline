if config["long_reads"]:
    rule nanofilt:
        input: fastq = "raw_fastq/{sample}_R1.fastq.gz"
        output: processed = "processed_fastq/{sample}_R1.fastq.gz"
        log: "logs/{sample}/preprocessing.log"
        threads: 2
        params:
            length = config["nanofilt_length"],
            quality = config["nanofilt_quality"], 
            maxlength = config["nanofilt_maxlength"],
            mingc = config["nanofilt_mingc"],
            maxgc = config["nanofilt_maxgc"],
            headcrop = config["nanofilt_headcrop"],
            tailcrop = config["nanofilt_tailcrop"]
        conda:  "../wrappers/nanofilt/env.yaml"
        script: "../wrappers/nanofilt/script.py"
else:
    rule trimmomatic:
        input: fastq = expand("raw_fastq/{{sample}}_{read_pair_tag}.fastq.gz", read_pair_tag = read_pair_tags)
        output: paired = expand("processed_fastq/{{sample}}_{read_pair_tag}.fastq.gz", read_pair_tag = read_pair_tags)
        log: "logs/{sample}/preprocessing.log"
        threads: 2
        params:
            trim_adapters = config["trim_adapters"],
            qc_trim_reads = config["trimmomatic_qc_trim_reads"],
            adapters = config["trimmomatic_adapters"],
            leading = config["trimmomatic_leading"],
            trailing = config["trimmomatic_trailing"],
            slidingwindow = config["trimmomatic_slidingwindow"],
            minlen = config["trimmomatic_minlen"]
        conda:  "../wrappers/trimmomatic/env.yaml"
        script: "../wrappers/trimmomatic/script.py"