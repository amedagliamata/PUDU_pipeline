rule emu_analysis:
    input: processed = "processed_fastq/{sample}_R1.fastq.gz",
           db = config["emu_db_path"] 
    output: result = "results/emu/{sample}_rel-abundance.tsv"
    log: "logs/{sample}/emu.log"
    threads: 30
    params: tax_lvl = config["organism_taxonomic_level"]
    conda: "../wrappers/emu/env.yaml"
    script: "../wrappers/emu/script.py"

rule emu_to_phyloseq:
       input: emu_files = expand("results/emu/{sample}_rel-abundance.tsv", sample=sample_name),
              experiment_design = "experiment_design/metadata.tsv"
       output: rds = "results/emu/emu_phyloseq.rds"
       log: "logs/all_samples/emu/emu_to_phyloseq.log"
       conda: "../wrappers/emu/env.yaml"
       script: "../wrappers/emu/script_to_phyloseq.py"