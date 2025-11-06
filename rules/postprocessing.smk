rule bracken_abundance:
    input: kraken_report = "results/{tool}/{sample}_kraken2_report.txt",
           database = config["kraken2_db"]
    output: ab_report = "results/{tool}/{sample}_kraken2_report_bracken.txt"
    log: "logs/{sample}/bracken_{tool}.log"
    params: kmer_length = config["kraken2_kmer"],
            read_length = config["kraken2_rlen"]
    threads: 1
    conda: "../wrappers/bracken/env.yaml"
    script: "../wrappers/bracken/script.py"

rule rarefaction_plot:
  input: kraken_report = "results/{tool}/{sample}_kraken2_report.txt"
  output: rare_plot = "results/{tool}/{sample}_rarefaction.jpg"
  log: "logs/{sample}/rarefaction_{tool}.log"
  params: tax_lvl = taxonomic_level
  threads: 1
  conda: "../wrappers/rarefaction_plot/env.yaml"
  script: "../wrappers/rarefaction_plot/script.py"

rule krona_graph_kraken:
    input: kraken_report = "results/{tool}/{sample}_kraken2_report.txt",
           database_updated = "results/krona_updated_db.txt"
    output: krona_graph = "results/{tool}/{sample}_krona_graph.html"
    log: "logs/{sample}/krona_graph_{tool}.log"
    threads: 1
    conda: "../wrappers/krona_graph/env.yaml"
    script: "../wrappers/krona_graph/script.py"

rule krona_update_db:
    input: kraken_report = expand("results/{tool}/{sample}_kraken2_report.txt", tool = tools, sample = sample_name)
    output: "results/krona_updated_db.txt"
    log: "logs/all_samples/krona_updated_db/krona_updated.log"
    conda: "../wrappers/krona_graph/env.yaml"
    shell: "ktUpdateTaxonomy.sh; touch {output} >> {log}"

rule create_otu_table:
    input: kraken_reports = expand("results/{{tool}}/{sample}_kraken2_report.txt", sample = sample_name)
    output: otu_table = "results/{tool}/otu_table_{taxlvl}.csv"
    log: "logs/all_samples/otu_tables/{tool}/{taxlvl}_otu_table.log"
    params: lv = taxonomic_level
    threads: 1
    conda: "../wrappers/create_otu_table/env.yaml"
    script: "../wrappers/create_otu_table/script.py"