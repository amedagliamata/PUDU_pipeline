# AGENT.md

Guidance for AI coding agents (Claude Code, Copilot, etc.) working in this repository.

## What this is

PUDU (Pipeline for Universal Diversity Unveiling) is a **Snakemake** bioinformatics pipeline for metagenomic and amplicon sequencing analysis. It supports Illumina short reads and Oxford Nanopore/PacBio long reads, running them through QC, preprocessing, and taxonomic classification (Kraken2, Centrifuger, DADA2, EMU).

See [README.md](README.md) for user-facing documentation (features, config options, output layout, troubleshooting). Read it before making changes - this file only covers things the README doesn't: how the code is structured and how to work on it safely.

## Repository layout

```
Snakefile                  # entry point: computes `all_input`, includes rule files
config.yaml                # all user-tunable parameters (read via `config["..."]`)
rules/*.smk                # Snakemake rule definitions, grouped by pipeline stage
wrappers/<tool>/            # one directory per tool
  env.yaml                  #   conda environment for that rule
  script.py                 #   Python glue invoked via `script:` (uses `snakemake.*` object)
  *.R                       #   R scripts for R-based tools (dada2, emu phyloseq conversion)
adapters/*.fa               # Trimmomatic adapter FASTAs
raw_fastq/                  # user-supplied input (not versioned)
experiment_design/metadata.tsv  # sample metadata for DADA2/EMU (not versioned)
qc_reports/, processed_fastq/, results/, logs/  # pipeline outputs (not versioned)
```

Rule files, by pipeline stage:
- `rules/quality_control.smk` — FastQC/NanoPlot/MultiQC on raw and processed reads
- `rules/preprocessing.smk` — Trimmomatic (short reads) / NanoFilt (long reads)
- `rules/short_reads.smk` — Kraken2, Centrifuger, DADA2 for short reads
- `rules/long_reads.smk` — Kraken2, Centrifuger, EMU for long reads
- `rules/postprocessing.smk` — Bracken, Krona graphs, rarefaction curves, OTU tables

## Conventions

- **Every wrapper follows the same shape**: a `script.py` that reads `snakemake.input/output/params/log/threads`, writes a header to the log file, logs the resolved conda environment (`conda list`), builds a shell command string, appends it to the log, then runs it via `shell(command)`. Follow this pattern exactly when adding a new tool wrapper - don't introduce a different invocation style (e.g. `subprocess.run` without logging) without a reason.
- **One conda env per tool**, declared in `wrappers/<tool>/env.yaml` and referenced from the rule via `conda: "../wrappers/<tool>/env.yaml"`. Never bundle multiple tools' dependencies into one env file.
- **All tunables live in `config.yaml`**, accessed in rules as `config["key_name"]`, never hardcoded in a rule or wrapper. When adding a parameter, add it to `config.yaml` with a comment explaining units/defaults, matching the existing section headers (`###...###` banners).
- **Wildcards** (`sample`, `tool`, `read_pair_tag`, `taxlvl`) are constrained in the `Snakefile` via `wildcard_constraints`. If you add a new tool or taxonomic level option, update the relevant constraint there too.
- **Validation lives in the Snakefile**, not buried in rules: incompatible config combinations (e.g. long reads + paired, DADA2 + long reads, EMU + short reads) raise `ValueError`/`warnings.warn` early in `Snakefile`, before any rule runs. Add new cross-parameter validation there.
- **Output paths are stage-based**, not tool-based, except under `results/<tool>/`: `qc_reports/`, `processed_fastq/`, `results/<tool>/`, `logs/<sample or all_samples>/`. Keep new rules consistent with this layout so `rule all` / `all_input()` in the `Snakefile` can assemble them via `expand(...)`.
- **`rule all`'s inputs are computed dynamically** in `all_input()` in the `Snakefile`, gated by `config["analysis_type"]` and the `use_<tool>` flags. When adding a new tool, add its outputs there behind its own `use_<tool>` flag rather than always requesting them.

## Databases and metadata are not managed by the pipeline

Kraken2, Centrifuger, DADA2, and EMU all require pre-built databases that the user downloads/builds *outside* this pipeline and points to via paths in `config.yaml` (`kraken2_db`, `centrifuger_db_path`, `dada2_train_set`/`dada2_sp_assign`, `emu_db_path`). No wired-in rule downloads or builds a database. `wrappers/kraken2_reference/script.py` (a `kraken2-build`/`bracken-build` wrapper) exists but is not referenced by any rule in `rules/*.smk` or the `Snakefile` — treat it as dead code, not a working "build the DB" entry point, unless you're the one wiring it in.

Likewise, `experiment_design/metadata.tsv` (a hardcoded input path used by the `dada2` rule in `rules/short_reads.smk` and the EMU-to-phyloseq step in `rules/long_reads.smk`) is real sample metadata the user must author, it's matched against sample IDs during DADA2's phyloseq assembly (`wrappers/dada2/filtering_step.R`). An agent has no basis for generating its contents (condition/group labels, etc.); if it's missing, ask the user for it or point them at the format shown in the README rather than fabricating values.

## Making changes

- Adding a new classification/QC tool: create `wrappers/<tool>/{env.yaml,script.py}`, add a rule in the appropriate `rules/*.smk` file, add a `use_<tool>` flag and its parameters to `config.yaml`, wire its outputs into `all_input()` in the `Snakefile`, and add it to the wildcard `tools` list logic if it participates in Bracken/Krona/rarefaction postprocessing.
- Keep short-read and long-read logic separate (`short_reads.smk` vs `long_reads.smk`) even when a tool (e.g. Kraken2) supports both — this mirrors the existing structure and the `is_paired`/`long_reads` config-driven branching.
- Never commit real data: `raw_fastq/`, `experiment_design/metadata.tsv`, databases, and pipeline outputs (`qc_reports/`, `processed_fastq/`, `results/`, `logs/`) are user-provided/generated and should stay out of version control.

## Validating changes

There's no unit test suite — this is a Snakemake workflow, so validation means running it.

```bash
# Static check: confirm the DAG builds without executing anything
snakemake -n --use-conda

# Visualize the DAG for a sanity check on rule wiring
snakemake --dag | dot -Tsvg > dag.svg

# Full run (requires conda/mamba and configured databases)
snakemake --use-conda --cores 8
```

When editing a single wrapper, prefer targeting its specific output file rather than running the whole pipeline, e.g.:

```bash
snakemake --use-conda --cores 4 results/kraken/<sample>_kraken2_report.txt
```

Always check `config.yaml` combinations relevant to your change (e.g. `is_paired`, `long_reads`, `use_<tool>`) - the `Snakefile` will raise on invalid combinations, which is expected behavior, not a bug to work around.
