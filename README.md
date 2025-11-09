# PUDU Pipeline

**P**ipeline for **U**niversal **D**iversity **U**nveiling

A comprehensive Snakemake-based bioinformatics pipeline for metagenomic and amplicon sequencing analysis, supporting both Illumina short reads and long reads.

Its name refers to the pudú (*Pudu puda*), the smallest deer in the world, native to the temperate forests of southern Chile. Just as the animal occupies a small ecological footprint while maintaining high adaptability, PUDU is intentionally designed as a compact and modular workflow, providing  lightweight but effective framework for microbiome and metagenomic analysis

## Overview

PUDU Pipeline is a flexible and modular workflow designed for taxonomic classification and quality control of sequencing data. It integrates multiple state-of-the-art tools and provides automated quality control, preprocessing, and taxonomic profiling for both short-read (Illumina) and long-read sequencing technologies.

### Key Features

- **Multi-platform support**: Works with Illumina (paired-end/single-end) and long reads
- **Multiple classifiers**: Kraken2, Centrifuger, DADA2, and EMU
- **Comprehensive QC**: FastQC, NanoPlot, and MultiQC reporting
- **Flexible preprocessing**: Trimmomatic for short reads, NanoFilt for long reads
- **Visualization**: Krona graphs and rarefaction curves
- **OTU table generation**: At any taxonomic level (Domain to Species)
- **Modular design**: Enable/disable tools as needed

## Supported Tools

### Quality Control
- **FastQC**: Quality assessment for Illumina reads
- **NanoPlot**: Quality assessment for long reads
- **MultiQC**: Aggregated quality control reports

### Preprocessing
- **Trimmomatic**: Adapter trimming and quality filtering for Illumina reads
- **NanoFilt**: Quality and length filtering for long reads

### Taxonomic Classification
- **Kraken2**: Fast k-mer based metagenomic classifier with Bracken abundance estimation
- **Centrifuger**: Compressed FM-index classifier (replacement for Centrifuge)
- **DADA2**: Amplicon sequence variant (ASV) inference for paired-end Illumina data
- **EMU**: Expectation-Maximization algorithm for 16S targeted long-read sequencing

### Visualization & Analysis
- **Krona**: Interactive taxonomic visualization
- **Rarefaction curves**: Assess sampling depth adequacy
- **OTU tables**: Taxonomic abundance matrices at custom levels
- **Phyloseq**: R objects for downstream analysis (EMU)

## Requirements

### Software Dependencies
- **Snakemake** ≥ 5.18.0
- **Conda/Mamba** (for environment management)

All tool-specific dependencies are managed through conda environments defined in `wrappers/*/env.yaml` files.

### Database Requirements

Depending on which tools you enable, you'll need:

- **Kraken2**: Pre-built Kraken2 database (e.g., Standard, PlusPF, PlusPFP)
- **Centrifuger**: Centrifuger index files (.cf format)
- **DADA2**: SILVA, GreenGenes, or RDP training sets
- **EMU**: EMU-specific database with taxonomy.tsv and species_taxid.fasta

## Installation

```bash
# Clone the repository
git clone https://github.com/amedagliamata/PUDU_pipeline.git
cd PUDU_pipeline

# Create conda environment with Snakemake
conda create -n pudu snakemake>=5.18.0
conda activate pudu
```

## Quick Start

### 1. Prepare Your Data

Place raw FASTQ files in the `raw_fastq/` directory:

**For paired-end Illumina:**
```
raw_fastq/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

**For single-end or long reads:**
```
raw_fastq/
├── sample1_R1.fastq.gz
└── sample2_R1.fastq.gz
```

### 2. Create Metadata File (Optional, for DADA2/EMU)

Create `experiment_design/metadata.tsv` with sample information:
```tsv
sample_name	condition	group
sample1	control	A
sample2	treatment	B
```

### 3. Configure the Pipeline

Edit `config.yaml` to specify:

```yaml
# Basic configuration
sample_name: ["sample1", "sample2"]
analysis_type: "full"  # or "raw_qc" for QC only

# Sequencing type
is_paired: true        # false for single-end or ONT
long_reads: false      # true for long reads

# Taxonomic level for OTU tables
organism_taxonomic_level: "G"  # D, P, C, O, F, G, or S

# Enable/disable tools
use_kraken: true
use_centrifuger: false
use_dada2: false
use_emu: false

# Database paths
kraken2_db: "/path/to/kraken2/db"
# ... other tool-specific parameters
```

### 4. Run the Pipeline


**Execute the pipeline:**
```bash
# Using conda environments
snakemake --use-conda --cores 8

# With mamba (faster)
snakemake --use-conda --conda-frontend mamba --cores 8
```

## Configuration Guide

### Analysis Types

- **`full`**: Complete analysis including QC, trimming, and classification
- **`raw_qc`**: Only quality control on raw reads (no classification)

### Sequencing Technology Settings

| Sequencing Type | `is_paired` | `long_reads` |
|----------------|-------------|--------------|
| Illumina Paired-End | `true` | `false` |
| Illumina Single-End | `false` | `false` |
| Long reads | `false` | `true` |

### Taxonomic Levels

Choose from: `D` (Domain), `P` (Phylum), `C` (Class), `O` (Order), `F` (Family), `G` (Genus), `S` (Species)

### Tool-Specific Configuration

#### Trimmomatic (Illumina)
```yaml
trim_adapters: true
trimmomatic_qc_trim_reads: true
trimmomatic_adapters: "NexteraPE-PE"  # or path to custom adapters
trimmomatic_leading: 3
trimmomatic_trailing: 3
trimmomatic_slidingwindow: "4:15"
trimmomatic_minlen: 36
```

#### NanoFilt (Long reads)
```yaml
nanofilt_length: 50        # Minimum read length
nanofilt_quality: 10       # Minimum average quality
nanofilt_maxlength: 99999999
nanofilt_headcrop: 0       # Trim from start
nanofilt_tailcrop: 0       # Trim from end
```

#### Kraken2
```yaml
kraken2_db: "/path/to/database"
kraken2_conf_interval: 0.2  # Confidence threshold
kraken2_rlen: 100           # Read length for Bracken
kraken2_kmer: 35            # K-mer length
```

#### DADA2
```yaml
dada2_train_set: "/path/to/silva_train.fa"
dada2_sp_assign: "/path/to/silva_species.fa"
dada2_for_trunc: 280
dada2_rev_trunc: 260
dada2_multithread: true     # false on Windows
```

## Output Structure

```
PUDU_pipeline/
├── qc_reports/
│   ├── raw_multiqc_report.html
│   ├── processed_multiqc_report.html
│   └── {sample}/
│       ├── raw_fastqc/
│       └── processed_fastqc/
├── processed_fastq/
│   └── {sample}_{R1,R2}.fastq.gz
├── results/
│   ├── kraken/
│   │   ├── {sample}_kraken2_report.txt
│   │   ├── {sample}_kraken2_report_bracken.txt
│   │   ├── {sample}_krona_graph.html
│   │   ├── {sample}_rarefaction.jpg
│   │   └── otu_table_{taxlvl}.csv
│   ├── centrifuger/
│   │   ├── {sample}_centrifuger.tsv
│   │   ├── {sample}_kraken2_report_bracken.txt
│   │   ├── {sample}_krona_graph.html
│   │   ├── {sample}_rarefaction.jpg
│   │   └── otu_table_{taxlvl}.csv
│   ├── dada2/
│   │   └── summary.tsv
│   └── emu/
│       ├── {sample}_rel-abundance.tsv
│       └── emu_phyloseq.rds
└── logs/
```

## Workflow Examples

### Example 1: Illumina Paired-End Metagenomics

```yaml
sample_name: ["sample1", "sample2"]
analysis_type: "full"
is_paired: true
long_reads: false
use_kraken: true
use_centrifuger: true
trimmomatic_qc_trim_reads: true
kraken2_db: "/databases/kraken2/standard"
```

### Example 2: Long reads 16S Sequencing

```yaml
sample_name: ["ont_sample1", "ont_sample2"]
analysis_type: "full"
is_paired: false
long_reads: true
use_kraken: true
use_emu: true
organism_taxonomic_level: "S"
nanofilt_length: 1000
nanofilt_quality: 10
emu_db_path: "/databases/emu/silva_db"
```

### Example 3: Illumina Amplicon with DADA2

```yaml
sample_name: ["amplicon1", "amplicon2"]
analysis_type: "full"
is_paired: true
long_reads: false
use_dada2: true
trim_adapters: true
trimmomatic_adapters: "TruSeq3-PE"
dada2_train_set: "/databases/silva/silva_train.fa"
dada2_sp_assign: "/databases/silva/silva_species.fa"
```

## Troubleshooting

### Common Issues

1. **"Long reads are enabled but the reads are paired"**
   - Set `is_paired: false` for Long reads data

2. **"DADA2 is enabled but the reads are not paired"**
   - DADA2 requires paired-end reads; set `is_paired: true`

3. **"EMU is applied only to long-reads"**
   - EMU is designed for long reads; set `long_reads: true` or disable EMU

4. **Conda environment issues**
   - Use `--conda-frontend mamba` for faster environment resolution
   - Try `snakemake --use-conda --conda-cleanup-pkgs cache`

### Performance Optimization

- Increase `--cores` for parallel execution
- Use `--resources mem_mb=<value>` to limit memory usage
- Enable tool-specific multithreading (e.g., `dada2_multithread: true`)

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contact

For questions or support, please open an issue on GitHub or contact the repository maintainer.

## Acknowledgments

This pipeline integrates tools developed by the scientific community. We thank all developers and maintainers of the integrated software.
