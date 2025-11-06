library(dada2)
library(DECIPHER)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(data.table)
library(tibble)


run_all <- function(experiment_design, for_trunc, rev_trunc, multithread, train_set, sp_assign, for_error, rev_error, files_list) {


    print(files_list)

    exp_des <- fread(experiment_design)
    exp_des <- column_to_rownames(exp_des, var = "sample_name")

    print("Loaded experiment design")


    fnFs <- sort(grep("_R1.fastq.gz$", files_list, value = TRUE))
    fnRs <- sort(grep("_R2.fastq.gz$", files_list, value = TRUE))
    
    print("Forward and reverse files ready")

    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    
    filtFs <- file.path(getwd(), "results/dada2/", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(getwd(), "results/dada2/", paste0(sample.names, "_R_filt.fastq.gz"))

    print("Filtered path")

    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                         truncLen = c(as.numeric(for_trunc),as.numeric(rev_trunc)), maxN = 0,
                         maxEE = c(as.numeric(for_error), as.numeric(rev_error)), truncQ = 2, rm.phix = TRUE,
                         compress = TRUE, multithread = multithread)

    print("Filtered")
    
    errF <- learnErrors(filtFs, multithread = multithread)
    errR <- learnErrors(filtRs, multithread = multithread)

    print("Learned errors")

    dadaFs <- dada(filtFs, err = errF, multithread = multithread)
    dadaRs <- dada(filtRs, err = errR, multithread = multithread)

    print("Denoised")

    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                          minOverlap = 12, maxMismatch = 1,
                          verbose = TRUE)

    print("merged")


    seqtab <- makeSequenceTable(mergers)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                        multithread = multithread,
                                        verbose = TRUE)
    print(rownames(seqtab.nochim))

    print("chimeras removed")

    sum(seqtab.nochim)/sum(seqtab)

    getN <- function(x) sum(getUniques(x))
    track <- as.data.table(cbind(out, sapply(dadaFs, getN),
                       sapply(dadaRs, getN),
                       sapply(mergers, getN),
                       rowSums(seqtab.nochim)))
    colnames(track) <- c("input", "filtered", "denoisedF",
                         "denoisedR", "merged", "nonchim")
    track$sample_name <- sample.names
    fwrite(track, paste0(getwd(), "/results/dada2/summary.tsv"))

    taxa <- assignTaxonomy(seqtab.nochim,
                           train_set,
                           multithread = multithread)

    print("taxonomy assigned")

    taxa <- addSpecies(taxa, sp_assign)

    exp_des <- as.data.frame(exp_des[match(rownames(seqtab.nochim), rownames(exp_des)), , drop = FALSE])
    
    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                   sample_data(exp_des),
                   tax_table(taxa))
    
    print("Sequence table rownames:")
    print(rownames(seqtab.nochim))
    print("Metadata rownames:")
    print(rownames(exp_des))
    print("phyloseq object created")

    dna <- DNAStringSet(taxa_names(ps))
    names(dna) <- taxa_names(ps)
    ps <- merge_phyloseq(ps, dna)
    taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))


    saveRDS(ps, file = paste0(getwd(), "/results/dada2/phyloseq.rds"))


}


script.dir <- dirname(gsub("--file=","",commandArgs()[grep("--file",commandArgs())]))


args <- commandArgs(trailingOnly = TRUE)

experiment_design <- args[1]
for_trunc <- as.numeric(args[2])
rev_trunc <- as.numeric(args[3])
multithread <- as.logical(toupper(args[4]))
train_set <- args[5]
sp_assign <- args[6]
for_error <- as.numeric(args[7])
rev_error <- as.numeric(args[8])
files_list <- tail(args, -8)

run_all(experiment_design, for_trunc, rev_trunc, multithread, train_set, sp_assign, for_error, rev_error, files_list)

