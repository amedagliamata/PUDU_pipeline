dada2_filtering <- function(data_dir, work_dir, trunclen, multithread, train_set, sp_assign, experiment_design) {

    fnFs <- sort(list.files(data_dir,
                            pattern = "_R1.fastq.gz",
                            full.names = TRUE))
    
    fnRs <- sort(list.files(data_dir,
                            pattern = "_R2.fastq.gz",
                            full.names = TRUE))

    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    
    filtFs <- file.path(work_dir, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(work_dir, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

    names(filtFs) <- sample.names
    names(filtRs) <- sample.names

    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                         truncLen = c(240,160), maxN = 0,
                         maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE,
                         compress = TRUE, multithread = F)
    
    errF <- learnErrors(filtFs, multithread = F)
    errR <- learnErrors(filtRs, multithread = F)

    dadaFs <- dada(filtFs, err = errF, multithread = F)
    dadaRs <- dada(filtRs, err = errR, multithread = F)

    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs,
                          minOverlap = 20, maxMismatch = 2,
                          verbose = TRUE)

    seqtab <- makeSequenceTable(mergers)
    seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                        multithread = multithread,
                                        verbose = TRUE)
    sum(seqtab.nochim)/sum(seqtab)

    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, sapply(dadaFs, getN),
                   sapply(dadaRs, getN),
                   sapply(mergers, getN),
                   rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF",
                         "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    fwrite(track, paste0(work_dir, "/results/dada2/summary.tsv"))

    taxa <- assignTaxonomy(seqtab.nochim,
                           train_set,
                           multithread = multithread)

    taxa <- addSpecies(taxa, sp_assign)

    experiment_design <- experiment_design[match(rownames(seqtab.nochim), rownames(experiment_design)), ]

    ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
                   sample_data(experiment_design),
                   tax_table(taxa))

    dna <- DNAStringSet(taxa_names(ps))
    names(dna) <- taxa_names(ps)
    ps <- merge_phyloseq(ps, dna)
    taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

    return(ps)

}