#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(phyloseq)
})

read_emu <- function(file) {
  sample_id <- sub("_.*$", "", basename(file))   
  dt <- fread(file)

  if ("abundance_est" %in% names(dt)) {
    dt[, Count := as.numeric(abundance_est)]
  } else if ("num_reads" %in% names(dt)) {
    dt[, Count := as.numeric(num_reads)]
  } else if ("abundance" %in% names(dt)) {
    dt[, Count := round(as.numeric(abundance) * 1e6, 0)]
  } else {
    stop("No abundance columns in: ", file)
  }

  dt <- dt[!is.na(tax_id)]
  dt[, tax_id := as.character(tax_id)]
  dt[, Sample := sample_id]

  ranks <- c("superkingdom","phylum","class","order","family","genus","species")
  miss <- setdiff(ranks, names(dt))
  if (length(miss) > 0) for (cc in miss) dt[, (cc) := NA_character_]

  dt[, c("tax_id","superkingdom","phylum","class","order","family","genus","species","Count","Sample"), with=FALSE]
}

run_all <- function(metadata, output_rds, emu_files) {
  if (length(emu_files) == 0) stop("No EMU files (*.tsv) provided.")

  message("Reading EMU: ", length(emu_files), " file(s)")
  emu_list <- lapply(emu_files, read_emu)
  emu_all  <- rbindlist(emu_list, use.names = TRUE, fill = TRUE)

  ranks <- c("superkingdom","phylum","class","order","family","genus","species")
  otu_wide <- data.table::dcast(
    emu_all,
    tax_id + superkingdom + phylum + class + order + family + genus + species ~ Sample,
    value.var = "Count", fill = 0
  )

  sample_cols <- setdiff(names(otu_wide), c("tax_id", ranks))
  if (length(sample_cols) == 0) {
    stop("No sample columns detected in dcast(). Check Sample IDs.")
  }

  otu_mat_taxa_rows <- as.matrix(otu_wide[, ..sample_cols])
  storage.mode(otu_mat_taxa_rows) <- "numeric"
  rownames(otu_mat_taxa_rows) <- otu_wide$tax_id

  keep_taxa <- rowSums(otu_mat_taxa_rows) > 0
  if (!any(keep_taxa)) {
    stop("All taxa have zero counts. Check the abundance column used (abundance_est/num_reads/abundance).")
  }
  otu_mat_taxa_rows <- otu_mat_taxa_rows[keep_taxa, , drop=FALSE]

  seqtab.nochim <- t(otu_mat_taxa_rows)

  # --- Metadata ---
  if (!is.null(metadata) && nzchar(metadata) && file.exists(metadata)) {
    exp_des <- as.data.frame(fread(metadata))
    if (!("SampleID" %in% names(exp_des))) {
      names(exp_des)[1] <- "SampleID"
    }
    rownames(exp_des) <- as.character(exp_des$SampleID)
  } else {
    exp_des <- data.frame(SampleID = rownames(seqtab.nochim),
                          row.names = rownames(seqtab.nochim),
                          stringsAsFactors = FALSE)
  }

  common <- intersect(rownames(seqtab.nochim), rownames(exp_des))
  if (length(common) == 0) {
    stop("sample_name column in metadata does not match sample names in EMU files.\n",
         "Ej (EMU): ", paste(head(rownames(seqtab.nochim)), collapse=", "), "\n",
         "Ej (META): ", if (exists("exp_des")) paste(head(rownames(exp_des)), collapse=", ") else "sin metadata")
  }
  seqtab.nochim <- seqtab.nochim[common, , drop=FALSE]
  exp_des <- exp_des[common, , drop=FALSE]

  if (nrow(seqtab.nochim) == 0 || ncol(seqtab.nochim) == 0) {
    stop("Empty OTU table after filtering.")
  }

  tax_tab <- as.data.table(otu_wide[keep_taxa, c("tax_id", ranks), with=FALSE])
  tax_mat <- as.matrix(tax_tab[, ranks, with=FALSE])
  storage.mode(tax_mat) <- "character"
  rownames(tax_mat) <- tax_tab$tax_id

  colnames(tax_mat) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  empties <- !is.na(tax_mat) & trimws(tax_mat) == ""
  if (any(empties)) tax_mat[empties] <- NA_character_

  ps <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    sample_data(exp_des),
    tax_table(tax_mat)
  )

  taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

  dir.create(dirname(output_rds), showWarnings = FALSE, recursive = TRUE)
  saveRDS(ps, file = output_rds)
  message("OK: phyloseq saved at: ", normalizePath(output_rds))


}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript emu_to_phyloseq.R <metadata.csv> <output.rds> <EMU1.tsv> [EMU2.tsv ...]")
}
metadata <- args[1]
output   <- args[2]
emu_files <- args[-c(1,2)]

run_all(metadata, output, emu_files)