#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_enriched_kmers.R target.fa background.fa output.csv")
}

target_file <- args[1]
background_file <- args[2]
output_file <- args[3]

# Read FASTA sequences
target_seqs <- readRNAStringSet(target_file)
background_seqs <- readRNAStringSet(background_file)

# Enrichment function
enriched_kmers <- function(target_seqs, background_seqs, k = 6) {
  target_kmer_counts <- oligonucleotideFrequency(target_seqs, width = k, step = 1)
  target_counts_sum <- colSums(target_kmer_counts)
  target_kmer_tot <- sum(target_counts_sum)

  bg_kmer_counts <- oligonucleotideFrequency(background_seqs, width = k, step = 1)
  bg_counts_sum <- colSums(bg_kmer_counts)
  bg_kmer_tot <- sum(bg_counts_sum)

  all_kmers <- union(names(target_counts_sum), names(bg_counts_sum))

  df_kmers <- data.frame(
    kmer = all_kmers,
    target = as.integer(target_counts_sum[all_kmers]),
    background = as.integer(bg_counts_sum[all_kmers]),
    stringsAsFactors = FALSE
  )

  df_kmers[is.na(df_kmers)] <- 0

  df_kmers$pvalue <- sapply(seq_len(nrow(df_kmers)), function(i) {
    t <- df_kmers$target[i]
    b <- df_kmers$background[i]
    mat <- matrix(c(t, target_kmer_tot - t, b, bg_kmer_tot - b), nrow = 2)
    fisher.test(mat, alternative = "two.sided")$p.value
  })

  df_kmers$padj <- p.adjust(df_kmers$pvalue, method = "BH")
  df_kmers <- df_kmers[order(df_kmers$padj), ]

  return(df_kmers)
}

# Run for multiple k values
k_values <- 5:10
all_results <- lapply(k_values, function(k) {
  res <- enriched_kmers(target_seqs, background_seqs, k = k)
  res$k <- k
  res
})
combined_kmers <- do.call(rbind, all_results)

fwrite(combined_kmers, output_file, row.names = FALSE)
