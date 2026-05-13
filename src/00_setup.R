# ==============================================================================
# 00_setup.R
# Shared configuration for Figure 5 analysis
#
# This script is sourced by every RMarkdown file in the project.
# It loads required packages, sets global options, defines reusable helper
# functions, and establishes color palettes used across all figures.
# ==============================================================================

# ---- Load packages -----------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(Biostrings)
  library(ggpubr)
  library(ggforce)
  library(biomaRt)
  library(msa)
  library(motifmatchr)
  library(UpSetR)
  library(stringdist)
  library(ggrepel)
  library(ggseqlogo)
  library(ggVennDiagram)
  library(eulerr)
  library(paletteer)
  library(DESeq2)
})

# ---- Global options ----------------------------------------------------------

options(scipen = 100)
theme_update(text = element_text(family = "Arial"))

# ---- File paths --------------------------------------------------------------

# All paths are relative to the project root.
# When sourcing from src/, set the working directory to the project root, or
# use here::here() for robust path resolution.

data_dir    <- file.path("..", "data")
results_dir <- file.path("..", "results")
figures_dir <- file.path(results_dir, "figures")

# Create output directories if they don't exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Color palettes ----------------------------------------------------------

#' Color palette for KD/control gene groups (g3, c1–c4, cth1)
colors_guide <- c(
  "c1"   = "#CC79A7",
  "c2"   = "#D991BA",
  "c3"   = "#E6A6C3",
  "c4"   = "#F2C1D1",
  "g3"   = "#1f78b4",
  "cth1" = "red"
)

#' Color palette for endogenous motif groups
highlight_colors <- c(
  "2A 8mer"            = "deepskyblue",
  "3A 8mer"            = "#0E34B2FF",
  "more than one site" = "purple3",
  "one site"           = "plum1",
  "Rest of genes"      = "gray40"
)

# ---- Helper functions --------------------------------------------------------

#' Wilcoxon rank-sum test against a baseline group
#'
#' For each non-baseline level of `class_col`, performs a two-sided Wilcoxon
#' test vs. the baseline group. Returns a data.frame with group-level summary
#' statistics, p-values, BH-adjusted p-values, and rank-biserial effect sizes.
#'
#' @param data       A data.frame.
#' @param class_col  Column name (string) containing group labels.
#' @param baseline   The reference group label.
#' @param value_col  Column name (string) containing the numeric values.
#' @param alternative Alternative hypothesis: "two.sided", "less", or "greater".
#' @return A data.frame with one row per non-baseline group.
wilcox_vs_baseline <- function(data,
                               class_col   = "group",
                               baseline    = "Rest of genes",
                               value_col   = "log2_fc",
                               alternative = "two.sided") {
  out <- list()
  for (g in setdiff(unique(data[[class_col]]), baseline)) {
    x <- data[data[[class_col]] == g, value_col, drop = TRUE]
    y <- data[data[[class_col]] == baseline, value_col, drop = TRUE]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]

    n1 <- length(x)
    n2 <- length(y)

    if (n1 > 0 && n2 > 0 && length(unique(c(x, y))) > 1) {
      wt  <- wilcox.test(x, y, alternative = alternative)
      W   <- unname(wt$statistic)
      U   <- W - n1 * (n1 + 1) / 2
      rrb <- (2 * U) / (n1 * n2) - 1
      out[[g]] <- data.frame(
        group               = g,
        n_group             = n1,
        n_base              = n2,
        med_group           = median(x),
        med_base            = median(y),
        stat_W              = W,
        p_value             = wt$p.value,
        effect_rank_biserial = rrb
      )
    } else {
      out[[g]] <- data.frame(
        group               = g,
        n_group             = n1,
        n_base              = n2,
        med_group           = if (n1 > 0) median(x) else NA,
        med_base            = if (n2 > 0) median(y) else NA,
        stat_W              = NA,
        p_value             = NA,
        effect_rank_biserial = NA
      )
    }
  }
  res       <- do.call(rbind, out)
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res
}


#' Batched biomaRt query to avoid Ensembl HTTP 500 errors
#'
#' Splits a large vector of IDs into smaller batches and queries biomaRt
#' sequentially, with a short pause between requests to avoid overloading
#' the server. This is especially necessary for heavy attributes like
#' `coding`, `3utr`, and `5utr` on archived Ensembl hosts.
#'
#' @param attributes Character vector of biomaRt attributes.
#' @param filters    Filter name (e.g., "ensembl_transcript_id").
#' @param values     Character vector of IDs to query.
#' @param mart       A biomaRt Mart object (from `useEnsembl()`).
#' @param batch_size Number of IDs per batch (default: 200).
#' @return A data.frame with all results combined.
batch_getBM <- function(attributes, filters, values, mart, batch_size = 200) {
  chunks  <- split(values, ceiling(seq_along(values) / batch_size))
  results <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    message("biomaRt batch ", i, " of ", length(chunks),
            " (", length(chunks[[i]]), " IDs)")
    results[[i]] <- tryCatch(
      getBM(attributes = attributes,
            filters    = filters,
            values     = chunks[[i]],
            mart       = mart),
      error = function(e) {
        message("  Batch ", i, " failed: ", conditionMessage(e),
                "\n  Retrying in 5 seconds...")
        Sys.sleep(5)
        getBM(attributes = attributes,
              filters    = filters,
              values     = chunks[[i]],
              mart       = mart)
      }
    )
    Sys.sleep(1)
  }

  bind_rows(results)
}


#' Find motif positions in 3'UTR sequences
#'
#' Searches for exact matches of RNA motifs in the 3'UTR sequences provided.
#' DNA sequences are internally converted to RNA (T → U) before matching.
#'
#' @param df     A data.frame with columns `ensembl_transcript_id` and `utr3`.
#' @param motifs A character vector of RNA motifs (e.g., "AUAAAUA").
#' @return A data.frame with columns: Transcript_ID, motif, start, end.
find_motif_positions <- function(df, motifs) {

  df <- df |>
    mutate(utr3_rna = str_replace_all(utr3, "T", "U"))

  seqs <- RNAStringSet(df$utr3_rna)
  names(seqs) <- df$ensembl_transcript_id

  results <- list()

  for (motif in motifs) {
    hits      <- vmatchPattern(motif, seqs, fixed = TRUE)
    flat_hits <- unlist(hits)

    if (length(flat_hits) > 0) {
      hit_names <- rep(names(hits), times = elementNROWS(hits))

      motif_df <- data.frame(
        Transcript_ID = hit_names,
        motif         = motif,
        start         = start(flat_hits),
        end           = end(flat_hits),
        stringsAsFactors = FALSE
      )
      results[[motif]] <- motif_df
    }
  }

  if (length(results) > 0) {
    return(do.call(rbind, results))
  } else {
    return(data.frame(
      Transcript_ID = character(),
      motif         = character(),
      start         = integer(),
      end           = integer(),
      stringsAsFactors = FALSE
    ))
  }
}


#' Combine p-values using Fisher's method
#'
#' @param pvals Numeric vector of p-values.
#' @return Combined p-value (chi-squared test).
combine_p_fisher <- function(pvals) {
  pvals <- pmax(pvals, 1e-300)
  X2    <- -2 * sum(log(pvals))
  pchisq(X2, df = 2 * length(pvals), lower.tail = FALSE)
}


#' Meta-analysis of k-mer enrichment across control iterations
#'
#' Combines Fisher's exact test p-values from multiple control comparisons
#' using Fisher's method, and computes log2 fold change of normalized
#' k-mer frequencies (target vs. mean background).
#'
#' @param df A data.frame with columns: kmer, pvalue, iteration, target, background.
#' @return A data.frame with columns: kmer, target_freq, bg_mean, log2fc, fisher_padj, k.
meta_analysis <- function(df) {

  pval_table <- df |>
    dplyr::select(kmer, pvalue, iteration) |>
    pivot_wider(names_from = iteration, values_from = pvalue, names_prefix = "pval_")

  pval_table$fisher_p <- apply(
    pval_table |> dplyr::select(starts_with("pval_")), 1, combine_p_fisher
  )
  pval_table$fisher_padj <- p.adjust(pval_table$fisher_p, method = "BH")

  fc_table <- df |>
    group_by(iteration) |>
    mutate(
      target_freq = target / sum(target),
      bg_freq     = background / sum(background)
    ) |>
    ungroup() |>
    dplyr::transmute(kmer, target_freq, bg_freq, iteration) |>
    pivot_wider(names_from = iteration, values_from = bg_freq, names_prefix = "bg_") |>
    transmute(
      kmer, target_freq,
      bg_mean = rowMeans(across(starts_with("bg_")))
    ) |>
    mutate(log2fc = log2((target_freq + 1e-6) / (bg_mean + 1e-6)))

  df_fc_p <- left_join(
    fc_table,
    pval_table |> dplyr::select(kmer, fisher_padj),
    by = "kmer"
  ) |>
    mutate(k = nchar(kmer))

  return(df_fc_p)
}


#' Keep only the most distal (non-overlapping) motif hits per transcript
#'
#' Starting from the 3' end, greedily keeps motifs that don't overlap.
#'
#' @param tab A data.frame with columns: motif, end (and optionally start).
#' @return Filtered data.frame with non-overlapping motif hits.
keep_most_distal <- function(tab) {
  tab <- tab |>
    mutate(width = nchar(motif),
           start = end - width + 1L) |>
    arrange(desc(end), desc(width))

  kept      <- logical(nrow(tab))
  min_start <- Inf

  for (i in seq_len(nrow(tab))) {
    if (tab$end[i] < min_start) {
      kept[i]   <- TRUE
      min_start <- tab$start[i]
    }
  }

  tab |> filter(kept)
}


#' Create a FASTA file from a data.frame of UTR sequences
#'
#' @param df       A data.frame with `ensembl_transcript_id` and `utr3` columns.
#' @param filename Output FASTA file path.
#' @param as_rna   If TRUE, convert T → U and write as RNA sequences.
create_fasta <- function(df, filename, as_rna = FALSE) {
  if (as_rna) {
    seqs_rna <- str_replace_all(df$utr3, "T", "U")
    seqs     <- RNAStringSet(seqs_rna)
  } else {
    seqs <- DNAStringSet(df$utr3)
  }
  names(seqs) <- df$ensembl_transcript_id
  writeXStringSet(seqs, filepath = filename)
}
