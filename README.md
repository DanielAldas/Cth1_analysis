# Cth1 is a spatio-temporally regulated maternal RNA decay factor essential for zebrafish development

This repository contains the R analysis code and data used to reproduce the figures in:

> **Cth1 is a spatio-temporally regulated maternal RNA decay factor essential for zebrafish development**  
> Gopal Kushawah, V. Daniel Aldas Bulos, Danielson Baia Amaral, Huzaifa Hassan, Stephanie H. Nowotarski, Ariel A. Bazzini.  
> DOI: [pending]

**Summary:** Maternal mRNA decay is essential for the maternal-to-zygotic transition (MZT), yet its spatial regulation has remained largely unexplored. We identify Cth1 as the first maternal RNA decay factor in zebrafish, acting independently of zygotic genome activation and regulated in a spatiotemporal manner by conserved 3′UTR cis-elements. Cas13d-mediated depletion of maternal *cth1* revealed early developmental defects that could not be assessed with Cas9 mutants due to infertility from impaired gametogenesis. Mechanistically, Cth1 recognizes AU-rich elements in maternal 3′UTRs to promote deadenylation and decay; deletion of these motifs or Cth1 knockdown stabilizes both reporter and endogenous transcripts. These findings establish Cth1 as a key component of the maternal program and provide the first direct evidence that spatially regulated mRNA decay, outside the germline, contributes to early vertebrate development.

---

## Repository structure

```
├── src/                          # Analysis scripts (run in numbered order)
│   ├── 00_setup.R                # Shared configuration, helpers, color palettes
│   ├── 01_define_groups.Rmd      # Define gene groups (g3, c1–c4)
│   ├── 02_kd_scatter_plots.Rmd   # Knockdown concordance and overlap plots
│   ├── 03_timecourse_expression.Rmd
│   ├── 04_decay_analysis.Rmd
│   ├── 05_characterize_groups.Rmd
│   ├── 06_motif_enrichment.Rmd
│   ├── 07_endogenous_motif.Rmd
│   ├── 08_polya_status.Rmd
│   └── motif_enrichment.R        # K-mer counting helper (sourced by 06)
├── data/                         # Input data files (see data/README.md)
├── results/
│   ├── figures/                  # Generated figures (PDF/PNG)
│   └── kmers/                    # K-mer enrichment tables
└── doc/
```

---

## Dependencies

All scripts are written in R. The following packages are required.

**CRAN:**

```r
install.packages(c(
  "tidyverse", "data.table", "ggplot2", "ggpubr", "ggforce",
  "ggrepel", "ggseqlogo", "ggVennDiagram", "eulerr",
  "paletteer", "UpSetR", "stringdist"
))
```

**Bioconductor:**

```r
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("Biostrings", "biomaRt", "DESeq2", "motifmatchr", "msa"))
```

Genome annotation uses **Ensembl release 106** (*Danio rerio*, GRCz11). The `biomaRt` queries in `05_characterize_groups.Rmd` connect to the Ensembl archive at `https://apr2022.archive.ensembl.org`.

---

## Data

Input files live in `data/`. See [`data/README.md`](data/README.md) for full provenance and download details.

| File | Description |
|------|-------------|
| `original_KD_DGE_results.csv` | DESeq2 results from Cas13d *cth1* knockdown (gRNA2 and gRNA7 at 2 and 4 hpf) |
| `DGE_new_OE.csv` | DESeq2 results from *cth1* overexpression and alpha-amanitin conditions |
| `combined_tpm.csv` | RSEM TPM values from knockdown experiments (Poly(A) and RiboZero libraries) |
| `santi_timecourse.csv` | Wild-type zebrafish timecourse TPM (0–8 hpf); from [Medina-Muñoz et al. 2021](https://doi.org/10.1186/s13059-020-02251-5) |
| `Zygotic.csv` | Pure zygotic gene IDs; from [Baia-Amaral et al. 2024](https://doi.org/10.1186/s13059-024-03197-8) |
| `utr3_protein_coding_longest_tx_Ariel.txt` | 3′UTR sequences for protein-coding genes (longest transcript per gene) |

Raw RNA-seq FASTQ files will be deposited in GEO (accession pending).

**Filtering criteria applied throughout:** 3′UTR length 100–800 nt; expression ≥ 5 CPM in Cas13d-alone controls.

---

## Running the analysis

Scripts are numbered and must be run **in order** from the `src/` directory. Each script saves an `.RData` checkpoint to `results/` that downstream scripts depend on.

To knit all RMarkdown files sequentially from R:

```r
library(rmarkdown)
scripts <- list.files("src", pattern = "^0[1-8].*\\.Rmd$", full.names = TRUE)
for (f in scripts) {
  render(f, knit_root_dir = normalizePath(".."))
}
```

Or knit a single script, setting the working directory to the project root:

```r
rmarkdown::render("src/03_timecourse_expression.Rmd", knit_root_dir = "..")
```

### Script summary

| Script | Description | Key output figures |
|--------|-------------|-------------------|
| `01_define_groups.Rmd` | Identifies *cth1*-responsive genes (g3: logFC ≥ 0.5, padj ≤ 0.05 for both gRNAs) and draws four random control sets (c1–c4, n = 500 each) | — |
| `02_kd_scatter_plots.Rmd` | gRNA2 vs. gRNA7 concordance; Venn/Euler diagrams of up- and down-regulated overlap | `scatter_KD_groups.pdf`, `venn_upregulated.png`, `venn_downregulated.png` |
| `03_timecourse_expression.Rmd` | Expression dynamics of g3 vs. controls across wild-type 0–8 hpf timecourse (RiboZero, Poly(A), alpha-amanitin) | `ribozero_timecourse_sinaplot.pdf`, `polya_timecourse_sinaplot.pdf`, `aamanitin_sina.pdf` |
| `04_decay_analysis.Rmd` | RNA decay rates (log2 FC relative to 2 hpf); Wilcoxon and ECDF comparisons | `decay_ribo_sina.pdf`, `decay_polya_sina.pdf`, `decay_alpha_sina.pdf` |
| `05_characterize_groups.Rmd` | CDS, 3′UTR, and 5′UTR length distributions queried from Ensembl v106 | `cds_sinaplot.pdf`, `3utr_sinaplot.pdf`, `5utr_sinaplot.pdf` |
| `06_motif_enrichment.Rmd` | K-mer enrichment (5–10 nt) in g3 3′UTRs vs. each control; Fisher's method meta-analysis; volcano plots and sequence logos | `volcano_{5–10}mers_g3.pdf`, `distance_plot.pdf`, `relative_position_plot.pdf` |
| `07_endogenous_motif.Rmd` | Classifies all filtered genes by enriched motif presence; links motif occupancy to KD response and RNA decay | `stack_endogenous_with_motif.pdf` |
| `08_polya_status.Rmd` | Poly(A) tail dynamics (log2 Poly(A)/RiboZero TPM); KS and Wilcoxon tests between KD and control | `endogenous_polystatus.pdf`, `polyA_FC_*.pdf` |

---

## Citation

If you use this code, please cite:

> Kushawah G, Aldas Bulos VD, Baia Amaral D, Hassan H, Nowotarski SH, Bazzini AA. Cth1 is a spatio-temporally regulated maternal RNA decay factor essential for zebrafish development. DOI: [pending]
