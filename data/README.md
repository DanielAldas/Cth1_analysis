# Data Directory

This directory contains the input data files required for reproducing the figures on our manuscript.

**These files are not tracked by git** — add `data/` to `.gitignore` if data files
contain sensitive or large datasets. Instead, document provenance here.

## Files

| File | Description | Source | Date |
|------|-------------|--------|------|
| `original_KD_DGE_results.csv` | DESeq2 DGE results from Cas13d cth1 knockdown experiments (gRNA2 and gRNA7 at 2 and 4 hpf) | Internal pipeline | [2025-07-21] |
| `DGE_new_OE.csv` | DESeq2 DGE results from overexpression experiments, including alpha-amanitin–treated conditions | Internal pipeline | [2025-07-21] |
| `utr3_protein_coding_longest_tx.txt` | 3'UTR lengths for protein-coding genes (longest transcript per gene); columns: Gene_ID, Transcript_ID, Len, seq_3utr | Ariel (collaborator) | [2025-05-29] |
| `santi_timecourse.csv` | Wild-type zebrafish timecourse TPM values (0–8 hpf) across RiboZero, PolyA, and alpha-amanitin–treated libraries | Medina-Muñoz et al. (2021) *Genome Biology* 22:14, Additional file 2: Table S1. DOI: [10.1186/s13059-020-02251-5](https://doi.org/10.1186/s13059-020-02251-5). Stowers Original Data Repository Accession: [LIBPB-1584](https://www.stowers.org/research/publications/libpb-1584) | 2021-01-05 |
| `Zygotic.csv` | Pure zygotic gene IDs (column: gene_id) | Baia-Amaral et al. (2024) *Genome Biology*  25:74, Additional file 2 DOI: [10.1186/s13059-024-03197-8](https://doi.org/10.1186/s13059-024-03197-8)| [2024-03-19] |
|`combined_tpm.csv`| RSEM TPM table from Cas13d cth1 knockdown experiments (gRNA2 and gRNA7 at 4 hpf) of Poly(A) enrichment and RiboZero RNA-seq libraries with samples treated (or not) with alpha-amanitin | Internal pipeline | [2026-03-19] |

## Notes

- DGE results use Ensembl gene IDs (ENSDARG format) as primary keys.
- 3'UTR length filtering: 100–800 nt; expression filtering: CPM ≥ 5 in cas13d-alone controls.
- Ensembl version 106 was used for all genome annotations and biomaRt queries.


