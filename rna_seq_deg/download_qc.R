# dependencies
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("recount3", "DESeq2", "edgeR"))
# install.packages(c("tidyverse", "data.table", "here"))

suppressPackageStartupMessages({
    library(recount3)
    library(DESeq2)
    library(edgeR)
    library(tidyverse)
    library(data.table)
    library(here)
})

# paths
output_dir <- here::here("outputs")
dir.create(output_dir, showWarnings = FALSE)

# download TCGA-BRCA via recount3
# recount3 provides raw counts directly, no manual file download needed
# project "BRCA" under "data_sources/tcga" gives ~1,100 breast tumor samples
message("Downloading TCGA-BRCA raw counts via recount3 (may take a few minutes)...")
rse <- create_rse_manual(
    project      = "BRCA",
    project_home = "data_sources/tcga",
    organism     = "human",
    annotation   = "gencode_v26",
    type         = "gene"
)
message(sprintf("RSE downloaded: %d genes x %d samples", nrow(rse), ncol(rse)))

# extract raw counts
# recount3 stores counts in the "raw_counts" assay as a matrix of integers
counts_raw <- assay(rse, "raw_counts")

# extract and clean metadata
# TCGA clinical variables are stored in colData(rse)
message("Cleaning metadata...")
meta_raw <- as.data.frame(colData(rse))

meta <- meta_raw %>%
    rownames_to_column("sample_id") %>%
    mutate(
        er_status         = tcga.xml_breast_carcinoma_estrogen_receptor_status,
        pr_status         = tcga.xml_breast_carcinoma_progesterone_receptor_status,
        her2_status       = tcga.xml_lab_proc_her2_neu_immunohistochemistry_receptor_status,
        age               = as.integer(tcga.xml_age_at_initial_pathologic_diagnosis),
        pathologic_stage  = tcga.cgc_case_pathologic_stage,
        histological_type = tcga.xml_histological_type
    ) %>%
    select(sample_id, er_status, pr_status, her2_status, age, pathologic_stage, histological_type) %>%
    filter(er_status %in% c("Positive", "Negative")) %>%
    filter(her2_status %in% c("Positive", "Negative")) %>%
    filter(!is.na(pathologic_stage)) %>%
    mutate(
        er_status         = factor(er_status,   levels = c("Negative", "Positive")),
        her2_status       = factor(her2_status, levels = c("Negative", "Positive")),
        pr_status         = factor(pr_status,   levels = c("Negative", "Positive")),
        age               = as.integer(age),
        pathologic_stage  = factor(pathologic_stage),
        histological_type = factor(histological_type)
    )

message(sprintf("Samples after metadata cleaning: %d", nrow(meta)))
message("ER status distribution:")
print(table(meta$er_status))
message("HER2 status distribution:")
print(table(meta$her2_status))
message("Pathologic stage distribution:")
print(table(meta$pathologic_stage))
message("Histological type distribution:")
print(table(meta$histological_type))

# align samples between counts and metadata
common_samples <- intersect(meta$sample_id, colnames(counts_raw))
message(sprintf("Samples in common (counts + metadata): %d", length(common_samples)))

counts_aligned <- counts_raw[, common_samples]
meta_aligned <- meta %>% filter(sample_id %in% common_samples) %>%
                arrange(match(sample_id, common_samples))

# sanity check
stopifnot(all(colnames(counts_aligned) == meta_aligned$sample_id))

# QC filtering
# remove genes with zero counts across all samples
nonzero_genes <- rowSums(counts_aligned > 0) > 0
counts_filtered <- counts_aligned[nonzero_genes, ]
message(sprintf("Genes after removing all-zero genes: %d", nrow(counts_filtered)))

# require expression in at least 10% of samples (edgeR filterByExpr approach)
library(edgeR)
dge_tmp <- DGEList(counts = counts_filtered)
keep <- filterByExpr(dge_tmp, group = meta_aligned$er_status,
                     min.count = 10, min.total.count = 15)
counts_filtered <- counts_filtered[keep, ]
message(sprintf("Genes after filterByExpr: %d", nrow(counts_filtered)))

# sample-level QC: flag low-library-size outliers (<3 SD below mean log lib size)
lib_sizes <- colSums(counts_filtered)
log_lib   <- log2(lib_sizes)
low_lib_flag <- log_lib < (mean(log_lib) - 3 * sd(log_lib))
if (any(low_lib_flag)) {
    message(sprintf("Removing %d low-library-size samples", sum(low_lib_flag)))
    counts_filtered <- counts_filtered[, !low_lib_flag]
    meta_aligned    <- meta_aligned[!low_lib_flag, ]
}

message(sprintf("Final dimensions: %d genes x %d samples",
                nrow(counts_filtered), ncol(counts_filtered)))

# save outputs
message("Saving cleaned files...")

fwrite(as.data.frame(counts_filtered), file.path(output_dir, "counts_filtered.csv"),
       row.names = TRUE)

write_csv(meta_aligned, file.path(output_dir, "metadata_clean.csv"))

message("QC complete. Files written to outputs/")
message("\nNext: run deseq2_analysis.R")