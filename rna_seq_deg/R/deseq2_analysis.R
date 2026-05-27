# run DESeq2 for two biologically meaningful contrasts
# dependencies
# BiocManager::install(c("DESeq2", "BiocParallel"))
# install.packages(c("tidyverse", "data.table", "here"))

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(tidyverse)
  library(data.table)
  library(here)
})

# SnowParam for windows, MultiCoreParam for linux
register(SnowParam(4))

output_dir <- here::here("outputs")

# load data
message("Loading filtered counts and metadata...")
counts <- fread(file.path(output_dir, "counts_filtered.csv"), data.table = FALSE)
rownames(counts) <- counts[[1]]
counts <- counts[,-1]
counts <- as.matrix(counts)
mode(counts) <- "integer"

meta <- read_csv(file.path(output_dir, "metadata_clean.csv"),
                 show_col_types = FALSE) %>%
                 column_to_rownames("sample_id")

counts <- counts[, rownames(meta)]
stopifnot(all(colnames(counts) == rownames(meta)))

# helper: run DESeq2 for a given contrast
run_deseq2 <- function(counts, meta, group_var, ref_level, contrast_level,
                       label, covariates = c("pathologic_stage")) {
    message(sprintf("\nRunning DESeq2: %s ", label))

    keep <- !is.na(meta[[group_var]])
    counts_sub <- counts[, keep]
    meta_sub <- meta[keep, , drop = FALSE]

    # drop unused factor levels in pathologic_stage to avoid deseq2 errors
    meta_sub$pathologic_stage <- droplevels(meta_sub$pathologic_stage)

    formula_str <- paste("~", paste(c(covariates, group_var), collapse = " + "))
    message(sprintf("  Formula: %s", formula_str))
    message(sprintf("  Samples: %d (%s=%s: %d | %s=%s: %d)",
        ncol(counts_sub),
        group_var, ref_level, sum(meta_sub[[group_var]] == ref_level, na.rm = TRUE),
        group_var, contrast_level, sum(meta_sub[[group_var]] == contrast_level, na.rm = TRUE)
    ))

    meta_sub[[group_var]] <- relevel(factor(meta_sub[[group_var]]), ref = ref_level)
    dds <- DESeqDataSetFromMatrix(countData = counts_sub,
                                  colData = meta_sub,
                                  design = as.formula(formula_str))
    dds <- DESeq(dds, parallel = TRUE)

    # apelgm requires the coef name, not a contrast vector
    # coef name format: groupvar_level_vs_reflevel

    coef_name <- paste0(group_var, "_", contrast_level, "_vs_", ref_level)
    message(sprintf("  Shrinkage coef: %s", coef_name))

    res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm", parallel = TRUE)
    res_df <- as.data.frame(res_shrunk) %>%
            rownames_to_column("gene_id") %>%
            mutate(
                significance = case_when(
                    is.na(padj) ~ "Filtered",
                    padj < 0.05 & log2FoldChange > 1 ~ "Up",
                    padj < 0.05 & log2FoldChange < -1 ~ "Down",
                    padj < 0.05 ~ "Sig (|LFC|<1)",
                    TRUE ~ "NS"
                ),
                abs_lfc = abs(log2FoldChange),
                neg_log10_padj = -log10(padj)
            ) %>%
            arrange(padj, desc(abs_lfc))
    
    message(sprintf("  Significant DEGs (padj < 0.05 & |LFC| > 1): Up=%d, Down=%d",
        sum(res_df$significance == "Up", na.rm = TRUE),
        sum(res_df$significance == "Down", na.rm = TRUE)
    ))

    list(dds = dds, results = res_df)
}

# contrast A: ER+ vs ER-
er_out <- run_deseq2(
    counts, meta,
    group_var = "er_status",
    ref_level = "Negative",
    contrast_level = "Positive",
    label = "ER+ vs ER-"
)