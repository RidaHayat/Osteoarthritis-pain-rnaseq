#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

counts_path <- file.path("data", Sys.getenv("COUNTS_FILE", "counts_geo.tsv"))
metadata_path <- file.path("data", "metadata.csv")
out_dir <- "results"
fig_dir <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path) || !file.exists(metadata_path)) {
  stop("Missing data/counts_geo.tsv or data/metadata.csv; see data/README.md.")
}

count_table <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(count_table) < 2L) stop("Count table needs gene IDs and sample count columns.")
raw_gene_ids <- count_table[[1]]
gene_ids <- sub("_.*$", "", raw_gene_ids)
gene_symbols <- ifelse(grepl("_", raw_gene_ids), sub("^[^_]+_", "", raw_gene_ids), NA_character_)
count_matrix <- as.matrix(count_table[-1])
rownames(count_matrix) <- gene_ids
suppressWarnings(storage.mode(count_matrix) <- "numeric")
if (anyNA(count_matrix) || any(count_matrix < 0)) {
  stop("Expected counts must be non-negative and have no missing values.")
}
if (anyDuplicated(rownames(count_matrix))) stop("Gene IDs are not unique after removing symbol suffixes.")
# The supplied gene-level expected-count matrix contains fractional values. DESeq2
# requires integer counts; round to the nearest integer and document this choice.
expected_count_matrix <- count_matrix
count_matrix <- round(expected_count_matrix)
storage.mode(count_matrix) <- "integer"

metadata <- read.csv(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Sample.ID", "Phenotype")
if (!all(required %in% names(metadata))) stop("Metadata needs Sample.ID and Phenotype columns.")
if (anyDuplicated(metadata$Sample.ID) || anyDuplicated(colnames(count_matrix))) stop("Sample IDs must be unique.")
if (!setequal(metadata$Sample.ID, colnames(count_matrix))) stop("Count and metadata sample IDs do not match.")
metadata <- metadata[match(colnames(count_matrix), metadata$Sample.ID), , drop = FALSE]
metadata$Phenotype <- factor(metadata$Phenotype, levels = c("Pain:High", "Pain:Low"))
if (anyNA(metadata$Phenotype)) stop("Phenotype must be Pain:High or Pain:Low.")
group_sizes <- table(metadata$Phenotype)
if (!identical(as.integer(group_sizes), c(5L, 5L))) stop("Expected five samples per pain group; inspect metadata.")

# Keep only genes with more than one count across all samples.
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ Phenotype
)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)

# Positive fold change means higher expression in low-pain samples.
res <- results(dds, contrast = c("Phenotype", "Pain:Low", "Pain:High"), alpha = 0.05)
res_table <- as.data.frame(res)
res_table$gene_id <- rownames(res_table)
res_table$gene_symbol <- gene_symbols[match(res_table$gene_id, gene_ids)]
res_table <- res_table[order(res_table$padj, na.last = TRUE), ]
write.csv(res_table, file.path(out_dir, "differential_expression.csv"), row.names = FALSE)

# Compare the nearest-integer choice with the truncation used in the original Rmd.
truncated_matrix <- trunc(expected_count_matrix)
storage.mode(truncated_matrix) <- "integer"
dds_truncated <- DESeqDataSetFromMatrix(
  countData = truncated_matrix,
  colData = metadata,
  design = ~ Phenotype
)
dds_truncated <- dds_truncated[rowSums(counts(dds_truncated)) > 1, ]
dds_truncated <- DESeq(dds_truncated, quiet = TRUE)
res_truncated <- results(
  dds_truncated,
  contrast = c("Phenotype", "Pain:Low", "Pain:High"),
  alpha = 0.05
)
rounding_sensitivity <- data.frame(
  conversion = c("round to nearest integer", "truncate toward zero"),
  genes_retained = c(nrow(res_table), nrow(res_truncated)),
  genes_with_pvalue = c(sum(!is.na(res_table$pvalue)), sum(!is.na(res_truncated$pvalue))),
  genes_BH_adjusted_p_below_0_05 = c(sum(res_table$padj < 0.05, na.rm = TRUE), sum(res_truncated$padj < 0.05, na.rm = TRUE))
)
write.csv(rounding_sensitivity, file.path(out_dir, "rounding_sensitivity.csv"), row.names = FALSE)

sig <- subset(res_table, !is.na(padj) & padj < 0.05)
n_up_low <- sum(sig$log2FoldChange > 0)
n_up_high <- sum(sig$log2FoldChange < 0)
writeLines(c(
  "Comparison: Pain:Low vs Pain:High",
  paste("Samples: Pain:High =", unname(group_sizes[["Pain:High"]]), "; Pain:Low =", unname(group_sizes[["Pain:Low"]])),
  paste("Genes retained after low-count filtering:", nrow(res_table)),
  paste("Genes with a non-missing test p-value:", sum(!is.na(res_table$pvalue))),
  paste("Genes with BH-adjusted p-value < 0.05:", nrow(sig)),
  paste("Higher in Pain:Low among significant genes:", n_up_low),
  paste("Higher in Pain:High among significant genes:", n_up_high),
  "Input: RSEM expected counts rounded to nearest integer for DESeq2.",
  paste("Truncation sensitivity: BH-adjusted p < 0.05 for", sum(res_truncated$padj < 0.05, na.rm = TRUE), "genes.")
), file.path(out_dir, "analysis_summary.txt"))

# PCA on variance-stabilized expression values.
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "Phenotype", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
p_pca <- ggplot(pca_data, aes(PC1, PC2, color = Phenotype)) +
  geom_point(size = 3) +
  labs(
    title = "Knee OA synovium: pain groups",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Pain group"
  ) + theme_classic()
ggsave(file.path(fig_dir, "pca_by_pain_group.png"), p_pca, width = 7, height = 5, dpi = 300)

# DESeq2 MA plot.
png(file.path(fig_dir, "ma_plot.png"), width = 2100, height = 1500, res = 300)
plotMA(res, ylim = c(-5, 5), main = "Low pain vs high pain")
dev.off()

# Volcano plot. DESeq2's adjusted p-values use the BH procedure by default.
volcano_data <- subset(res_table, !is.na(padj) & is.finite(log2FoldChange))
volcano_data$neg_log10_padj <- -log10(pmax(volcano_data$padj, .Machine$double.xmin))
volcano_data$category <- ifelse(!is.na(volcano_data$padj) & volcano_data$padj < 0.05, "Adjusted p < 0.05", "Not significant")
p_volcano <- ggplot(volcano_data, aes(log2FoldChange, neg_log10_padj, color = category)) +
  geom_point(alpha = 0.65, size = 1) +
  scale_color_manual(values = c("Adjusted p < 0.05" = "#C44E52", "Not significant" = "#8C8C8C")) +
  labs(
    title = "Differential expression: low pain vs high pain",
    x = "Log2 fold change (low pain / high pain)",
    y = "-log10 adjusted p-value",
    color = NULL
  ) + theme_classic()
ggsave(file.path(fig_dir, "volcano.png"), p_volcano, width = 7, height = 5, dpi = 300)

cat(readLines(file.path(out_dir, "analysis_summary.txt")), sep = "\n")
cat("\nTop genes by adjusted p-value:\n")
print(head(res_table[, c("gene_id", "gene_symbol", "log2FoldChange", "pvalue", "padj")], 10))
