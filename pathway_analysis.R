#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(ggplot2)
})

result_file <- file.path("results", "differential_expression.csv")
out_dir <- file.path("results", "pathways")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(result_file)) stop("Run analysis.R before pathway_analysis.R.")

results <- read.csv(result_file, stringsAsFactors = FALSE)
if (!all(c("gene_id", "pvalue", "padj") %in% names(results))) {
  stop("Differential-expression table is missing required columns.")
}
query <- unique(results$gene_id[!is.na(results$padj) & results$padj < 0.05])
background <- unique(results$gene_id[!is.na(results$pvalue)])
if (length(query) == 0L) stop("No genes pass the adjusted p-value threshold.")

# Reuse the recorded result by default so the report is reproducible without a
# network request. Set REFRESH_GPROFILER=1 to query the current g:Profiler release.
cache_path <- file.path(out_dir, "gprofiler_response.json")
refresh <- identical(Sys.getenv("REFRESH_GPROFILER", "0"), "1")
if (file.exists(cache_path) && !refresh) {
  body <- paste(readLines(cache_path, warn = FALSE), collapse = "\n")
} else {
  payload <- list(
    organism = "hsapiens",
    query = unname(query),
    sources = c("GO:BP", "REAC"),
    user_threshold = 0.05,
    significance_threshold_method = "fdr",
    domain_scope = "custom_annotated",
    background = unname(background),
    output = "json",
    no_evidences = FALSE
  )
  request <- request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/") |>
    req_method("POST") |>
    req_user_agent("OsteoarthritisPainRNAseq/1.0") |>
    req_body_json(payload)
  response <- req_perform(request)
  body <- resp_body_string(response)
  writeLines(body, cache_path)
}
parsed <- fromJSON(body, flatten = TRUE)

if (is.null(parsed$result) || length(parsed$result) == 0L) {
  enrichment <- data.frame()
} else {
  enrichment <- as.data.frame(parsed$result, stringsAsFactors = FALSE)
  enrichment <- enrichment[, c("source", "native", "name", "p_value", "term_size",
                               "effective_domain_size", "query_size", "intersection_size",
                               "precision", "recall"), drop = FALSE]
  enrichment <- enrichment[order(enrichment$p_value), , drop = FALSE]
}
write.csv(enrichment, file.path(out_dir, "go_reactome_enrichment.csv"), row.names = FALSE)

run_date <- if (!is.null(parsed$timestamp)) substr(parsed$timestamp, 1, 10) else as.character(Sys.Date())
run_version <- if (!is.null(parsed$meta$version)) parsed$meta$version else "not returned"
writeLines(c(
  paste("g:Profiler date:", run_date),
  paste("g:Profiler version:", run_version),
  paste("DE genes submitted:", length(query)),
  paste("Tested genes used as background:", length(background)),
  paste("Significant GO:BP / Reactome terms (FDR < 0.05):", nrow(enrichment))
), file.path(out_dir, "enrichment_summary.txt"))

if (nrow(enrichment) > 0L) {
  top <- head(enrichment, 15)
  top$name <- factor(top$name, levels = rev(top$name))
  p <- ggplot(top, aes(x = -log10(p_value), y = name, color = source)) +
    geom_point(aes(size = intersection_size)) +
    labs(x = "-log10 FDR-adjusted p-value", y = NULL, size = "Genes in overlap", color = "Database") +
    theme_classic()
  ggsave(file.path(out_dir, "top_pathways.png"), p, width = 9, height = 6, dpi = 300)
}
cat(readLines(file.path(out_dir, "enrichment_summary.txt")), sep = "\n")
if (nrow(enrichment) > 0L) print(head(enrichment, 15))
