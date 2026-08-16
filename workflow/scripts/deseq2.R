#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 10) {
  stop("Usage: deseq2.R counts.tsv samples.tsv gene_map.tsv outdir reference test min_count_per_sample paper_padj paper_abs_fc modern_fdr")
}
counts_file <- args[1]
samples_file <- args[2]
gene_map_file <- args[3]
outdir <- args[4]
reference_level <- args[5]
test_level <- args[6]
min_count_per_sample <- as.integer(args[7])
paper_padj <- as.numeric(args[8])
paper_abs_fc <- as.numeric(args[9])
modern_fdr <- as.numeric(args[10])

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts_df <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(colnames(counts_df)[1] == "gene_id")
gene_ids <- counts_df$gene_id
count_mat <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(count_mat) <- gene_ids
storage.mode(count_mat) <- "integer"

meta <- read.delim(samples_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(colnames(count_mat) %in% meta$sample)) stop("Not all count columns occur in samples.tsv")
meta <- meta[match(colnames(count_mat), meta$sample), , drop = FALSE]
rownames(meta) <- meta$sample
if (any(is.na(meta$condition))) stop("Missing condition in sample metadata")
meta$condition <- relevel(factor(meta$condition), ref = reference_level)

# Current DESeq2-style prefilter: require a modest count in at least the size of
# the smallest experimental group. DESeq2 would also perform independent filtering
# later; this step mainly removes features with essentially no information.
smallest_group_size <- min(table(meta$condition))
keep <- rowSums(count_mat >= min_count_per_sample) >= smallest_group_size
message("Keeping ", sum(keep), " / ", length(keep), " genes with count >= ",
        min_count_per_sample, " in at least ", smallest_group_size, " samples")
dds <- DESeqDataSetFromMatrix(countData = count_mat[keep, , drop = FALSE], colData = meta, design = ~ condition)
dds <- DESeq(dds)

# Raw Wald-test result carries p-values/padj.  Shrunk LFC is added for stable ranking/visualization.
res_raw <- results(dds, contrast = c("condition", test_level, reference_level), alpha = modern_fdr)
coef_candidates <- resultsNames(dds)
coef_name <- grep(paste0("^condition_", test_level, "_vs_", reference_level, "$"), coef_candidates, value = TRUE)
if (length(coef_name) != 1) {
  # Defensive fallback for syntactically altered level names.
  coef_name <- grep("^condition_.*_vs_.*$", coef_candidates, value = TRUE)
}
if (length(coef_name) != 1) stop("Could not uniquely identify condition coefficient: ", paste(coef_candidates, collapse = ", "))
res_shr <- lfcShrink(dds, coef = coef_name, type = "apeglm")

res <- data.frame(
  gene_id = rownames(res_raw),
  baseMean = res_raw$baseMean,
  log2FoldChange_raw = res_raw$log2FoldChange,
  lfcSE = res_raw$lfcSE,
  stat = res_raw$stat,
  pvalue = res_raw$pvalue,
  padj = res_raw$padj,
  log2FoldChange_shrunken = res_shr$log2FoldChange,
  stringsAsFactors = FALSE
)
gene_map <- read.delim(gene_map_file, stringsAsFactors = FALSE, check.names = FALSE)
mi <- match(res$gene_id, gene_map$gene_id)
res$gene_name <- gene_map$gene_name[mi]
res$gene_type <- gene_map$gene_type[mi]
# Put human-readable annotation next to the stable identifier.
res <- res[, c("gene_id", "gene_name", "gene_type", setdiff(colnames(res), c("gene_id", "gene_name", "gene_type")))]
res$fold_change_raw <- 2^res$log2FoldChange_raw
res$direction <- ifelse(res$log2FoldChange_raw > 0, "up", ifelse(res$log2FoldChange_raw < 0, "down", "flat"))
res <- res[order(res$padj, res$pvalue, na.last = TRUE), ]
write.table(res, file.path(outdir, "deseq2_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

paper_lfc <- log2(paper_abs_fc)
paper_sig <- !is.na(res$padj) & res$padj < paper_padj & abs(res$log2FoldChange_raw) > paper_lfc
modern_sig <- !is.na(res$padj) & res$padj < modern_fdr

write.table(res[paper_sig, ], file.path(outdir, "deseq2_paper_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[paper_sig & res$log2FoldChange_raw > 0, ], file.path(outdir, "deseq2_paper_up.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[paper_sig & res$log2FoldChange_raw < 0, ], file.path(outdir, "deseq2_paper_down.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[modern_sig, ], file.path(outdir, "deseq2_fdr_significant.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[modern_sig & res$log2FoldChange_raw > 0, ], file.path(outdir, "deseq2_fdr_up.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[modern_sig & res$log2FoldChange_raw < 0, ], file.path(outdir, "deseq2_fdr_down.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

norm <- counts(dds, normalized = TRUE)
write.table(data.frame(gene_id = rownames(norm), norm, check.names = FALSE), file.path(outdir, "normalized_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
write.table(data.frame(gene_id = rownames(assay(vsd)), assay(vsd), check.names = FALSE), file.path(outdir, "vst_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# PCA
pca_df <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
p <- ggplot(pca_df, aes(PC1, PC2, label = name, shape = condition)) +
  geom_point(size = 3) + geom_text(vjust = -0.8, check_overlap = TRUE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size = 12)
ggsave(file.path(outdir, "PCA.pdf"), p, width = 7, height = 5)

# Sample-distance heatmap
sample_dists <- dist(t(assay(vsd)))
pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 7, height = 6)
pheatmap(as.matrix(sample_dists), annotation_col = data.frame(condition = meta$condition, row.names = rownames(meta)), border_color = NA)
dev.off()

# MA plot using shrunken effect sizes
pdf(file.path(outdir, "MA_apeglm.pdf"), width = 7, height = 5)
plotMA(res_shr, alpha = modern_fdr, ylim = c(-4, 4))
dev.off()

# Volcano plot: shrunken effect size, raw DESeq2 adjusted P value.
volc <- res
volc$minuslog10padj <- -log10(pmax(volc$padj, .Machine$double.xmin))
volc$significant <- ifelse(modern_sig, paste0("FDR<", modern_fdr), "not significant")
p <- ggplot(volc, aes(log2FoldChange_shrunken, minuslog10padj, shape = significant)) +
  geom_point(alpha = 0.55, size = 1.4, na.rm = TRUE) +
  geom_vline(xintercept = c(-paper_lfc, paper_lfc), linetype = "dashed") +
  geom_hline(yintercept = -log10(modern_fdr), linetype = "dotted") +
  xlab("log2 fold change (apeglm-shrunken)") + ylab("-log10 adjusted P value") +
  theme_bw(base_size = 12)
ggsave(file.path(outdir, "volcano.pdf"), p, width = 7, height = 5)

saveRDS(dds, file.path(outdir, "dds.rds"))
saveRDS(vsd, file.path(outdir, "vsd.rds"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

summary_lines <- c(
  paste("Samples:", ncol(dds)),
  paste("Genes tested:", nrow(dds)),
  paste0("Paper-like threshold: padj < ", paper_padj, " and |FC| > ", paper_abs_fc, " -> ", sum(paper_sig), " genes"),
  paste0("Modern reporting threshold: FDR < ", modern_fdr, " -> ", sum(modern_sig), " genes"),
  paste("Coefficient shrunk with apeglm:", coef_name)
)
writeLines(summary_lines, file.path(outdir, "summary.txt"))
