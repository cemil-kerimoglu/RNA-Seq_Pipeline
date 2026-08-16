#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(topGO)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9) {
  stop("Usage: topgo.R de_all.tsv selected.tsv counts.tsv out.tsv ontology universe_mode node_size paper_weight_p label")
}
de_file <- args[1]
selected_file <- args[2]
counts_file <- args[3]
out_file <- args[4]
ontology <- args[5]
universe_mode <- args[6]
node_size <- as.integer(args[7])
paper_weight_p <- as.numeric(args[8])
label <- args[9]

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

de <- read.delim(de_file, stringsAsFactors = FALSE, check.names = FALSE)
sel <- read.delim(selected_file, stringsAsFactors = FALSE, check.names = FALSE)
counts_df <- read.delim(counts_file, stringsAsFactors = FALSE, check.names = FALSE)

strip_version <- function(x) sub("\\.[0-9]+$", "", x)
de_ids <- unique(strip_version(de$gene_id))
sel_ids <- unique(strip_version(sel$gene_id))
annotated_count_ids <- unique(strip_version(counts_df$gene_id))

# org.Mm.eg.db provides the Ensembl-to-GO mapping used by topGO. Restricting to
# mapped identifiers prevents unannotated IDs from entering the GO graph.
ens_keys <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
if (universe_mode == "tested") {
  universe <- intersect(de_ids, ens_keys)
} else if (universe_mode == "all_annotated") {
  # This is the closest switch in this workflow to the paper's stated
  # "all genes in the genome" universe: all genes in the pinned count annotation
  # that can be mapped by the installed mouse annotation package.
  universe <- intersect(annotated_count_ids, ens_keys)
} else {
  stop("universe_mode must be 'tested' or 'all_annotated'")
}
selected <- intersect(sel_ids, universe)

if (length(universe) == 0) stop("No universe genes map to org.Mm.eg.db")

gene_list <- as.integer(universe %in% selected)
names(gene_list) <- universe

GOdata <- new(
  "topGOdata",
  description = label,
  ontology = ontology,
  allGenes = gene_list,
  geneSel = function(x) x == 1L,
  annot = annFUN.org,
  mapping = "org.Mm.eg.db",
  ID = "ensembl",
  nodeSize = node_size
)

classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
weight <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

n_nodes <- length(usedGO(GOdata))
tab <- GenTable(
  GOdata,
  classicFisher = classic,
  weightFisher = weight,
  weight01Fisher = weight01,
  orderBy = "weight01Fisher",
  ranksOf = "weight01Fisher",
  topNodes = n_nodes
)

parse_p <- function(x) {
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}
tab$classicFisher_numeric <- parse_p(tab$classicFisher)
tab$weightFisher_numeric <- parse_p(tab$weightFisher)
tab$weight01Fisher_numeric <- parse_p(tab$weight01Fisher)
# topGO's graph-aware tests are not independent; this BH column is supplied as
# a descriptive convenience, not as a replacement for interpreting topGO's
# graph-aware algorithms.
tab$weight01_BH_descriptive <- p.adjust(tab$weight01Fisher_numeric, method = "BH")
tab$paper_like_weight_pass <- !is.na(tab$weightFisher_numeric) & tab$weightFisher_numeric < paper_weight_p
tab$paper_like_weight01_pass <- !is.na(tab$weight01Fisher_numeric) & tab$weight01Fisher_numeric < paper_weight_p
tab$gene_set <- label
tab$universe_mode <- universe_mode
tab$universe_size <- length(universe)
tab$selected_size <- length(selected)

write.table(tab, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  c(
    paste("Gene set:", label),
    paste("Ontology:", ontology),
    paste("Universe mode:", universe_mode),
    paste("Universe genes with org.Mm.eg.db mapping:", length(universe)),
    paste("Selected genes with mapping:", length(selected)),
    paste("GO nodes tested:", n_nodes),
    paste0("Paper numeric cutoff flags: weight and weight01 Fisher P < ", paper_weight_p)
  ),
  sub("\\.tsv$", ".summary.txt", out_file)
)
