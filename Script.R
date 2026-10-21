############################################################
# Hub gene pipeline (REVISED): GSE161969 (decidua) + GSE306548 (placenta)
# DEGs (FDR) -> same-direction overlap + enrichment test
#   -> GO/KEGG -> PPI -> top 10 hubs
#
# Changes from the original "Both dataset.R":
#  1. Decidua DEGs are called on padj (FDR), not raw p.
#  2. Decidua log2FC is shrunk (apeglm) and low-count genes are dropped,
#     so extreme fold changes from near-zero counts do not enter the overlap.
#  3. Placenta can use YOUR OWN full DESeq2 table (placenta_file) instead of
#     the authors' precomputed one.
#  4. Overlap = same-direction genes only; discordant genes are saved
#     separately and never go into enrichment / PPI / hubs.
#  5. Overlap is tested (hypergeometric) against the genes tested in BOTH
#     datasets, overall and per direction, plus a threshold-free correlation.
#  6. GO / KEGG use that same tested-gene background as the universe.
#  7. Duplicate / debugging code removed.
############################################################

## ---- 0. Packages ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc <- c("GEOquery","limma","DESeq2","edgeR","clusterProfiler","org.Hs.eg.db",
          "STRINGdb","enrichplot","apeglm")
cran <- c("tidyverse","pheatmap","ggrepel","igraph","ggVennDiagram","ggraph","ashr")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), ask = FALSE)
install.packages(setdiff(cran, rownames(installed.packages())))

library(GEOquery); library(limma); library(DESeq2); library(edgeR)
library(tidyverse); library(pheatmap); library(ggrepel)
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(STRINGdb); library(igraph); library(ggVennDiagram)
library(ggraph)

############################################################
## ---- 1. SETTINGS ----
############################################################
options(timeout = 1200)
dir.create("data",    showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

species <- 9606                         # human

# Primary DEG thresholds: FDR < 0.05 and |log2FC| >= 1 in BOTH datasets.
deg_cfg <- list(
  GSE161969 = list(stat = "padj", p = 0.05, fc = 1),
  GSE306548 = list(stat = "padj", p = 0.05, fc = 1)
)

# Decidua (4 vs 3 samples) quality controls
use_shrunk_lfc <- TRUE   # TRUE: use apeglm-shrunk log2FC for the |log2FC| cutoff
basemean_min   <- 20     # drop decidua genes with DESeq2 baseMean below this (0 = off)

# Sensitivity analysis only (printed + saved, never used downstream):
exploratory_p <- 0.05    # raw p cutoff for decidua

# PLACENTA input.
#  NULL  -> download the authors' precomputed DESeq2 table from GEO
#  path  -> YOUR full DESeq2 results (ALL tested genes, first column = gene symbol), e.g.
#           "C:/RPL project/Placenta_GSE306548/Result306548.csv"
placenta_file <- NULL

# OPTIONAL manual overrides (leave NULL to auto-detect groups from GEO metadata)
manual_case <- NULL     # e.g. c("GSM4929001", ...)  (RPL)
manual_ctrl <- NULL     # e.g. c("GSM4929005", ...)  (control)

# OPTIONAL: force a specific counts file for GSE161969
counts_file_override <- NULL

############################################################
## ---- 2. HELPER FUNCTIONS ----
############################################################
read_any <- function(f) {
  if (grepl("\\.xlsx?$", f, ignore.case = TRUE)) {
    if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
    as.data.frame(readxl::read_excel(f))
  } else {
    sep <- if (grepl("\\.csv(\\.gz)?$", f, ignore.case = TRUE)) "," else "\t"
    read.delim(f, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  }
}

# Convert Ensembl IDs to symbols if needed
to_symbols <- function(ids) {
  ids <- as.character(ids)
  if (mean(grepl("^ENSG", ids)) > 0.8) {
    ids <- sub("\\..*", "", ids)
    sym <- mapIds(org.Hs.eg.db, ids, "SYMBOL", "ENSEMBL", multiVals = "first")
    return(unname(sym))
  }
  ids
}

pick_col <- function(nm, pattern) {
  hit <- grep(pattern, nm, ignore.case = TRUE, value = TRUE)
  if (length(hit)) hit[1] else NA_character_
}

# One-sided hypergeometric test for overlap of gene sets A and B within N tested genes
hyper_overlap <- function(a, b, N, label) {
  a <- unique(a); b <- unique(b)
  k <- length(intersect(a, b))
  e <- length(a) * length(b) / N
  data.frame(test = label, n_A = length(a), n_B = length(b), overlap = k,
             expected = round(e, 1), fold_enrichment = round(k / e, 2),
             p_value = phyper(k - 1, length(b), N - length(b), length(a), lower.tail = FALSE))
}

############################################################
## ---- 3. GSE161969 (RNA-seq, RPL vs control decidua) ----
############################################################
## 3a. Sample metadata
gse <- getGEO("GSE161969", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

txt_cols <- unique(c(intersect(c("title","source_name_ch1","description"), colnames(meta)),
                     grep("^characteristics", colnames(meta), value = TRUE)))
meta$txt <- apply(meta[, txt_cols, drop = FALSE], 1, paste, collapse = " | ")

case_pat <- "\\bRPL\\b|\\bRSA\\b|\\bRM\\b|\\bURSA\\b|recurrent|miscarriage|abortion|pregnancy loss"
ctrl_pat <- "control|normal|healthy|fertile|\\bNC\\b"

if (!is.null(manual_case) && !is.null(manual_ctrl)) {
  meta$group <- NA_character_
  meta$group[meta$geo_accession %in% manual_case] <- "case"
  meta$group[meta$geo_accession %in% manual_ctrl] <- "ctrl"
} else {
  is_case <- grepl(case_pat, meta$txt, ignore.case = TRUE)
  is_ctrl <- grepl(ctrl_pat, meta$txt, ignore.case = TRUE)
  meta$group <- ifelse(is_case & !is_ctrl, "case", ifelse(is_ctrl & !is_case, "ctrl", NA))
}

cat("\n--- GSE161969 sample groups (CHECK THIS) ---\n")
print(meta[, c("geo_accession", "title", "group")])
print(table(meta$group, useNA = "ifany"))

if (anyNA(meta$group) || min(table(meta$group)) < 2) {
  cat("\nFull text used for matching:\n"); print(meta$txt)
  stop("Could not assign groups reliably. Set manual_case / manual_ctrl (GSM IDs) in SETTINGS, ",
       "then re-run from section 3.")
}

## 3b. Download supplementary file(s) with the count data
supp_dir <- "data/GSE161969"
if (length(list.files(supp_dir)) == 0) {
  tryCatch(getGEOSuppFiles("GSE161969", baseDir = "data"),
           error = function(e) message("Download failed: ", conditionMessage(e)))
}
files <- list.files(supp_dir, full.names = TRUE, recursive = TRUE)
cat("\nFiles in", supp_dir, ":\n"); print(basename(files))
if (length(files) == 0) {
  stop("No files downloaded. Download the supplementary file(s) for GSE161969 manually ",
       "from its GEO page and place them in data/GSE161969/, then re-run section 3b.")
}

cand <- files[grepl("\\.(txt|tsv|csv|xls|xlsx)(\\.gz)?$", files, ignore.case = TRUE) &
                !grepl("RAW\\.tar", files)]
if (!is.null(counts_file_override)) {
  cand <- counts_file_override
} else if (length(cand) == 0) {
  stop("Only a RAW.tar (per-sample files) or unsupported format found.")
} else if (length(cand) > 1) {
  pref <- cand[grepl("count|read", basename(cand), ignore.case = TRUE)]
  if (length(pref) == 0) pref <- cand[grepl("fpkm|tpm|expr|matrix|gene", basename(cand), ignore.case = TRUE)]
  if (length(pref) == 0) pref <- cand
  cand <- pref[1]
}
counts_file <- cand[1]
cat("\nUsing file:", counts_file, "\n")

## 3c. Read matrix and match columns to samples
raw <- read_any(counts_file)
cat("Dimensions:", dim(raw), "\n")

# sample number is the leading number in the title, e.g. "42 (control 1)" -> "42"
meta$num <- sub("^\\s*(\\d+).*$", "\\1", meta$title)
count_cols <- paste0("X", meta$num, "_count")
print(data.frame(title = meta$title, gsm = meta$geo_accession, column = count_cols))

if (!all(count_cols %in% colnames(raw))) {
  stop("Missing columns: ", paste(setdiff(count_cols, colnames(raw)), collapse = ", "))
}

mat <- as.matrix(raw[, count_cols, drop = FALSE])
mode(mat) <- "numeric"
colnames(mat) <- meta$geo_accession

genes <- as.character(raw$gene_name)
ok <- !is.na(genes) & genes != "" & genes != "-" & !grepl("^MSTRG", genes) &
  !is.na(rowSums(mat))
mat <- rowsum(mat[ok, ], genes[ok])          # sum duplicate gene symbols
group <- factor(meta$group, levels = c("ctrl", "case"))

cat("\nGenes in matrix:", nrow(mat), "\n")
print(table(group))
print(setNames(as.character(group), colnames(mat)))   # sanity check of group assignment

## 3d. Differential expression
is_count <- all(abs(mat - round(mat)) < 1e-6)
if (is_count) {
  cat("Integer values detected -> DESeq2\n")
  mat <- round(mat)
  mat <- mat[rowSums(edgeR::cpm(mat) > 1) >= min(table(group)), ]
  dds <- DESeqDataSetFromMatrix(mat, data.frame(group, row.names = colnames(mat)), ~ group)
  dds <- DESeq(dds)
  
  res_raw <- results(dds, contrast = c("group", "case", "ctrl"))
  res161969 <- as.data.frame(res_raw) |>
    tibble::rownames_to_column("gene") |>
    dplyr::rename(log2FC_unshrunk = log2FoldChange)
  
  if (use_shrunk_lfc) {
    coef_name <- "group_case_vs_ctrl"
    stopifnot(coef_name %in% resultsNames(dds))
    shr <- tryCatch(
      lfcShrink(dds, coef = coef_name, type = "apeglm"),
      error = function(e) {
        message("apeglm failed (", conditionMessage(e), ") -> using ashr")
        lfcShrink(dds, coef = coef_name, type = "ashr")
      })
    # keep Wald pvalue/padj from the standard results; take only the LFC from the shrunk fit
    res161969$log2FC <- shr$log2FoldChange[match(res161969$gene, rownames(shr))]
  } else {
    res161969$log2FC <- res161969$log2FC_unshrunk
  }
  expr161969 <- assay(vst(dds, blind = TRUE))
} else {
  cat("Non-integer values -> treating as FPKM/TPM, using limma on log2(x+1)\n")
  expr161969 <- log2(mat + 1)
  expr161969 <- expr161969[rowMeans(expr161969) > 0.5, ]
  fit <- eBayes(lmFit(expr161969, model.matrix(~ group)), trend = TRUE, robust = TRUE)
  res161969 <- tibble::rownames_to_column(topTable(fit, coef = 2, number = Inf), "gene")
  res161969 <- dplyr::rename(res161969, log2FC = logFC, padj = adj.P.Val, pvalue = P.Value)
  res161969$log2FC_unshrunk <- res161969$log2FC
}

res161969 <- res161969[!is.na(res161969$padj) & !is.na(res161969$pvalue) &
                         !is.na(res161969$log2FC), ]
if (basemean_min > 0 && "baseMean" %in% colnames(res161969)) {
  n0 <- nrow(res161969)
  res161969 <- res161969[res161969$baseMean >= basemean_min, ]
  cat("Decidua: dropped", n0 - nrow(res161969), "genes with baseMean <", basemean_min, "\n")
}
cat("Decidua genes tested:", nrow(res161969), "\n")

############################################################
## ---- 4. GSE306548 (placenta, RPL vs Normal) ----
############################################################
if (is.null(placenta_file)) {
  placenta_file <- "data/GSE306548/GSE306548_DESeq2Out_RPL_vs_Normal.csv.gz"
  dir.create(dirname(placenta_file), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(placenta_file)) {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE306nnn/GSE306548/suppl/",
                  basename(placenta_file))
    tryCatch(download.file(url, placenta_file, mode = "wb", method = "libcurl"),
             error = function(e) message("Download failed: ", conditionMessage(e)))
  }
}
if (!file.exists(placenta_file) || file.size(placenta_file) < 1000) {
  stop("Placenta results file not found: ", placenta_file,
       "\nSet placenta_file in SETTINGS to your own full DESeq2 table, or download the ",
       "authors' file from the GEO page (Supplementary files).")
}

x  <- read_any(placenta_file)
nm <- colnames(x)
cat("\nPlacenta file columns:\n"); print(nm)

c_fc <- pick_col(nm, "^log2FoldChange$|^log2FC$|^logFC$|log2.*fold")
c_pa <- pick_col(nm, "^padj$|adj.*p|^FDR$|qval")
c_pv <- pick_col(nm, "^pvalue$|^p[._ ]?value$|^PValue$")
c_gn <- pick_col(nm, "^symbol$|gene_?name|gene.symbol|^gene$")
if (is.na(c_gn)) c_gn <- nm[1]          # blank/first column = gene IDs
if (is.na(c_fc) || is.na(c_pa)) {
  stop("Could not identify log2FC / padj columns in the placenta file; check the column names above.")
}
cat("Using columns -> gene:", c_gn, "| log2FC:", c_fc, "| padj:", c_pa, "| pvalue:", c_pv, "\n")


nm[1] <- "gene"; colnames(x) <- nm
c_gn <- "gene"

res306548 <- data.frame(gene   = to_symbols(x[[c_gn]]),
                        log2FC = as.numeric(x[[c_fc]]),
                        pvalue = if (!is.na(c_pv)) as.numeric(x[[c_pv]]) else NA_real_,
                        padj   = as.numeric(x[[c_pa]]),
                        stringsAsFactors = FALSE)
res306548 <- res306548[!is.na(res306548$gene) & res306548$gene != "" &
                         !is.na(res306548$log2FC) & !is.na(res306548$padj), ]
res306548 <- res306548[order(res306548$padj), ]
res306548 <- res306548[!duplicated(res306548$gene), ]
cat("Placenta genes tested:", nrow(res306548), "\n")
cat("Placenta DEGs (padj<0.05, |log2FC|>=1):",
    sum(res306548$padj < 0.05 & abs(res306548$log2FC) >= 1), "\n")

# NOTE: positive log2FC must mean higher in RPL in BOTH datasets. Decidua is case vs ctrl (RPL vs control);
# for placenta "RPL_vs_Normal" in DESeq2 means positive = higher in RPL. Verify on the GEO page if unsure.

############################################################
## ---- 5. CALL DEGs + VOLCANO PLOTS ----
############################################################
out <- list(GSE161969 = list(res = res161969, expr = expr161969, group = group),
            GSE306548 = list(res = res306548, expr = NULL,       group = NULL))

for (id in names(out)) {
  cfg_ <- deg_cfg[[id]]
  r <- out[[id]]$res
  r$stat_used <- r[[cfg_$stat]]
  r$DEG <- ifelse(r$stat_used < cfg_$p & r$log2FC >=  cfg_$fc, "Up",
                  ifelse(r$stat_used < cfg_$p & r$log2FC <= -cfg_$fc, "Down", "NS"))
  out[[id]]$res <- r
  write.csv(r, file.path("results", paste0(id, "_DE_all.csv")), row.names = FALSE)
  cat("\n", id, "(", cfg_$stat, "<", cfg_$p, ", |log2FC| >=", cfg_$fc, ")\n"); print(table(r$DEG))
  
  p <- ggplot(r, aes(log2FC, -log10(stat_used), colour = DEG)) +
    geom_point(alpha = .6, size = 1) +
    scale_colour_manual(values = c(Up = "firebrick", Down = "steelblue", NS = "grey80")) +
    geom_vline(xintercept = c(-cfg_$fc, cfg_$fc), linetype = 2) +
    geom_hline(yintercept = -log10(cfg_$p), linetype = 2) +
    labs(y = paste0("-log10(", cfg_$stat, ")"), title = id) + theme_bw()
  ggsave(file.path("results", paste0(id, "_volcano.pdf")), p, width = 6, height = 5)
}

############################################################
## ---- 6. OVERLAP: same-direction DEGs + statistical test ----
############################################################
r1 <- out$GSE161969$res
r2 <- out$GSE306548$res

# Background = genes tested in BOTH datasets (this is what makes the test valid)
tested_both <- intersect(r1$gene, r2$gene)
N <- length(tested_both)
r1 <- r1[r1$gene %in% tested_both, ]
r2 <- r2[r2$gene %in% tested_both, ]
cat("\nGenes tested in both datasets (background N):", N, "\n")

g <- function(r, d) r$gene[r$DEG == d]
d1_all <- r1$gene[r1$DEG != "NS"]; d2_all <- r2$gene[r2$DEG != "NS"]

up   <- intersect(g(r1, "Up"),   g(r2, "Up"))
down <- intersect(g(r1, "Down"), g(r2, "Down"))
common <- c(up, down)                                   # SAME-direction genes only
discordant <- c(intersect(g(r1, "Up"),   g(r2, "Down")),
                intersect(g(r1, "Down"), g(r2, "Up")))
cat("Common UP:", length(up), "| Common DOWN:", length(down),
    "| Same-direction total:", length(common),
    "| Discordant (excluded):", length(discordant), "\n")

# Enrichment of the overlap versus chance
ov_stats <- rbind(
  hyper_overlap(d1_all,           d2_all,           N, "Any direction"),
  hyper_overlap(g(r1, "Up"),      g(r2, "Up"),      N, "Up x Up"),
  hyper_overlap(g(r1, "Down"),    g(r2, "Down"),    N, "Down x Down"),
  hyper_overlap(g(r1, "Up"),      g(r2, "Down"),    N, "Up x Down (discordant)"),
  hyper_overlap(g(r1, "Down"),    g(r2, "Up"),      N, "Down x Up (discordant)")
)
print(ov_stats, row.names = FALSE)
write.csv(ov_stats, "results/overlap_hypergeometric_tests.csv", row.names = FALSE)

# Threshold-free agreement across all commonly tested genes
lfc_both <- merge(r1[, c("gene", "log2FC")], r2[, c("gene", "log2FC")], by = "gene",
                  suffixes = c("_decidua", "_placenta"))
rho <- suppressWarnings(cor.test(lfc_both$log2FC_decidua, lfc_both$log2FC_placenta, method = "spearman"))
cat(sprintf("Spearman correlation of log2FC across %d tested genes: rho = %.3f (p = %.2g)\n",
            nrow(lfc_both), rho$estimate, rho$p.value))

# Sensitivity (exploratory, NOT used downstream): decidua at raw p instead of FDR
d1_raw <- r1[r1$pvalue < exploratory_p & abs(r1$log2FC) >= deg_cfg$GSE161969$fc, ]
same_raw <- sum(sign(d1_raw$log2FC[d1_raw$gene %in% d2_all]) ==
                  sign(r2$log2FC[match(d1_raw$gene[d1_raw$gene %in% d2_all], r2$gene)]))
cat(sprintf("Sensitivity: decidua raw p<%.2f -> %d DEGs, %d same-direction overlaps (exploratory only)\n",
            exploratory_p, nrow(d1_raw), same_raw))

# Tables
shared <- data.frame(gene = intersect(d1_all, d2_all))
shared$log2FC_decidua  <- r1$log2FC[match(shared$gene, r1$gene)]
shared$padj_decidua    <- r1$padj[match(shared$gene, r1$gene)]
shared$log2FC_placenta <- r2$log2FC[match(shared$gene, r2$gene)]
shared$padj_placenta   <- r2$padj[match(shared$gene, r2$gene)]
shared$direction <- ifelse(shared$log2FC_decidua > 0 & shared$log2FC_placenta > 0, "Up in both",
                           ifelse(shared$log2FC_decidua < 0 & shared$log2FC_placenta < 0, "Down in both",
                                  "Discordant"))
write.csv(shared, "results/shared_DEGs_all_with_direction.csv", row.names = FALSE)
write.csv(data.frame(gene = common, direction = rep(c("Up", "Down"), c(length(up), length(down)))),
          "results/common_DEGs.csv", row.names = FALSE)
write.csv(shared[shared$direction == "Discordant", ], "results/discordant_DEGs.csv", row.names = FALSE)

ggsave("results/venn_DEGs.pdf",
       ggVennDiagram(list(`161969 Up` = g(r1, "Up"), `306548 Up` = g(r2, "Up"),
                          `161969 Down` = g(r1, "Down"), `306548 Down` = g(r2, "Down"))),
       width = 7, height = 6)

if (length(common) < 10) {
  stop("Fewer than 10 same-direction DEGs. Consider |log2FC| >= 0.585 in deg_cfg and re-run from section 5.")
}

############################################################
## ---- 7. ENRICHMENT (GO + KEGG), background = tested genes ----
############################################################
ids      <- bitr(common,       "SYMBOL", "ENTREZID", org.Hs.eg.db)
universe <- bitr(tested_both,  "SYMBOL", "ENTREZID", org.Hs.eg.db)$ENTREZID

ego <- enrichGO(ids$ENTREZID, universe = universe, OrgDb = org.Hs.eg.db, ont = "ALL",
                pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  write.csv(as.data.frame(ego), "results/GO_enrichment.csv", row.names = FALSE)
  ggsave("results/GO_dotplot.pdf",
         dotplot(ego, split = "ONTOLOGY", showCategory = 8) + facet_grid(ONTOLOGY ~ ., scales = "free"),
         width = 8, height = 10)
} else message("No significant GO terms.")

ekegg <- tryCatch(enrichKEGG(ids$ENTREZID, universe = universe, organism = "hsa", pvalueCutoff = 0.05),
                  error = function(e) { message("KEGG failed (needs internet): ", conditionMessage(e)); NULL })
if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  write.csv(as.data.frame(ekegg), "results/KEGG_enrichment.csv", row.names = FALSE)
  ggsave("results/KEGG_dotplot.pdf", dotplot(ekegg, showCategory = 15), width = 7, height = 6)
} else message("No significant KEGG pathways.")

############################################################
## ---- 8. PPI NETWORK (STRING) ----
############################################################
string_db <- STRINGdb$new(version = "12.0", species = species,
                          score_threshold = 400, network_type = "full", input_directory = "")
mapped <- string_db$map(data.frame(gene = common), "gene", removeUnmappedRows = TRUE)
edges  <- string_db$get_interactions(mapped$STRING_id)
edges$from <- mapped$gene[match(edges$from, mapped$STRING_id)]
edges$to   <- mapped$gene[match(edges$to,   mapped$STRING_id)]
edges <- unique(edges[edges$from != edges$to & !is.na(edges$from) & !is.na(edges$to), ])
write.csv(edges, "results/PPI_edges.csv", row.names = FALSE)    # importable into Cytoscape

net <- simplify(graph_from_data_frame(edges[, c("from", "to")], directed = FALSE))
net <- induced_subgraph(net, which(degree(net) > 0))
cat("\nPPI network:", vcount(net), "nodes,", ecount(net), "edges\n")
if (vcount(net) < 10) stop("PPI network too small. Lower score_threshold to 150 or relax DEG thresholds.")

############################################################
## ---- 9. HUB GENES (cytoHubba-style, top 10) ----
############################################################
mcc <- function(gr) {                      # Maximal Clique Centrality
  cl <- max_cliques(gr, min = 2)
  sc <- setNames(numeric(vcount(gr)), V(gr)$name)
  for (cq in cl) { nmv <- as_ids(cq); sc[nmv] <- sc[nmv] + factorial(length(nmv) - 1) }
  sc
}
mnc <- function(gr) {                      # Maximum Neighborhood Component
  sapply(V(gr)$name, function(v) {
    nb <- neighbors(gr, v)
    if (length(nb) == 0) return(0)
    max(components(induced_subgraph(gr, nb))$csize)
  })
}

metrics <- data.frame(
  gene        = V(net)$name,
  Degree      = degree(net),
  MCC         = as.numeric(mcc(net)[V(net)$name]),
  MNC         = as.numeric(mnc(net)),
  Closeness   = suppressWarnings(closeness(net, normalized = TRUE)),
  Betweenness = betweenness(net, normalized = TRUE),
  stringsAsFactors = FALSE
)
methods <- c("Degree", "MCC", "MNC", "Closeness", "Betweenness")

top10_each <- lapply(methods, function(m) metrics$gene[order(metrics[[m]], decreasing = TRUE)][1:10])
names(top10_each) <- methods

ranks <- sapply(methods, function(m) rank(-metrics[[m]], ties.method = "min"))
metrics$MeanRank <- rowMeans(ranks)
metrics$N_top10  <- sapply(metrics$gene, function(x) sum(sapply(top10_each, function(v) x %in% v)))

hub <- head(metrics[order(-metrics$N_top10, metrics$MeanRank), ], 10)
cat("\n===== TOP 10 HUB GENES =====\n"); print(hub, row.names = FALSE)

write.csv(metrics, "results/network_metrics_all.csv", row.names = FALSE)
write.csv(hub,     "results/top10_hub_genes.csv",      row.names = FALSE)
ggsave("results/hub_methods_venn.pdf", ggVennDiagram(top10_each), width = 8, height = 6)

sub <- induced_subgraph(net, hub$gene)
V(sub)$deg <- degree(sub)
ggsave("results/hub_network.pdf",
       ggraph(sub, layout = "circle") + geom_edge_link(alpha = .4) +
         geom_node_point(aes(size = deg), colour = "firebrick") +
         geom_node_text(aes(label = name), repel = TRUE) + theme_void(),
       width = 6, height = 6)

############################################################
## ---- 10. VALIDATION OF HUB GENES ----
############################################################
hub_genes <- hub$gene

# Direction / effect-size concordance across datasets (FDR in both)
fc_tab <- data.frame(
  gene = hub_genes,
  log2FC_GSE161969 = r1$log2FC[match(hub_genes, r1$gene)],
  log2FC_GSE161969_unshrunk = r1$log2FC_unshrunk[match(hub_genes, r1$gene)],
  log2FC_GSE306548 = r2$log2FC[match(hub_genes, r2$gene)],
  padj_GSE161969   = r1$padj[match(hub_genes, r1$gene)],
  padj_GSE306548   = r2$padj[match(hub_genes, r2$gene)],
  baseMean_GSE161969 = if ("baseMean" %in% colnames(r1)) r1$baseMean[match(hub_genes, r1$gene)] else NA
)
fc_tab$same_direction <- sign(fc_tab$log2FC_GSE161969) == sign(fc_tab$log2FC_GSE306548)
print(fc_tab); write.csv(fc_tab, "results/hub_gene_concordance.csv", row.names = FALSE)

# Heatmap (only GSE161969 has sample-level data)
m <- out$GSE161969$expr[intersect(hub_genes, rownames(out$GSE161969$expr)), , drop = FALSE]
ann <- data.frame(Group = out$GSE161969$group, row.names = colnames(m))
pdf("results/GSE161969_hub_heatmap.pdf", 7, 5)
pheatmap(m, scale = "row", annotation_col = ann, show_colnames = TRUE)
dev.off()

cat("\nDone. All outputs are in the 'results' folder.\n")
