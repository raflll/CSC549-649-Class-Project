# Part 2: Statistics
# Which genes differ between T2D and non-diabetic (ND) islets?
#
# Run from Part2_statistics/:
#   Rscript scripts/part2_statistics.R
#
# Inputs  (Part1): gene_count_matrix.csv, Table1_donors.csv, mapping_table.csv
# Outputs: results/ and figures/

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(limma)
})

set.seed(1)
in_dir  <- "../Part1_data_pipeline"
res_dir <- "results"
fig_dir <- "figures"
dir.create(res_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

PADJ <- 0.05   # FDR cut-off for the up/down lists

# ---------------------------------------------------------------- data
counts_raw <- read.csv(file.path(in_dir, "gene_count_matrix.csv"),
                       check.names = FALSE)
donors  <- read.csv(file.path(in_dir, "Table1_donors.csv"))
mapping <- read.csv(file.path(in_dir, "mapping_table.csv"))

meta <- merge(donors, mapping[, c("run", "flowcell", "overall_pct", "unique_pct")],
              by = "run")
meta$batch <- factor(ifelse(meta$flowcell == "C11YEACXX", "A", "B"))
meta$group <- factor(meta$group, levels = c("ND", "T2D"))   # ND = reference
meta$sex   <- factor(meta$sex)
rownames(meta) <- meta$run

# gene ids: "ENSG...|SYMBOL" or "ENSG..." only
ids     <- counts_raw$gene_id
ensembl <- sub("\\|.*", "", ids)
symbol  <- ifelse(grepl("\\|", ids), sub("^[^|]*\\|", "", ids), NA)
counts  <- as.matrix(counts_raw[, -1])
rownames(counts) <- ensembl
mode(counts) <- "integer"
genes <- data.frame(ensembl = ensembl, symbol = symbol, row.names = ensembl)

meta   <- meta[colnames(counts), ]                 # same order as the matrix
stopifnot(identical(rownames(meta), colnames(counts)))
label  <- paste0(meta$pair, "_", meta$group)       # e.g. P1_T2D
names(label) <- meta$run

# ---------------------------------------------------------------- QC
lib  <- colSums(counts)
cpm  <- t(t(counts) / lib) * 1e6
sym_cpm <- function(s) {
  i <- which(genes$symbol %in% s)
  if (length(i) == 0) return(rep(0, ncol(cpm)))
  colSums(cpm[i, , drop = FALSE])
}

# purity: endocrine hormones vs exocrine (acinar) and ductal markers
endocrine <- c("INS", "GCG", "SST", "PPY", "GHRL")
acinar    <- c("PRSS1", "PRSS2", "CPA1", "CPA2", "CELA3A", "CELA3B", "CTRB1",
               "CTRB2", "PNLIP", "CPB1", "CEL", "CLPS", "REG1A")
ductal    <- c("KRT19", "KRT7", "CFTR", "SPP1")

qc <- data.frame(
  run = meta$run, label = label[meta$run], pair = meta$pair, group = meta$group,
  sex = meta$sex, batch = meta$batch, library_size = lib,
  INS_cpm = sym_cpm("INS"), GCG_cpm = sym_cpm("GCG"), SST_cpm = sym_cpm("SST"),
  endocrine_pct = sym_cpm(endocrine) / 1e4,
  acinar_pct    = sym_cpm(acinar) / 1e4,
  ductal_pct    = sym_cpm(ductal) / 1e4
)
qc$endo_to_exo <- qc$endocrine_pct / (qc$acinar_pct + qc$ductal_pct)

# sex check: XIST (female) vs chrY genes (male)
ygenes <- c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY")
qc$XIST_cpm <- sym_cpm("XIST")
qc$chrY_cpm <- sym_cpm(ygenes)
qc$sex_from_expr <- ifelse(qc$XIST_cpm > qc$chrY_cpm, "Female", "Male")
qc$sex_match <- qc$sex_from_expr == as.character(qc$sex)

# low purity flag: endocrine share below half the median of all donors
qc$low_purity <- qc$endocrine_pct < 0.5 * median(qc$endocrine_pct)
write.csv(qc, file.path(res_dir, "qc_islet_purity_sex.csv"), row.names = FALSE)

qc_long <- rbind(
  data.frame(label = qc$label, group = qc$group, type = "endocrine (INS, GCG, SST, PPY, GHRL)", pct = qc$endocrine_pct),
  data.frame(label = qc$label, group = qc$group, type = "acinar (PRSS1, CPA1, CELA3A, ...)",   pct = qc$acinar_pct),
  data.frame(label = qc$label, group = qc$group, type = "ductal (KRT19, KRT7, CFTR, SPP1)",    pct = qc$ductal_pct))
qc_long$label <- factor(qc_long$label, levels = qc$label[order(qc$group, qc$pair)])
p <- ggplot(qc_long, aes(label, pct, fill = type)) +
  geom_col() +
  facet_grid(. ~ group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("#7f7f7f", "#c9a227", "#2c7fb8")) +
  labs(x = NULL, y = "% of all reads (CPM / 10,000)", fill = NULL,
       title = "Islet purity: share of reads from endocrine vs exocrine marker genes") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom",
        legend.direction = "vertical")
ggsave(file.path(fig_dir, "qc_purity.png"), p, width = 8, height = 5.5, dpi = 200)

p <- ggplot(qc, aes(log10(XIST_cpm + 1), log10(chrY_cpm + 1), colour = sex, shape = group)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = label), size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(Female = "#d95f02", Male = "#1b9e77")) +
  labs(x = "log10(XIST CPM + 1)", y = "log10(chrY genes CPM + 1)",
       colour = "reported sex", shape = "group",
       title = "Sex check: reported sex vs XIST and chrY expression") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "qc_sex.png"), p, width = 6.5, height = 5, dpi = 200)

# ---------------------------------------------------------------- DESeq2
# keep genes with >= 10 reads in >= 8 donors (8 = size of the smallest group)
keep <- rowSums(counts >= 10) >= 8
dds <- DESeqDataSetFromMatrix(counts[keep, ], colData = meta,
                              design = ~ batch + sex + group)
dds <- DESeq(dds, quiet = TRUE)

# ---------------------------------------------------------------- PCA
vsd <- vst(dds, blind = TRUE)
pca_plot <- function(mat, title, file) {
  top <- head(order(apply(mat, 1, var), decreasing = TRUE), 500)
  pc  <- prcomp(t(mat[top, ]))
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  d   <- data.frame(pc$x[, 1:3], meta, label = label[rownames(pc$x)])
  p <- ggplot(d, aes(PC1, PC2, colour = group, shape = batch)) +
    geom_point(size = 3.5) +
    geom_text_repel(aes(label = label), size = 3, show.legend = FALSE) +
    scale_colour_manual(values = c(ND = "#2c7fb8", T2D = "#d7301f")) +
    labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"),
         title = title, subtitle = "VST counts, 500 most variable genes") +
    theme_bw(base_size = 11)
  ggsave(file, p, width = 7, height = 5.5, dpi = 200)
  list(pc = pc, ve = ve, d = d)
}
pca_raw <- pca_plot(assay(vsd), "PCA of the 16 donors",
                    file.path(fig_dir, "pca.png"))
vst_bc <- removeBatchEffect(assay(vsd), batch = meta$batch, batch2 = meta$sex,
                            design = model.matrix(~ group, meta))
pca_bc <- pca_plot(vst_bc, "PCA after removing batch and sex (limma)",
                   file.path(fig_dir, "pca_batch_sex_removed.png"))

# how much does each PC relate to group / batch / sex / purity?
pc_assoc <- do.call(rbind, lapply(1:5, function(k) {
  x <- pca_raw$pc$x[, k]
  data.frame(PC = paste0("PC", k), var_pct = pca_raw$ve[k],
             p_group = t.test(x ~ meta$group)$p.value,
             p_batch = t.test(x ~ meta$batch)$p.value,
             p_sex   = t.test(x ~ meta$sex)$p.value,
             r_endocrine_pct = cor(x, qc$endocrine_pct[match(rownames(pca_raw$pc$x), qc$run)]))
}))
write.csv(pc_assoc, file.path(res_dir, "pca_associations.csv"), row.names = FALSE)

# ---------------------------------------------------------------- DE results
res  <- results(dds, contrast = c("group", "T2D", "ND"), alpha = PADJ)
shr  <- lfcShrink(dds, coef = "group_T2D_vs_ND", type = "apeglm",
                 apeMethod = "nbinomC*", quiet = TRUE)   # default optimiser left 101 genes unfitted
out  <- data.frame(genes[rownames(res), ],
                   baseMean = res$baseMean,
                   log2FoldChange = res$log2FoldChange,
                   lfcSE = res$lfcSE,
                   log2FoldChange_shrunk = shr$log2FoldChange,
                   stat = res$stat, pvalue = res$pvalue, padj = res$padj)
out <- out[order(out$pvalue), ]
write.csv(out, file.path(res_dir, "de_results.csv"), row.names = FALSE, na = "")

sig  <- subset(out, !is.na(padj) & padj < PADJ)
up   <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)
name_of <- function(d) ifelse(is.na(d$symbol), d$ensembl, d$symbol)
writeLines(name_of(up),   file.path(res_dir, "up_genes.txt"))
writeLines(name_of(down), file.path(res_dir, "down_genes.txt"))
writeLines(name_of(out),  file.path(res_dir, "expressed_genes.txt"))  # background
write.csv(up,   file.path(res_dir, "up_genes.csv"),   row.names = FALSE, na = "")
write.csv(down, file.path(res_dir, "down_genes.csv"), row.names = FALSE, na = "")

# ---------------------------------------------------------------- volcano
v <- out[!is.na(out$padj), ]
v$class <- "not significant"
v$class[v$padj < PADJ & v$log2FoldChange > 0] <- "up in T2D"
v$class[v$padj < PADJ & v$log2FoldChange < 0] <- "down in T2D"
v$name <- name_of(v)
lab <- head(v[v$padj < PADJ, ], 20)
if (nrow(lab) < 10) lab <- head(v, 10)        # label something even if few hits
p <- ggplot(v, aes(log2FoldChange, -log10(pvalue), colour = class)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text_repel(data = lab, aes(label = name), size = 3, colour = "black",
                  max.overlaps = 30, min.segment.length = 0) +
  scale_colour_manual(values = c("up in T2D" = "#d7301f", "down in T2D" = "#2c7fb8",
                                 "not significant" = "grey70")) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(x = "log2 fold change (T2D / ND)", y = "-log10(p value)", colour = NULL,
       title = "Differential expression, T2D vs ND islets",
       subtitle = sprintf("DESeq2, ~ batch + sex + group; %d up, %d down at FDR < %.2f",
                          nrow(up), nrow(down), PADJ)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "volcano.png"), p, width = 7, height = 6, dpi = 200)

png(file.path(fig_dir, "ma_plot.png"), width = 1400, height = 1100, res = 200)
DESeq2::plotMA(shr, ylim = c(-4, 4), alpha = PADJ,
       main = "MA plot, shrunken log2 fold change (apeglm)")
dev.off()

# ---------------------------------------------------------------- sensitivity
# 1) the model recommended in Part 1 (no sex term)
# 2) also adjusting for islet purity (endocrine share, scaled)
dds2 <- DESeq(DESeqDataSetFromMatrix(counts[keep, ], meta, ~ batch + group), quiet = TRUE)
meta$purity <- as.numeric(scale(qc$endocrine_pct[match(meta$run, qc$run)]))
dds3 <- DESeq(DESeqDataSetFromMatrix(counts[keep, ], meta, ~ batch + sex + purity + group),
              quiet = TRUE)
# 3) without the two batch-B T2D donors (P7_T2D, P8_T2D), which drive the top hits
drop <- meta$run %in% c("SRR957965", "SRR957967")
dds4 <- DESeq(DESeqDataSetFromMatrix(counts[keep, !drop], droplevels(meta[!drop, ]),
                                     ~ batch + sex + group), quiet = TRUE)
n_sig <- function(r) sum(!is.na(r$padj) & r$padj < PADJ)
r2 <- results(dds2, contrast = c("group", "T2D", "ND"), alpha = PADJ)
r3 <- results(dds3, contrast = c("group", "T2D", "ND"), alpha = PADJ)
s1 <- rownames(sig)
s2 <- rownames(r2)[which(r2$padj < PADJ)]
s3 <- rownames(r3)[which(r3$padj < PADJ)]
r4 <- results(dds4, contrast = c("group", "T2D", "ND"), alpha = PADJ)
s4 <- rownames(r4)[which(r4$padj < PADJ)]
sens <- data.frame(
  model = c("~ batch + sex + group (main)", "~ batch + group", "~ batch + sex + purity + group",
            "main, without P7_T2D and P8_T2D"),
  n_fdr05 = c(length(s1), length(s2), length(s3), length(s4)),
  overlap_with_main = c(length(s1), length(intersect(s1, s2)), length(intersect(s1, s3)),
                        length(intersect(s1, s4))))
loo <- data.frame(genes[s1, ], padj_main = res[s1, "padj"],
                  padj_batch_group = r2[s1, "padj"], padj_with_purity = r3[s1, "padj"],
                  padj_without_P7P8_T2D = r4[s1, "padj"], pvalue_without_P7P8_T2D = r4[s1, "pvalue"])
write.csv(loo, file.path(res_dir, "sensitivity_top_genes.csv"), row.names = FALSE, na = "")
write.csv(sens, file.path(res_dir, "sensitivity_models.csv"), row.names = FALSE)

# genes reported in earlier T2D islet studies (Fadista 2014 and others)
known <- c("IAPP", "SLC2A2", "CHL1", "LRFN2", "RASGRP1", "PPP1R1A", "CAPN13",
           "GLRA1", "FFAR4", "ARG2", "CCDC141", "SLC30A8", "TCF7L2", "PDX1",
           "MAFA", "INS", "GCG", "G6PC2", "ABCC8", "UNC5D", "SFRP4")
kn <- out[out$symbol %in% known, c("symbol", "baseMean", "log2FoldChange", "pvalue", "padj")]
kn <- kn[order(kn$pvalue), ]
write.csv(kn, file.path(res_dir, "known_t2d_genes.csv"), row.names = FALSE, na = "")

# ---------------------------------------------------------------- summary
sink(file.path(res_dir, "summary.txt"))
cat("Genes in matrix:", nrow(counts), "\n")
cat("Genes tested (>= 10 reads in >= 8 donors):", sum(keep), "\n")
cat("FDR <", PADJ, ":", nrow(sig), "genes (", nrow(up), "up,", nrow(down), "down in T2D )\n")
cat("FDR < 0.1:", sum(out$padj < 0.1, na.rm = TRUE), "\n")
cat("Raw p < 0.01:", sum(out$pvalue < 0.01, na.rm = TRUE), "\n\n")
cat("Top 25 genes by p value:\n")
print(head(out[, c("symbol", "baseMean", "log2FoldChange", "log2FoldChange_shrunk", "pvalue", "padj")], 25),
      row.names = FALSE, digits = 3)
cat("\nPCA variance (%):", pca_raw$ve[1:5], "\n")
print(pc_assoc, row.names = FALSE, digits = 3)
cat("\nQC:\n")
print(qc[, c("label", "batch", "library_size", "endocrine_pct", "acinar_pct", "ductal_pct",
             "sex", "sex_from_expr", "sex_match", "low_purity")], row.names = FALSE, digits = 3)
cat("\nSensitivity:\n"); print(sens, row.names = FALSE); print(loo, row.names = FALSE, digits = 3)
cat("\nKnown T2D islet genes:\n"); print(kn, row.names = FALSE, digits = 3)
cat("\n"); print(sessionInfo())
sink()
cat("done\n")
