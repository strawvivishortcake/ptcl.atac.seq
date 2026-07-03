#### 0. Packages ####
library(monaLisa)
library(SummarizedExperiment)
library(Biostrings)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2024)
library(RSQLite)
library(DESeq2)
library(motifmatchr)     # NEW - for per-sample motif matching
library(ComplexHeatmap)
library(circlize)
library(dplyr)

data_dir <- "C:/Users/saanv/Downloads/SE_ATAC_count_matrix/SE_ATAC_count_matrix"
fig_dir  <- "C:/Users/saanv/Downloads/SE_ATAC_count_matrix/Figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

genome <- BSgenome.Hsapiens.UCSC.hg38
JASPARConnect <- JASPAR2024::JASPAR2024()
motifs_pfm <- TFBSTools::getMatrixSet(
  RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPARConnect)),
  opts = list(species = 9606, collection = "CORE")
)
motifs <- TFBSTools::toPWM(motifs_pfm)

matrix_files <- c(
  "ALK_neg_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ALK_pos_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ANKL_vs_NK_se_ATAC_count_matrix.tsv",
  "ATLL_vs_CD4_se_ATAC_count_matrix.tsv",
  "BIA_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "CTCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ENKTL_vs_NK_se_ATAC_count_matrix.tsv",
  "nTFHL_vs_TFH_se_ATAC_count_matrix.tsv",
  "PTCL-NOS_vs_CD4_se_ATAC_count_matrix.tsv",
  "T-PLL_vs_TcmCD4_se_ATAC_count_matrix.tsv"
)

ref_prefix_overrides <- c("T-PLL_vs_TcmCD4" = "CD4_Tcm")

#### 1. monaLisa helper (same as your existing runmonaLisa, unchanged) ####
runmonaLisa <- function(DAR, motifs, peaks, genome, nBins = 11, minAbsLfc = 0.3,
                         background = c("otherBins", "zeroBin"), minMatchScore = 10) {
  peak_FC <- jitter(DAR$logFC, 0.001)
  values(peaks) <- DataFrame(peak_FC)
  fc2 <- peak_FC[which(abs(peak_FC) >= minAbsLfc)]
  effectiveBins <- min(nBins, max(3, floor(length(fc2) / 5)))
  if (effectiveBins < nBins) nBins <- effectiveBins
  nElements <- ceiling(length(fc2) / nBins)

  bins <- monaLisa::bin(x = peaks$peak_FC, binmode = "equalN",
                         nElements = nElements, minAbsX = minAbsLfc)
  DARseqs <- getSeq(genome, peaks)
  if (!is.null(names(peaks))) names(DARseqs) <- names(peaks)

  has_N <- vcountPattern("N", DARseqs) > 0
  if (any(has_N)) {
    DARseqs <- DARseqs[!has_N]; peaks <- peaks[!has_N]; bins <- bins[!has_N]
  }

  se <- calcBinnedMotifEnrR(seqs = DARseqs, bins = bins, pwmL = motifs,
                             min.score = minMatchScore,
                             BPPARAM = BiocParallel::MulticoreParam(8),
                             background = match.arg(background))

  simes <- function(pval) min((length(pval) * pval[order(pval)]) / seq_along(pval))
  ML <- se
  zerobin <- which(colData(ML)$bin.nochange)
  ML <- ML[, -zerobin]
  MLp <- 10^-assays(ML)$negLog10P
  MLsimes <- apply(MLp, 1, simes)
  MLdf <- data.frame(p = MLsimes, padj = p.adjust(MLsimes))
  MLdf <- MLdf[order(MLdf$p), ]
  MLdf$rank <- seq_along(row.names(MLdf))
  list(res = se, df = MLdf)
}

#### 2. Per-sample motif score ####
# For each comparison: match every motif against every peak once (motifmatchr),
# then for each sample compute log2( mean normalized accessibility in that
# motif's matched peaks  /  mean normalized accessibility across ALL peaks in
# that comparison ). This mirrors monaLisa's own log2-enrichment logic but
# swaps "bin of peaks" for "one sample's normalized signal".
score_samples_for_comparison <- function(counts_norm, peaks, motifs, p_cutoff = 5e-5) {
  matches <- matchMotifs(motifs, peaks, genome = genome,
                          out = "matches", p.cutoff = p_cutoff)
  match_mat <- motifMatches(matches)              # peaks x motifs, logical

  sample_ids <- colnames(counts_norm)
  motif_names <- colnames(match_mat)

  score_mat <- matrix(NA_real_, nrow = length(motif_names), ncol = length(sample_ids),
                       dimnames = list(motif_names, sample_ids))

  global_mean <- colMeans(counts_norm)             # per-sample mean over ALL peaks

  for (m in motif_names) {
    in_motif <- match_mat[, m]
    if (sum(in_motif) < 5) next                    # skip motifs with too few matched peaks
    motif_mean <- colMeans(counts_norm[in_motif, , drop = FALSE])
    score_mat[m, ] <- log2((motif_mean + 1e-6) / (global_mean + 1e-6))
  }
  score_mat
}

#### 3. Main loop: run monaLisa (for TF significance) + per-sample scoring ####
all_monaLisa_df   <- list()
all_sample_scores <- list()
all_metadata      <- list()

for (file_name in matrix_files) {
  message("Processing: ", file_name)
  comparison  <- gsub("_se_ATAC_count_matrix.tsv", "", file_name)
  group_names <- strsplit(comparison, "_vs_")[[1]]

  counts_full <- read.delim(file.path(data_dir, file_name), row.names = 1, check.names = FALSE)

  ref_prefix <- if (comparison %in% names(ref_prefix_overrides)) {
    ref_prefix_overrides[[comparison]]
  } else group_names[2]

  grp2_cols <- grep(paste0("^", ref_prefix, "_"), colnames(counts_full), value = TRUE)
  grp1_cols <- setdiff(colnames(counts_full), grp2_cols)
  counts <- counts_full[, c(grp1_cols, grp2_cols)]

  if (length(grp1_cols) == 0 || length(grp2_cols) == 0) {
    stop("Column split failed for ", file_name)
  }

  colData <- DataFrame(condition = factor(
    c(rep(group_names[1], length(grp1_cols)), rep(group_names[2], length(grp2_cols))),
    levels = c(group_names[2], group_names[1])))

  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                 colData = colData, design = ~condition)
  dds <- DESeq(dds)
  DAR <- as.data.frame(results(dds, name = resultsNames(dds)[2]))
  DAR$logFC <- DAR$log2FoldChange

  valid <- !is.na(DAR$logFC) & is.finite(DAR$logFC)
  DAR <- DAR[valid, ]
  counts <- counts[valid, ]

  peaks <- GRanges(rownames(counts))
  names(peaks) <- rownames(counts)

  # -- monaLisa: identify significant TFs for this comparison --
  ml <- runmonaLisa(DAR = DAR, motifs = motifs, peaks = peaks, genome = genome)
  all_monaLisa_df[[comparison]] <- ml$df

  # -- per-sample scoring, using DESeq2 size-factor-normalized counts --
  counts_norm <- counts(dds, normalized = TRUE)
  counts_norm <- counts_norm[valid, , drop = FALSE]
  rownames(counts_norm) <- rownames(counts)

  all_sample_scores[[comparison]] <- score_samples_for_comparison(counts_norm, peaks, motifs)

  # -- metadata for annotation columns (reuses your existing logic) --
  for (col_name in colnames(counts)) {
    is_normal <- grepl("CD4|NK|TFH|TcmCD4|control|Normal", col_name, ignore.case = TRUE)
    all_metadata[[col_name]] <- data.frame(
      SampleID    = col_name,
      Subtype     = ifelse(is_normal, "Normal", comparison),
      Compartment = ifelse(is_normal, "Normal", "Disease"),
      SampleType  = ifelse(is_normal, "control", "tumor-cell_line"),
      Batch       = comparison,
      stringsAsFactors = FALSE
    )
  }
}

saveRDS(list(monaLisa = all_monaLisa_df, scores = all_sample_scores, metadata = all_metadata),
        file.path(fig_dir, "combined_pipeline_results.rds"))

#### 4. Assemble the combined matrix ####
metadata <- do.call(rbind, all_metadata)

# Choose which TFs to show: top-N significant per comparison, union across all
top_n_per_comparison <- 15
sig_tfs <- unique(unlist(lapply(all_monaLisa_df, function(df) {
  head(rownames(df)[order(df$padj)], top_n_per_comparison)
})))

# Build one big motif x sample matrix, aligning columns to `metadata`
combined_mat <- matrix(NA_real_, nrow = length(sig_tfs), ncol = nrow(metadata),
                        dimnames = list(sig_tfs, metadata$SampleID))

for (comparison in names(all_sample_scores)) {
  score_mat <- all_sample_scores[[comparison]]
  common_tfs <- intersect(sig_tfs, rownames(score_mat))
  common_samples <- intersect(colnames(score_mat), colnames(combined_mat))
  combined_mat[common_tfs, common_samples] <- score_mat[common_tfs, common_samples]
}

# Some samples (e.g. "Normal" controls) appear once per comparison batch they were
# paired with - if you want ONE de-duplicated Normal column set instead of repeats
# per batch, dedupe metadata/combined_mat by SampleID before this point.

combined_mat[is.na(combined_mat)] <- 0
matrix_scaled <- t(scale(t(combined_mat)))
matrix_scaled[is.na(matrix_scaled)] <- 0
matrix_scaled[matrix_scaled > 3.5]  <- 3.5
matrix_scaled[matrix_scaled < -3.5] <- -3.5

#### 5. TF grouping (reuse your existing tf_groups logic, subset to sig_tfs) ####
tf_groups <- rep("Other", length(sig_tfs))
names(tf_groups) <- sig_tfs
tf_groups[grepl("JUN|FOS|BATF|JDP2|Atf3", sig_tfs, ignore.case = TRUE)] <- "AP-1 (FOS/JUN/ATF/BATF)"
tf_groups[grepl("MAF|BACH|NFE2", sig_tfs, ignore.case = TRUE)]         <- "bZIP MAF/NRF2/BACH"
tf_groups[grepl("REL|NFKB", sig_tfs, ignore.case = TRUE)]              <- "NF-kB"
tf_groups[grepl("STAT|IRF", sig_tfs, ignore.case = TRUE)]              <- "JAK-STAT"
tf_groups[grepl("CEBP|ROR|TP53|BCL11|NR2|Nr1|HOX", sig_tfs, ignore.case = TRUE)] <- "Nuclear receptor/CEBP"
tf_groups[grepl("MYC|MAX|SNAI|ZEB|TCF|Ptf1A|FIGLA|TEAD|MNT|Mlxip|Arnt", sig_tfs, ignore.case = TRUE)] <- "bHLH/MYC/EMT"
tf_groups[grepl("GATA|TRPS1", sig_tfs, ignore.case = TRUE)]            <- "GATA"
tf_groups[grepl("CTCF", sig_tfs, ignore.case = TRUE)]                  <- "CTCF/insulator"
tf_groups[grepl("KLF|SP|ZNF148|ZNF281", sig_tfs, ignore.case = TRUE)]  <- "KLF/SP"
tf_groups[grepl("RFX|ZNF175|ZNF384|ZNF362|ZFP28|Hnf1A|NFATC", sig_tfs, ignore.case = TRUE)] <- "RFX / Zinc Finger"
tf_groups[grepl("RUNX|TBR|TBX|EOMES", sig_tfs, ignore.case = TRUE)]    <- "RUNX/CBF & T-box"
tf_groups[grepl("ETS|ETV|EHF|ELF|ELK|GABPA|ZBTB11|ERF|Erg", sig_tfs, ignore.case = TRUE)] <- "ETS family"
tf_groups[grepl("FOX|IKZF|Spi1|SPIB", sig_tfs, ignore.case = TRUE)]    <- "Forkhead / Lymphoid ZF"

#### 6. Annotations + heatmap (same styling as your mock-up script) ####
col_compartment <- c("Disease" = "#D95F02", "Normal" = "#1B9E77")
col_sampletype  <- c("control" = "#7570B3", "tumor-cell_line" = "#E7298A")

unique_subtypes <- unique(metadata$Subtype)
col_subtype <- setNames(RColorBrewer::brewer.pal(max(3, length(unique_subtypes)), "Set1")[1:length(unique_subtypes)], unique_subtypes)
col_subtype["Normal"] <- "#1B9E77"

unique_batches <- unique(metadata$Batch)
col_batch <- setNames(RColorBrewer::brewer.pal(max(3, length(unique_batches)), "Pastel1")[1:length(unique_batches)], unique_batches)

top_ann <- HeatmapAnnotation(
  Compartment = metadata$Compartment,
  SampleType  = metadata$SampleType,
  Subtype     = metadata$Subtype,
  Batch       = metadata$Batch,
  col = list(Compartment = col_compartment, SampleType = col_sampletype,
             Subtype = col_subtype, Batch = col_batch),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold")
)

ht_plot <- Heatmap(
  matrix_scaled,
  name = "log2 enrich",
  col = colorRamp2(c(-3.5, 0, 3.5), c("#313695", "white", "#a50026")),
  column_split = metadata$Subtype,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  top_annotation = top_ann,
  column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
  row_split = factor(tf_groups),
  cluster_rows = TRUE,
  show_row_names = TRUE,
  row_title_rot = 0,
  row_names_gp = gpar(fontsize = 5),
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  border = TRUE,
  gap = unit(1.2, "mm"),
  row_gap = unit(0.8, "mm")
)

pdf(file.path(fig_dir, "combined_persample_TF_heatmap.pdf"), width = 15, height = 14)
draw(ht_plot, annotation_legend_side = "right", heatmap_legend_side = "right")
dev.off()
