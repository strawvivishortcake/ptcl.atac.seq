#Step 1: Load Libraries
library(monaLisa)
library(SummarizedExperiment)
library(Biostrings)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)  # P1: genome sequence source
library(JASPAR2024)                    # P1: motif PWM source
library(RSQLite)
library(DESeq2)                        # needed to compute logFC from count matrices (P6)

#Step 2: Run monaLisa
runmonaLisa <- function(DAR, motifs, peaks, genome, nBins=11, minAbsLfc=0.3,
                        background = c("otherBins","zeroBin"), minMatchScore=10){

  # Inject small jitter to avoid non-unique breaks during binning
  peak_FC <- jitter(DAR$logFC, 0.001)
  values(peaks) <- DataFrame(peak_FC)

  fc2 <- peak_FC[which(abs(peak_FC) >= minAbsLfc)]

  # Guard against comparisons with very few significant peaks (e.g. T-PLL,
  # n=3 vs 3): with a fixed nBins, edge bins can end up with too few peaks
  # for reliable GC-content background matching, which triggers NA errors
  # downstream. Cap nBins so each bin has at least 5 peaks.
  effectiveBins <- min(nBins, max(3, floor(length(fc2) / 5)))
  if (effectiveBins < nBins) {
    message("  Only ", length(fc2), " peaks pass minAbsLfc=", minAbsLfc,
            " - reducing nBins from ", nBins, " to ", effectiveBins)
    nBins <- effectiveBins
  }
  nElements <- ceiling(length(fc2) / nBins)

  ptm <- proc.time()

  # Bin peaks by Log Fold Change
  bins <- monaLisa::bin(x = peaks$peak_FC,
                        binmode = "equalN",
                        nElements = nElements,
                        minAbsX = minAbsLfc)

  DARseqs <- getSeq(genome, peaks)
  if(!is.null(names(peaks))) names(DARseqs) <- names(peaks)

  # Drop any peaks whose sequence contains N bases (unmapped/gap regions) -
  # these break calcBinnedMotifEnrR's GC-content background weighting and
  # cause "NAs are not allowed in subscripted assignments".
  has_N <- vcountPattern("N", DARseqs) > 0
  if (any(has_N)) {
    message("  Dropping ", sum(has_N), " of ", length(DARseqs), " peaks with N bases in sequence")
    DARseqs <- DARseqs[!has_N]
    peaks   <- peaks[!has_N]
    bins    <- bins[!has_N]
  }

  # Calculate binned motif enrichment
  se <- calcBinnedMotifEnrR(seqs = DARseqs,
                            bins = bins,
                            pwmL = motifs,
                            min.score = minMatchScore,
                            BPPARAM = BiocParallel::MulticoreParam(8),
                            background = match.arg(background))

  # Simes method p-value aggregation function
  simes <- function(pval){
    min((length(pval) * pval[order(pval)]) / seq_along(pval))
  }

  ML <- se
  zerobin <- which(colData(ML)$bin.nochange)
  ML <- ML[, -zerobin]
  MLp <- 10^-assays(ML)$negLog10P

  MLsimes <- apply(MLp, 1, simes)

  # Structure output data frame with global adjusted p-values (FDR)
  MLdf <- data.frame(p = MLsimes, padj = p.adjust(MLsimes))
  MLdf <- MLdf[order(MLdf$p), ]
  MLdf$rank <- seq_along(row.names(MLdf))

  # Track trajectory trend using directional Spearman correlation across bins
  cors <- cor(t(assays(ML)$log2enr), seq_len(ncol(ML)), method = "spearman")[, 1]
  names(cors) <- row.names(ML)
  MLdf$binSpearman <- cors[row.names(MLdf)]

  runtime <- proc.time() - ptm

  return(list(res = se, runtime = runtime, df = MLdf))
}

#Step 3: Looping Files

# P2: absolute paths (R needs forward slashes, or escaped backslashes, on Windows)
data_dir <- "C:/Users/saanv/Downloads/SE_ATAC_count_matrix/SE_ATAC_count_matrix"
fig_dir  <- "C:/Users/saanv/Downloads/SE_ATAC_count_matrix/Figures"

matrix_files <- c(
  "ALK_neg_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ALK_pos_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ANKL_vs_NK_se_ATAC_count_matrix.tsv",
  "ATLL_vs_CD4_se_ATAC_count_matrix.tsv",
  "BIA_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "CTCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ENKTL_vs_NK_se_ATAC_count_matrix.tsv",
  "nTFHL_vs_TFH_se_ATAC_count_matrix.tsv",
  "PTCL-NOS_vs_CD4_se_ATAC_count_matrix.tsv",  # P3: was missing - confirm the correct comparator group for PTCL-NOS
  "T-PLL_vs_TcmCD4_se_ATAC_count_matrix.tsv"
)

# Some files don't name their reference-group columns after the group name
# in the "X_vs_Y" filename (e.g. T-PLL_vs_TcmCD4's reference columns are
# actually prefixed "CD4_Tcm_", not "TcmCD4_"). Add overrides here as they
# turn up rather than relying on the filename pattern alone.
ref_prefix_overrides <- c(
  "T-PLL_vs_TcmCD4" = "CD4_Tcm"
  # add more as needed, e.g.:
  # "nTFHL_vs_TFH" = "actual_prefix_here",
  # "ANKL_vs_NK"   = "actual_prefix_here"
)

# P4: load motifs + genome ONCE, outside the loop
if(!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

genome <- BSgenome.Hsapiens.UCSC.hg38
JASPARConnect <- JASPAR2024::JASPAR2024()
motifs_pfm <- TFBSTools::getMatrixSet(
  RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPARConnect)),
  opts = list(species = 9606, collection = "CORE")
)
motifs <- TFBSTools::toPWM(motifs_pfm)  # calcBinnedMotifEnrR requires PWMatrixList, not PFMatrixList

# P5/P6: the per-file body now runs inside a for loop, and actually reads the
# count matrix, computes logFC via DESeq2, and builds the peaks GRanges before
# calling runmonaLisa()
for (file_name in matrix_files) {

  message("Processing: ", file_name)
  output_pdf_name <- file.path(fig_dir,
                    paste0(gsub("_se_ATAC_count_matrix.tsv", "", file_name), "_heatmap.pdf"))

  # --- Read count matrix ---
  # Expected format: first column = peak ID as "chr:start-end", remaining
  # columns = per-sample raw counts, with column-name suffixes identifying
  # the two groups being compared (the "X_vs_Y" pair in the file name).
  counts_full <- read.delim(file.path(data_dir, file_name), row.names = 1, check.names = FALSE)

  comparison  <- gsub("_se_ATAC_count_matrix.tsv", "", file_name)
  group_names <- strsplit(comparison, "_vs_")[[1]]   # e.g. c("ALK_neg_ALCL", "CD4")

  # Only the reference/normal group is consistently prefixed by its own name
  # (e.g. "CD4_36", "CD4_37"). Disease samples are named by cell line /
  # sample ID instead (e.g. "DL-40_5", "FEPD_7"), so we identify them as
  # "everything that is NOT a reference-group column". Some files use a
  # different reference prefix than the filename implies - see
  # ref_prefix_overrides above.
  ref_prefix <- if (comparison %in% names(ref_prefix_overrides)) {
    ref_prefix_overrides[[comparison]]
  } else {
    group_names[2]
  }
  grp2_cols <- grep(paste0("^", ref_prefix, "_"), colnames(counts_full), value = TRUE)
  grp1_cols <- setdiff(colnames(counts_full), grp2_cols)
  counts <- counts_full[, c(grp1_cols, grp2_cols)]

  # Diagnostic: confirm the split looks right before running DESeq2
  message("  ", group_names[1], " (n=", length(grp1_cols), "): ",
          paste(grp1_cols, collapse = ", "))
  message("  ", group_names[2], " (n=", length(grp2_cols), "): ",
          paste(grp2_cols, collapse = ", "))
  if (length(grp1_cols) == 0 || length(grp2_cols) == 0) {
    stop("One of the groups has 0 columns for ", file_name,
        " - check that group_names[2] = '", group_names[2],
        "' actually matches a column prefix in this file.")
  }

  colData <- DataFrame(
    condition = factor(c(rep(group_names[1], length(grp1_cols)),
                          rep(group_names[2], length(grp2_cols))),
                        levels = c(group_names[2], group_names[1]))  # group2 = reference
  )

  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(counts)),
                                colData = colData,
                                design = ~ condition)
  dds <- DESeq(dds)
  DAR <- as.data.frame(results(dds, name = resultsNames(dds)[2]))
  DAR$logFC <- DAR$log2FoldChange

  # Drop peaks with NA/Inf logFC (common for low-count peaks, especially in
  # smaller-n comparisons like T-PLL) - these otherwise propagate NAs into
  # calcBinnedMotifEnrR's background GC-weight correction and crash it.
  valid <- !is.na(DAR$logFC) & is.finite(DAR$logFC)
  message("  Dropping ", sum(!valid), " of ", nrow(DAR), " peaks with NA/Inf logFC")
  DAR <- DAR[valid, ]
  counts <- counts[valid, ]

  # --- Build peaks GRanges from row names ("chr:start-end") ---
  peaks <- GRanges(rownames(counts))
  names(peaks) <- rownames(counts)

  # --- Run monaLisa ---
  result <- runmonaLisa(DAR = DAR, motifs = motifs, peaks = peaks, genome = genome)
  se <- result$res

  #Step 4: Filter for top 100 significance
  max_sig_per_row <- suppressWarnings(apply(assays(se)$negLog10Padj, 1, max, na.rm = TRUE))
  max_sig_per_row[is.infinite(max_sig_per_row)] <- 0

  top_100_indices <- order(max_sig_per_row, decreasing = TRUE)[1:100]
  seSel <- se[top_100_indices, ]

  #Step 5: Recalculate Structural Distances
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  hcl <- hclust(as.dist(1 - SimMatSel), method = "average")

  #Step 6: Save to PDF
  pdf(output_pdf_name, width = 10, height = 15)
  plotMotifHeatmaps(x = seSel,
                    which.plots = c("log2enr", "negLog10Padj"),
                    width = 1.8,
                    cluster = hcl,
                    maxEnr = 2,
                    maxSig = 10,
                    show_dendrogram = TRUE,
                    show_seqlogo = TRUE,
                    width.seqlogo = 1.2)
  dev.off()
}
