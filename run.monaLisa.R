#Step 1: Loading Libraries
library(monaLisa)
library(SummarizedExperiment)
library(Biostrings)
library(TFBSTools)
#Step 2: Run monaLisa
runmonaLisa <- function(DAR, motifs, peaks, genome, nBins=11, minAbsLfc=0.3,
                        background = c("otherBins","zeroBin"), minMatchScore=10){
  
  # Inject small jitter to avoid non-unique breaks during binning
  peak_FC <- jitter(DAR$logFC, 0.001) 
  values(peaks) <- DataFrame(peak_FC)

  fc2 <- peak_FC[which(abs(peak_FC) >= minAbsLfc)]
  nElements <- ceiling(length(fc2) / nBins)
  
  ptm <- proc.time()
  
  # Bin peaks by Log Fold Change
  bins <- monaLisa::bin(x = peaks$peak_FC, 
                        binmode = "equalN", 
                        nElements = nElements,
                        minAbsX = minAbsLfc)
  
  DARseqs <- getSeq(genome, peaks)
  if(!is.null(names(peaks))) names(DARseqs) <- names(peaks)
  
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
matrix_files <- c(
  "ALK_neg_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ALK_pos_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ANKL_vs_NK_se_ATAC_count_matrix.tsv",
  "ATLL_vs_CD4_se_ATAC_count_matrix.tsv",
  "BIA_ALCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "CTCL_vs_CD4_se_ATAC_count_matrix.tsv",
  "ENKTL_vs_NK_se_ATAC_count_matrix.tsv",
  "nTFHL_vs_TFH_se_ATAC_count_matrix.tsv",
  "T-PLL_vs_TcmCD4_se_ATAC_count_matrix.tsv"
)
# Check and initialize output directory for your figures
if(!dir.exists("Figures")) dir.create("Figures")

# Execution loop
output_pdf_name <- paste0(gsub("_se_ATAC_count_matrix.tsv", "", file_name), "_heatmap.pdf")

#Step 4: Filter for top 100 significance
# Quantify max row significance while safely suppressing and resolving un-calculated row edge-cases
  max_sig_per_row <- suppressWarnings(apply(assays(se)$negLog10Padj, 1, max, na.rm = TRUE))
  max_sig_per_row[is.infinite(max_sig_per_row)] <- 0 
  
  # Sort indices descending to extract the top 100
  top_100_indices <- order(max_sig_per_row, decreasing = TRUE)[1:100]
  seSel <- se[top_100_indices, ]

#Step 5: Recalculate Structural Distances
# Extract structural PFMs for the 100 isolated transcription factors
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)
  
  # Hierarchical clustering (1 - Pearson correlation matrix)
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
  message("SUCCESSFULLY COMPLETED: ", output_pdf_name)
}
