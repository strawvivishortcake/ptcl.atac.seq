library(data.table)
library(ComplexHeatmap)
library(circlize)
#Step 1: Defined paths
files <- c(
  ALK_neg = "C:/Users/Public/ALK_neg_ALCL_vs_CD4_se_ATAC_count_HOMER_TF/ALK_neg_knownResults.txt",
  ALK_pos = "C:/Users/Public/ALK_pos_ALCL_vs_CD4_se_ATAC_count_HOMER_TF/ALK_pos_knownResults.txt",
  ANKL    = "C:/Users/Public/ANKL_vs_NK_se_ATAC_count_HOMER_TF/ANKL_knownResults.txt",
  ATLL    = "C:/Users/Public/ATLL_vs_CD4_se_ATAC_count_HOMER_TF/ATLL_knownResults.txt",
  BIA     = "C:/Users/Public/BIA_ALCL_vs_CD4_se_ATAC_count_HOMER_TF/BIA_knownResults.txt",
  CTCL    = "C:/Users/Public/CTCL_vs_CD4_se_ATAC_count_HOMER_TF/CTCL_knownResults.txt",
  ENKTL   = "C:/Users/Public/ENKTL_vs_NK_se_ATAC_count_HOMER_TF/ENKTL_knownResults.txt",
  nTFHL   = "C:/Users/Public/nTFHL_vs_TFH_se_ATAC_count_HOMER_TF/nTFHL_knownResults.txt",
  TPLL    = "C:/Users/Public/T-PLL_vs_TcmCD4_se_ATAC_count_HOMER_TF/TPLL_knownResults.txt"
)
#Helper function: Convert character tags to logical binary metrics
parse_binary_hits <- function(x) {
  as.integer(x != "" & !is.na(x))
}
#Step 2: Extract
process_subtype_to_model <- function(file_path, subtype_name, top_n = 25) {
  if (!file.exists(file_path)) {
    cat("  Warning: File path does not exist for", subtype_name, "- Skipping.\n")
    return(NULL)
  }
  
  #Load peak file
  df <- fread(file_path, header = TRUE, sep = "\t")
  total_peaks <- nrow(df)

  # Grab all motif columns by searching for Distance From Peak
  motif_cols <- grep("Distance From Peak", colnames(df), value = TRUE)
  if (length(motif_cols) == 0) {
    cat("  Warning: No motif columns found for", subtype_name, "- Skipping.\n")
    return(NULL)
  }
  cat("Found", length(motif_cols), "motifs across", total_peaks)

  # Convert columns into 1 or 0
  motif_mat <- sapply(df[, ..motif_cols], parse_binary_hits)
  motif_mat <- as.matrix(motif_mat)
  
  #Isolate the bare motif labels and remove unnecessary quotes
  clean_names <- gsub("Distance From Peak.*", "", motif_cols)
  clean_names <- gsub('[\r\n]', "", clean_names) 
  clean_names <- gsub('"', "", clean_names) 
  clean_names <- trimws(clean_names)
  colnames(motif_mat) <- clean_names

  #Count total peaks with hit per motif column
  peaks_with_motif <- colSums(motif_mat)
  hit_rate <- peaks_with_motif / total_peaks
  
  #Establish 5% random occurrence
  bg_rate <- 0.05
  
  #Create a 2X25 Format
  log2_enrichment <- log2((hit_rate + 0.001) / (bg_rate + 0.001))
  
  p_values <- sapply(seq_along(peaks_with_motif), function(i) {
    res <- binom.test(x = peaks_with_motif[i], n = total_peaks, p = bg_rate, alternative = "greater")
    return(res$p.value)
  })
  
  neglog10p <- -log10(p_values)
  neglog10p[is.na(neglog10p) | is.infinite(neglog10p)] <- 300 
  
  #Create matrix
  mat <- cbind(
    log2Enrichment = log2_enrichment,
    NegLog10P      = neglog10p
  )
  rownames(mat) <- colnames(motif_mat)
  
  #Sort rows by highest statistical significance
  mat <- mat[order(mat[, "NegLog10P"], decreasing = TRUE), , drop = FALSE]
  
  #Get top 25 motifs
  selected_n <- min(top_n, nrow(mat))
  mat <- mat[1:selected_n, , drop = FALSE]
  
  return(mat)
}

#Step 3: Create loop
mat <- process_subtype_to_model(files[subtype], subtype, top_n = 25)
  if (is.null(mat)) next
  
  #Scale columns and rows
  mat_scaled <- scale(mat)
  mat_scaled[is.na(mat_scaled)] <- 0 # Defend against NaNs if variance drops
  
  #Clustering
  hc <- hclust(dist(mat_scaled), method = "average")
  
  #Create pdfs
  output_filename <- paste0(subtype, "_HOMER_MotifHeatmap.pdf")
  pdf(output_filename, width = 8, height = 10)
  
  ht <- Heatmap(
    mat_scaled,
    name = "Z-score",
    cluster_rows = hc,
    cluster_columns = FALSE,
    
    col = colorRamp2(
      c(-2, 0, 2),
      c("blue", "white", "red")
    ),
    
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    column_title = paste(subtype, "HOMER Motifs"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold")
  )
  
  draw(ht)
  dev.off()
