#  This is the script for finding the Degs of the mice sample 
# install the packages

pak::pkg_install("mogene10sttranscriptcluster.db")

# Load the packages 

#  This is the script for finding the Degs of the mice sample 
# install the packages
pak::pkg_install("mogene10sttranscriptcluster.db")
# Load the packages 
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ggpubr)
library(EnhancedVolcano)
library(tidyverse)
library(dplyr)
library(DESeq2)
library(sva)
library(mogene10sttranscriptcluster.db)
# Load the counts data 
counts_17 <- read.table("GSE179717_series_matrix.txt.gz", 
                        header = TRUE, 
                        row.names = 1, 
                        sep = "\t", 
                        comment.char = "!", 
                        check.names = FALSE) # series matrix is manually downloaded from the ncbi geo
# Build the full metadata
meta_17 <- data.frame(
  sample    = colnames(counts_17),
  dataset   = "GSE179717",
  condition = c(
    rep("control", 7),        
    rep("STZ", 6),            
    rep("STZ_bezafibrate", 5)   
  )
)

# Subset for control and STZ only 

meta_17 <-meta_17[meta_17$condition%in%c("control","STZ"),]
counts_17 <- counts_17[,meta_17$sample]  
rownames(meta_17)<-meta_17$sample  

meta_17$condition <- factor(meta_17$condition, levels = c("control","STZ"))  


# Annotate the genes from the count data

gene_symbols <- mapIds(
  mogene10sttranscriptcluster.db,
  keys      = rownames(counts_17),
  column    = "SYMBOL",
  keytype   = "PROBEID",
  multiVals = "first"
)

annot_map <- data.frame(
  probe_id   = names(gene_symbols),
  mgi_symbol = unname(gene_symbols)
)
sum(!is.na(annot_map$mgi_symbol))

# Bind the annot_map and counts_17

counts_annot<- cbind(annot_map,counts_17)


# save the counts with gene symbols

write.csv(counts_annot, "GSE179717_raw_counts_symbo.csv", row.names = FALSE)
  
 
# Filterout the probes with NA symbols 

annot_clean<-annot_map |> 
  filter(!is.na(mgi_symbol)&mgi_symbol != "")

# Filter counts matrix to keep only annotated probes 

counts_annotated <- counts_17[annot_clean$probe_id,]

# Calculate the mean expression per probe to break duplicate gene ties 

probe_means <- rowMeans(counts_annotated)
  
# For duplicate gene symbols, keep the probe with the highest average signals

annot_unique <-annot_clean |> 
  mutate(mean_expr=probe_means[probe_id]) |> 
  group_by(mgi_symbol) |> 
  slice_max(order_by = mean_expr, n=1, with_ties = FALSE) |> 
  ungroup()

# Final filtered counts matrix

counts_filtered <- counts_annotated[annot_unique$probe_id,]
rownames(counts_filtered) <- annot_unique$mgi_symbol

# save clean annotated counts matrix 

counts_save <- cbind(annot_unique[,c("probe_id","mgi_symbol")], counts_filtered)
write.csv(counts_save,"GSE179717_filtered_count_symbols.csv", row.names = FALSE)

# Run limma on cleaned data 

design<- model.matrix(~0 + meta_17$condition)
colnames(design) <- levels(meta_17$condition)

# Fit model using filtrate matrix

fit <- lmFit(counts_filtered, design)

contrast_matrix <- makeContrasts(
  STZ_vs_control = STZ-control, levels = design)

fit2 <-contrasts.fit(fit, contrast_matrix)
fit2 <-eBayes(fit2)

# Extract DEGs table 
res_df <- topTable(fit2, number = Inf, adjust.method = "BH")

write.csv(res_df,"GSE179717_DEGs_1.csv")


res_df$Gene_Symbol <- rownames(res_df)


write.csv(res_df,"GSE179717_DEGs_limma_filtered_main.csv", row.names = FALSE)

## Making plots

volc <- EnhancedVolcano(
  res_df,
  lab = res_df$Gene_Symbol,
  x = 'logFC',
  y = 'adj.P.Val',                  # Standard for publications (FDR)
  pCutoff = 0.05,                   # FDR threshold
  FCcutoff = 1.0,                   # Log2 Fold Change threshold
  
  # Select top genes to label (prevents messy text overlap)
  selectLab = head(res_df[order(res_df$adj.P.Val), "Gene_Symbol"], 15),
  
  # Sizing & Transparency
  pointSize = 2.5,
  labSize = 4.0,
  labCol = 'black',
  labFace = 'bold',
  boxedLabels = TRUE,               # Adds clean white boxes around gene names
  drawConnectors = TRUE,            # Draws neat lines connecting label to dot
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  
  # Custom Publication Colors
  col = c('grey70', 'forestgreen', 'royalblue', 'firebrick3'),
  colAlpha = 0.8,
  
  # Titles & Labels
  title = 'GSE179717: STZ vs Control',
  subtitle = 'Differential Expression Analysis',
  caption = 'Cutoffs: |log2FC| > 1.0 & FDR < 0.05',
  xlab = bquote(~Log[2]~ 'Fold Change'),
  ylab = bquote(~-Log[10]~ 'Adjusted P-Value'),
  
  # Layout & Grid Clean-Up
  legendPosition = 'right',
  legendLabSize = 10,
  legendIconSize = 4.0,
  gridlines.major = FALSE,          # Removes distracting grid lines
  gridlines.minor = FALSE
) + 
  # Fix x-axis overlap by widening axis breaks and plot aspect ratio
  ggplot2::scale_x_continuous(breaks = seq(-6, 6, by = 2)) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5),
    axis.title = ggplot2::element_text(size = 12, face = "bold"),
    axis.text = ggplot2::element_text(size = 10, color = "black")
  )

# Display plot
print(volc)

ggsave("GSE179717_volcano_plot.png", plot = volc, height=8, width= 15,dpi = 1000)

# Heat maps

# select top 50 unique genes by adjusted p-value 
top_50_df <- res_df[order(res_df$adj.P.Val),][1:50,]
mat_top50 <- counts_filtered[top_50_df$Gene_Symbol,]

# Row-wise Z-score scaling 

mat_scaled <- t(scale(t(mat_top50)))

# Annotation metadata setup 

annotation_col <- data.frame(
  condition= meta_17$condition,
  row.names = meta_17$sample
)

annotation_colors <-list(
  condition= c(control="#1F77B4",STZ="#D62728")
)

# Plot clean heatmap 
 heat <- pheatmap(
   mat_scaled,
   cluster_rows = TRUE,
   cluster_cols = TRUE,
   show_rownames = TRUE,
   show_colnames = TRUE,
   annotation_col = annotation_col,
   annotation_colors = annotation_colors,
   color = colorRampPalette(c("navy","white","firebrick3"))(100),
   main = "Top 50 Differentially Expressed Genes of (Z-score) for GSE179717 Dataset"
 )

ggsave("GSE179717_heatmap.png", plot = heat, height=10, width= 15,dpi = 1000)










































