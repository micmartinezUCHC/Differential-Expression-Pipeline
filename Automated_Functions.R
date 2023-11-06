
#Define functions to be used in the automated DESeq2 script 
suppressWarnings({
  #Data manipulation
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(dplyr))
  
  #Differential gene expression
  suppressPackageStartupMessages(library(DESeq2))
  suppressPackageStartupMessages(library(ashr))
  
  #Graphics and visualizations
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(patchwork))
  
  #GSEA analysis
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(org.Rn.eg.db))
  suppressPackageStartupMessages(library(org.Mm.eg.db))
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(babelgene))
})

#-----Mapping ENSEMBL IDs to gene Symbols in an organism dependent manner
#Function to map gene Symbols to ENSEMBL IDs for either mouse or rat
getSymbols <- function(RAWCOUNTS, ORGANISM) {
  if (ORGANISM == "rno") {
    RAWCOUNTS$Symbol <- mapIds(org.Rn.eg.db, key = RAWCOUNTS$Gene, column ="SYMBOL",
                               keytype = "ENSEMBL", multiVals = "first")
  } else {
    RAWCOUNTS$Symbol <- mapIds(org.Mm.eg.db, key = RAWCOUNTS$Gene, column = "SYMBOL",
                               keytype ="ENSEMBL", multiVals = "first")
  }
  raw <- RAWCOUNTS
  return(raw)
}


#-----Filtering and Preparation
preFilt <- function(RAWCOUNTS) {
  #Set variable
  raw <- RAWCOUNTS
  
  #Omit any gene with no symbol annotation or unknown genes
  raw <- raw[!is.na(raw$Symbol),]
  raw <- raw[!grepl("^LOC\\d+$", raw$Symbol),]
  raw <- raw[!grepl("^RGD\\d+$", raw$Symbol),]
  print("Data successfully loaded")
  print(paste("Number of columns in raw data = ", ncol(raw), sep = " "))
  print(colnames(raw))
  
  #Format the rownames to join Ensembl ID and symbol
  raw$geneIDs <- paste(raw$Gene, raw$Symbol, sep = " - ")
  #print(nrow(raw))
  rownames(raw) <- raw$geneIDs
  raw$Gene <- NULL
  raw$X <- NULL
  raw$Symbol <- NULL
  raw$geneIDs <- NULL
  print("Gene symbols successfully acquired")
  
  #Filter
  minimumCountpergene <- 10
  MinSampleWithminimumgeneCounts <- 5
  
  #Filter out low read counts for Normal vs control, Naproxen vs control, etc...
  raw <- raw[rowSums(data.frame(raw>minimumCountpergene)) > MinSampleWithminimumgeneCounts,]
  print("Filtering completed")
  
  raw <- return(raw)
}


#-----Generating Summary heatmap from the DEGs
degHeatmap <- function(TOPDEGS) {
  
  #Define variable
  top.degs <- TOPDEGS
  
  #Check if there are 100 or more rows. If true, results will be a glimpse at the top 50 and bottom 50
  if (nrow(top.degs) > 100) {
    #---Specify output tiff file and resolution
    tiff(paste("Figures", paste(ref, "vs", treatment, "heatmap.tiff", sep = "_"), sep = "/"), width = 10, height = 8, pointsize = 16, units = "in", res=1400)
    print("More than 100 DEGs found: Taking top 50 up-regulated and top 50 down-regulated")
    
    #Take the top 50 and bottom 50
    total_rows <- nrow(top.degs)
    top.degs.keep <- c(1:50, (total_rows - 49):total_rows)
    top.degs.subset <- top.degs[top.degs.keep,]
    
    #Pull the baseMean column and Log2FC columns as lists
    top.degs.log2fc <- as.matrix(top.degs.subset$log2FoldChange)
    colnames(top.degs.log2fc) <- "Log2 FC"
    top.degs.baseMean <- as.matrix(log(top.degs.subset$baseMean))
    colnames(top.degs.baseMean) <- "Log Base Mean"
    
    #Extract just the first column and the counts columns
    top.degs.subset <- top.degs.subset[,c(1,8:ncol(top.degs.subset))]
    rownames(top.degs.subset) <- top.degs.subset$Row.names
    top.degs.subset$Row.names <- NULL
    
    #Transpose, center, and scale the normalized counts
    scaled <- t(apply(top.degs.subset, 1, scale))
    colnames(scaled) <- colnames(top.degs.subset)
    
    #Map colors to values
    l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
                                0,
                                max(top.degs.log2fc)),
                              c("blue", "white", "red"))
    
    mean.colors <- colorRamp2(c(quantile(top.degs.baseMean)[1],
                                quantile(top.degs.baseMean)[4]),
                              c("white", "red"))
    
    #Isolate gene symbols for labeling
    labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(top.degs.subset))
    top.degs.subset$Symbols <- labels
    
    
    # Create a named vector for the colors
    color_vector <- c("tan", "black")
    names(color_vector) <- c(ref, treatment)
    
    #Set heatmap splitting pattern
    hmSplit <- rep(1:2, c(n, m))
    
    #Define the number of slices in the heatmap
    slices <- n+m
    
    #Create a heatmap annotation
    hmAnno <- HeatmapAnnotation(
      empty = anno_empty(border = FALSE),
      Group = anno_block(gp = gpar(fill = 2:slices), labels = c(ref, treatment)
      ))
    
    #Optional side annotation
    ha <- rowAnnotation(logBaseMean = anno_barplot(top.degs.baseMean, height = unit(4, "cm"), 
                                                   gp = gpar(fontsize = 7)),
                        annotation_name_rot = 90,
                        annotation_name_gp = gpar(fontsize = 10))
    
    #Heatmap for scaled data
    hmScaled <- Heatmap(scaled,
                        column_labels = colnames(scaled), 
                        name = "Z-score",
                        # top_annotation = HeatmapAnnotation(group = anno_block(gp = gpar(1:2, c(n, m))),
                        #                                    labels = rep(c(ref, treatment), times = c(n, m))),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        top_annotation = hmAnno,
                        #right_annotation = ha,
                        column_split = hmSplit,
                        row_names_gp = gpar(fontsize = 5),
                        column_title = paste(ref, "vs", treatment, sep = " "))
    #Heatmap for log2FC values
    hml2FC <- Heatmap(top.degs.log2fc,
                      row_labels = labels,
                      row_names_gp = gpar(fontsize = 5),
                      cluster_rows = FALSE,
                      name = "log2FC",
                      col = l2FC.colors,
                      cell_fun = function(j,i, x, y, w, h, col) {
                        grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                  gp = gpar(fontsize = 5, 
                                            col = "black"))})
    #Heatmap for average expression
    hmMean <- Heatmap(top.degs.baseMean,
                      row_labels = labels,
                      row_names_gp = gpar(fontsize = 5),
                      cluster_rows = FALSE,
                      name = "Log Base Mean",
                      col = mean.colors,
                      cell_fun = function(j, i, x, y, w, h, col) {
                        grid.text(round(top.degs.baseMean[i, j],1), x, y,
                                  gp = gpar(fontsize = 5,
                                            col = "black"))})
    #Draw the final heatmap
    HM <- hmScaled + hml2FC + hmMean
    
    #Draw the complex heatmap
    draw(HM)
    #Close the PDF device
    dev.off()
    
  } else {
    print("There is less than 100 DEGs present! Using all DEGs to generate heatmap")
    tiff(paste("Figures", paste(ref, "vs", treatment, "heatmap.tiff", sep = "_"), sep = "/"), width = 10, height = 8, pointsize = 16, units = "in", res=1400)
    
    #Redefine top degs
    top.degs.subset <- top.degs[order(top.degs$log2FoldChange, decreasing = TRUE),]
    
    #Pull the baseMean column and Log2FC columns as lists
    top.degs.log2fc <- as.matrix(top.degs.subset$log2FoldChange)
    colnames(top.degs.log2fc) <- "Log2FC"
    top.degs.baseMean <- as.matrix(log(top.degs.subset$baseMean))
    colnames(top.degs.baseMean) <- "BaseMean"
    
    #Extract just the first column and the counts columns
    top.degs.subset <- top.degs.subset[,c(1,8:ncol(top.degs.subset))]
    rownames(top.degs.subset) <- top.degs.subset$Row.names
    top.degs.subset$Row.names <- NULL
    
    #Transpose, center, and scale the normalized counts
    scaled <- t(apply(top.degs.subset, 1, scale))
    colnames(scaled) <- colnames(top.degs.subset)
    
    #Map colors to values
    l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
                                0,
                                max(top.degs.log2fc)),
                              c("blue", "white", "red"))
    
    mean.colors <- colorRamp2(c(quantile(top.degs.baseMean)[1],
                                quantile(top.degs.baseMean)[4]),
                              c("white", "red"))
    
    # l2FC.colors <- colorRamp2(c(min(top.degs.log2fc),
    #                             0,
    #                             max(top.degs.log2fc)),
    #                           c("blue", "white", "red"))
    
    #Isolate gene symbols for labeling
    labels <- gsub("^[^-]+-(.*)$", "\\1", rownames(top.degs.subset))
    top.degs.subset$Symbols <- labels
    
    
    # Create a named vector for the colors
    color_vector <- c("tan", "black")
    names(color_vector) <- c(ref, treatment)
    
    
    #Annotation column
    # hmAnno <- HeatmapAnnotation(group = design$Group,
    #                                name = "", 
    #                                show_annotation_name = FALSE,
    #                                col = list(group = color_vector))
    
    
    
    #Set heatmap splitting pattern
    hmSplit <- rep(1:2, c(n, m))
    
    #Define the number of slices in the heatmap
    slices <- n+m
    
    #Create a heatmap annotation
    hmAnno <- HeatmapAnnotation(
      empty = anno_empty(border = FALSE),
      Group = anno_block(gp = gpar(fill = 2:slices), labels = c(ref, treatment)
      ))
    
    ha <- rowAnnotation(logBaseMean = anno_barplot(top.degs.baseMean, height = unit(4, "cm"), 
                                                   gp = gpar(fontsize = 7)),
                        annotation_name_rot = 90,
                        annotation_name_gp = gpar(fontsize = 10))
    
    
    #Heatmap for scaled data
    hmScaled <- Heatmap(scaled,
                        column_labels = colnames(scaled), 
                        name = "Z-score",
                        # top_annotation = HeatmapAnnotation(group = anno_block(gp = gpar(1:2, c(n, m))),
                        #                                    labels = rep(c(ref, treatment), times = c(n, m))),
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        top_annotation = hmAnno,
                        right_annotation = ha,
                        column_split = hmSplit,
                        row_names_gp = gpar(fontsize = 5),
                        column_title = paste(ref, "vs", treatment, sep = " "))
    #Heatmap for log2FC values
    hml2FC <- Heatmap(top.degs.log2fc,
                      row_labels = labels,
                      row_names_gp = gpar(fontsize = 6),
                      cluster_rows = FALSE,
                      name = "log2FC",
                      col = l2FC.colors,
                      cell_fun = function(j,i, x, y, w, h, col) {
                        grid.text(round(top.degs.log2fc[i, j],1), x, y, 
                                  gp = gpar(fontsize = 5, 
                                            col = "black"))})
    #Heatmap for average expression
    # hmMean <- Heatmap(top.degs.baseMean,
    #                    row_labels = labels,
    #                    row_names_gp = gpar(fontsize = 5),
    #                    cluster_rows = FALSE,
    #                    name = "Log Base Mean",
    #                    col = mean.colors,
    #                    cell_fun = function(j, i, x, y, w, h, col) {
    #                      grid.text(round(top.degs.baseMean[i, j],1), x, y,
    #                                gp = gpar(fontsize = 5,
    #                                          col = "black"))})
    #Draw the final heatmap
    HM <- hmScaled + hml2FC #+ hmMean
    
    #Draw the complex heatmap
    draw(HM)
    #Close the PDF device
    dev.off()
    
  }
  heatmap <- return(HM)
}








