#Version 2 of the automated GSEA script. Now includes functionality to run on either rar or mouse samples
#Load libraries

#record script run-time start
start_time <- Sys.time()

suppressWarnings({
  #Data manipulation
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(argparse))
  
  #Graphics and visualizations
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ggh4x))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(patchwork))
  
  #GSEA analysis
  suppressPackageStartupMessages(library(clusterProfiler))
  suppressPackageStartupMessages(library(AnnotationDbi))
  suppressPackageStartupMessages(library(msigdbr))
  suppressPackageStartupMessages(library(enrichplot))
  suppressPackageStartupMessages(library(babelgene))
})

#Add arguments
parser <- ArgumentParser()
parser$add_argument("--DEGs", type = "character", default = NULL, help = "Path to DEG table with counts")
parser$add_argument("--name", type = "character", default = NULL, help = "Name of Analysis") #This should be kept the same from the first script
parser$add_argument("--ref", type = "character", default = NULL, help = "Reference group for DESeq2")
parser$add_argument("--treatment", type = "character", default = NULL, help = "Name of second group")
parser$add_argument("--out", type = "character", default = NULL, help = "Path to GSEA output folder") #Note: this should be ONE level within the directory
parser$add_argument("--org", type = "character", default = NULL)

#Parse arguments
parsed_args <- parser$parse_args()

#Define argument names to use within the script
DEGs <- parsed_args$DEGs
name <- parsed_args$name
ref <- parsed_args$ref
treatment <- parsed_args$treatment
out <- parsed_args$out
org <- parsed_args$org
DB <- parsed_args$DB
#ortho <- parsed_args$ortho

#Setting up output directory
output_dir <- out

#Check if output directory already exists
if (!file.exists(output_dir)) {
  dir.create(output_dir)
  cat("Output directory successfully created. \n")
} else {
  cat("Output directory already exists. \n")
}

#Set working directory
setwd(output_dir)

#Associate Organism Code to organism
if (org == "rno") {
  species <- "Rattus norvegicus"
  library(org.Rn.eg.db)
  DB = "org.Rn.eg.db"

} else {
  species <- "Mus musculus"
  library(org.Mm.eg.db)
  DB = "org.Mm.eg.db"
}

print("#----- Mapping Ensembl to Entrez IDs -----#")
print(species)
print(DB)


#Read in the DEG table
DEG.results <- read.csv(DEGs, header = TRUE, sep = ",")

#Map ENTREZ IDs to the ENSEMBL IDs
DEG.results$Entrez <- mapIds(get(DB), key = DEG.results$Ensembl,
                             column = "ENTREZID", keytype = "ENSEMBL",
                             multiVals = "first")

#Get a list of the genes
res.ordered <- DEG.results[order(DEG.results$log2FoldChange, decreasing = TRUE),]
res.ordered.genes <- res.ordered$log2FoldChange

#Assign Entrez IDs as names for the genes
names(res.ordered.genes) <- res.ordered$Entrez

#Remove duplicated Entrez IDs and their corresponding values
unique_entrez_genes <- names(res.ordered.genes[!duplicated(names(res.ordered.genes))])
unique_genes <- res.ordered.genes[unique_entrez_genes]
unique_genes <- sort(unique_genes, decreasing = TRUE)

print("#----- Running gseGO -----#")
#Run gseGO
gse <- gseGO(unique_genes,
             ont = "all",
             OrgDb = get(DB),
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             eps = 1e-300,
             verbose = TRUE,
             by = "fgsea")
saveRDS(gse, file = paste(name, "gseGO.rds", sep = "_"))

if (length(gse@result$Description > 0)) {
  gse.readable <- as.data.frame(setReadable(gse, OrgDb = get(DB), keyType = "ENTREZID"))
  write.csv(gse.readable, file = paste(name, "GO_enrichment_results.csv", sep = "_"))
  
  # Custom facet label mapping
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = paste("Enriched in", treatment, sep = " "), suppressed = paste("Enriched in", ref, sep = " "))
  )
  
  # Create the dotplot with custom facet labels
  GO.dotplot <- dotplot(gse, 
                        showCategory = 15, 
                        split = ".sign", 
                        font.size = 5.0, 
                        label_format = 50,
                        title = paste(name, "GSEA: GO", sep = " ")) + 
    facet_nested(~.sign + ONTOLOGY~., labeller = custom_labels, scales = "free")
  ggsave(paste(name, "GO_enrichment_dotplot.pdf", sep = "_"), GO.dotplot, width = 12, height = 10)
} else {
  print("No GO terms enriched")
}

print("#----- Running gseKEGG -----#")
kegg <- gseKEGG(unique_genes,
                organism = org,
                keyType = "kegg",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH",
                eps = 1e-300,
                verbose = TRUE,
                by = "fgsea")
saveRDS(kegg, file = paste(name, "gseKEGG.rds", sep = "_"))

if (length(kegg@result$Description > 0)) {
  kegg.readable <- as.data.frame(setReadable(kegg, OrgDb = get(DB), keyType = "ENTREZID"))
  write.csv(kegg.readable, file = paste(name, "KEGG_enrichment_results.csv", sep = "_"))
  
  # Custom facet label mapping
  custom_labels <- labeller(
    .default = label_value,
    .sign = c(activated = paste("Enriched in", treatment, sep = " "), suppressed = paste("Enriched in", ref, sep = " "))
  )
  
  # Create the dotplot with custom facet labels
  kegg.dotplot <- dotplot(kegg, 
                          showCategory = 15, 
                          split = ".sign", 
                          font.size = 7.5, 
                          label_format = 50,
                          title = paste(name, "GSEA: KEGG", sep = " ")) + 
    facet_grid(~.sign, labeller = custom_labels)
  ggsave(paste(name, "KEGG_enrichment_dotplot.pdf", sep = "_"), kegg.dotplot, width = 12, height = 10)
}


#------End
#Record end time
end_time <- Sys.time()

#Calculate elapsed time
elapsed_time <- end_time - start_time
formmated_time <- sprintf("%.2f", elapsed_time)

cat("GSEA completed in:", formmated_time, "Seconds", sep = " ")
cat("\n")





