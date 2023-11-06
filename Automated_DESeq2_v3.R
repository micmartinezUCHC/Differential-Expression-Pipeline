#The purpose of this script is to break the DEG pipeline into smaller easier to run code sections

######### USERS GUIDE #########
#Options#
# --counts <path/to/counts/file>
# --name <base name for all output files>
# --ref <name of reference group>
# --n <number of observations in reference group>
# --treatment <name of treatment group>
# --m <number of observations in treatment group>
# --out <path/to/output/directory>

#####REQUIRED FILE: You need the file "Cell_marker_Human.csv"
#Available from this link: http://117.50.127.228/CellMarker/CellMarker_download.html

#---------------# Parsing Arguments
#Set command line arguments
args <- commandArgs(trailingOnly = TRUE)

#record script run-time start
start_time <- Sys.time()

#Initialize an empty list to store names of output files and figures
output_files <- list()
output_figures <- list()

suppressPackageStartupMessages(library(argparse))

#Add arguments
parser <- ArgumentParser()
parser$add_argument("--counts", type = "character", default = NULL, help = "Path to raw counts matrix")
parser$add_argument("--name", type = "character", default = NULL, help = "Basename for all output files")
parser$add_argument("--ref", type = "character", default = NULL, help = "Reference group for DESeq2")
parser$add_argument("--n", type = "integer", default = NULL, help = "Number of observations for reference group")
parser$add_argument("--treatment", type = "character", default = NULL, help = "Name of second group")
parser$add_argument("--m", type = "integer", default = NULL, help = "Number of observations for second group")
parser$add_argument("--out", type = "character", default = NULL, help = "Output directory path")
parser$add_argument("--org", type = "character", default = NULL, help = "Organism code: either mmu for Mus musculus or rno for Rattus norvegicus")
#parser$add_argument("--ortho", type = "character", default = NULL, help = "Path to human cell-type data")

#Parse arguments
parsed_args <- parser$parse_args()

#Define argument names to use within the script
counts <- parsed_args$counts
name <- parsed_args$name
ref <- parsed_args$ref
n <- parsed_args$n
treatment <- parsed_args$treatment
m <- parsed_args$m
out <- parsed_args$out
org <- parsed_args$org
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
dir.create("DESeq2_Results")
dir.create("Figures")

#---------------# Beginning

#Load libraries
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

# --------------- Pre-processing -------------------#

#-----SOURCE FUNCTIONS-----#
source("/Users/mikemartinez/Desktop/DEG_Pipeline/Automated_Functions.R")

print("#----- Initial Count Processing -----#")
#Read in the data
raw <- read.csv(counts, header = TRUE, sep = ",")

#Obtain gene symbols (automated function: see Automated_Functions.R)
#---Filter out unknowns, filter out low read counts
raw <- preFilt(getSymbols(raw, org))
system("clear")


# --------------- DESeq2 -------------------#

print("#----- Creating Design Matrix -----#")
print(paste("Reference group:", ref, sep = " "))
print(paste("2nd Group:", treatment, sep = " "))

#Un-comment if something weird happens
#print(colnames(raw))

design <- data.frame(Sample = rep(c(ref, treatment),
                                  c(n,m)),
                     Group = rep(c(ref, treatment),
                                 c(n,m)))
rownames(design) <- colnames(raw)
print(design)

#Sanity check: are all colnames in raw present in the design table?
check1 <- all(colnames(raw) %in% rownames(design))

if (check1 != "TRUE") {
  stop("Error: Not all samples in counts file are present in the design table.")
}

check2 <- all(colnames(raw) == rownames(design))

if (check2 != TRUE) {
  stop("Error: Samples are not in the same order between counts file and design matrix")
}


print("#----- Running Deseq2 Algorithm -----#")
dds <- DESeqDataSetFromMatrix(countData = raw,
                              colData = design,
                              design = ~ Group)
  
#Set a factor level
dds$Group <- relevel(dds$Group, ref = ref)

#Set random seed
set.seed(03061999)

#Run DESeq2 algorithm
dds <- DESeq(dds)

#Save dds object
saveRDS(dds, paste(name, "DifferentialExpression.RDS", sep = "_"))
print("DESeq2 successful !!!!!")

#View the results
res <- results(dds)
res.df <- as.data.frame(res)

#Get the normalixed counts
counts <- counts(dds, normalized = TRUE)

res.ordered <- res.df[order(res.df$log2FoldChange, decreasing = TRUE),]
res.ordered <- merge(res.ordered, counts, by = 0, all = TRUE)
res.ordered$Ensembl <- gsub("^(.*?) - .*", "\\1", res.ordered$Row.names)
res.ordered$Symbols <- gsub("^[^-]+-(.*)$", "\\1", res.ordered$Row.names)
write.csv(res.ordered, file = paste("DESeq2_Results", paste(name, "All_DEGs.csv", sep = "_"), sep = "/"))


# --------------- PCA -------------------#

print("#----- Performing Variance Stabilizing Transformating and plotting first 2 principal components -----#")

#Create a summarized experiment
vsd <- vst(dds)

#Plot PCA
PCA <- plotPCA(vsd, intgroup = "Group") +
  geom_text_repel(aes(label = rownames(design))) +
  ggtitle(paste(ref, "vs", treatment, sep = " ")) +
  stat_ellipse(geom = "polygon", type = "norm", level = 0.90, alpha = 0.05, aes(fill = group)) +
  theme_bw()
ggsave(paste("Figures", paste(name, "PCA.tiff", sep = "_"), sep = "/"), dpi = 800)


system("clear")


print("#----- Generating Heatmap -----#")

# --------------- Heatmap Generation -------------------#

#Filter the data frame based on padj and log2FC
res.top <- res.ordered[res.ordered$padj < 0.05 & abs(res.ordered$log2FoldChange) > 2 & res.ordered$baseMean > 50,]
res.mid <- res.ordered[res.ordered$padj < 0.05 & abs(res.ordered$log2FoldChange) >= 1,]
res.top <- na.omit(res.top)
res.mid <- na.omit(res.mid)

#Order the genes by log2FC
res.top <- res.top[order(res.top$log2FoldChange, decreasing = TRUE),]
write.csv(res.top, file = paste("DESeq2_Results", paste(name, "filteredDEGs_padj005_log2FC2.csv", sep = "_"), sep = "/"))


#Heatmap fo top 50 and bottom 50 
#Order top.degs df by decreasing log2FC
top.degs <- res.top[order(res.top$log2FoldChange, decreasing = TRUE),]
top.degs$Ensembl <- NULL
top.degs$Symbols <- NULL

HM <- degHeatmap(top.degs)

# Define the file path for the summary file
summary_file <- "Summary.txt"

# Define the information to be appended to summary file
Upreg <- res.top[res.top$log2FoldChange >= 2,]
Dnreg <- res.top[res.top$log2FoldChange <= -2,]
mid_Upreg <- res.mid[res.mid$log2FoldChange >= 1,]
mid_Dnreg <- res.mid[res.mid$log2FoldChange <= -1,]


#------End
#Record end time
end_time <- Sys.time()

#Calculate elapsed time
elapsed_time <- end_time - start_time
formmated_time <- sprintf("%.2f", elapsed_time)

# Print the script run-time
cat("------------------------------\n")
cat("Differential Expression Analysis Completed in:", formmated_time, "Seconds", sep = " ")
cat("\n")

#Append information to summary file
info <- c(
  paste("Total number of Genes(unfiltered):", nrow(res.ordered)),
  paste("Total number of filtered DEGs p < 0.05, l2FC > 2:", nrow(res.top)),
  paste("Total number of filtered up-regulated DEGs (l2FC > 2):", nrow(Upreg)),
  paste("Total number of filtered down-regulated DEGs (l2FC > 2):", nrow(Dnreg)),
  paste("Total number of filtered DEGs p < 0.05, l2FC > 1:", nrow(res.mid)),
  paste("Total number of filtered up-regulated DEGs (l2FC > 1):", nrow(mid_Upreg)),
  paste("Total number of filtered down-regulated DEGs (l2FC > 1):", nrow(mid_Dnreg)),
  paste("Differential Expression Analysis Completed in:", formmated_time, "Seconds", sep = " ")
)

# Write the information to the summary file
cat(info, file = summary_file, sep = "\n")



