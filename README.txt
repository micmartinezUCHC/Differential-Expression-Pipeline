Mike Martinez: Differential Expression Pipeline (DEPl)
Rosenberg Lab, UConn Health Center (UCHC)


DEPl is a semi-automated RNASeq analysis pipeline designed for simple differential expression experiments containing one factor designs.
DEPl runs differential expression analysis and geneset enrichment analyses using the popular Bioconductor packages DESeq2 (Love, Huber, Anders 2014) and clusterProfiler (Yu, Wang, Han, He 2012). 

This pipeline has some pre-requisites
1.) You must create a folder named "Pipeline_Results" in a location where you want your output files to be directed. You can manually change this base output directory name
	directly from the DEPl.sh script on line 26. Pipeline results are stored in date/time stamped folders.
2.) You must have R installed on your computer and install the required packages (see the RScripts for details)

This pipeline outputs the following files
1.) A folder called "DESeq2_Results"
	a.) A csv file containing all DEGs
	b.) a csv file containing DEGs filtered by padj < 0.05 and abs(log2FC) > 2
2.) A .rds file containing the DESeq2 object
3.) A folder called "Figures"
	a.) a heatmap showing the top 50 and bottom 50 differentially expressed genes (if there are more than 100 DEGs passing padj and fold change thresholds), or all DEGs if number of DEGs < 100.
	b.) a PCA of the firt 2 principal components
4.) A folder called "GSEA_Results"
	a.) Gene Ontology enrichment tabular results
	b.) Gene Ontology dotplot
	c.) KEGG enrichment tabular results
	d.) KEGG enrichment dotplot
	e.) Gene Ontology gseGO .rds object
	f.) KEGG gseKEGG .rds object
5.) Summary.txt file containing gene expression information summary and runtime summary

The script requires user-input to run:
1.) Run mode:
	1: Only run DESeq2
	2: Run DESeq2 and GSEA
2.) Path to counts file
	Must specify the absolute path. Do not use "". Ex: /path/to/counts.csv
	Counts file must have the first column names "Gene" and must be Ensembl IDs
3.) Name of Analysis
	This will be used as the basename for all output files. Do not include spaces or hyphens. Use CamelCase, snakeCase or underscores
4.) Reference level
	Used as the denominator level for analysis
5.) Number of observations of reference level 
6.) Numerator level for analysis
7.) Number of observations of numerator level
8.) Organism
	Either rno or mmu (for rat and mice)
	DO NOT capitalize or use quotes. 
	

	
	



