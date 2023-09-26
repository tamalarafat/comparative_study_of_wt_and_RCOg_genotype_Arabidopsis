# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_marker_identification", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/Functions_matrix_manipulation", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Col0.leaf.vs.protoplast_table_2pseudorep_final_August_2022.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# WT A. thaliana
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
COL_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_1ST_2/filtered_feature_bc_matrix/")

# Load data - WT OX 2nd Experiment (leaf 6 and 7)
COL_data_2E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col_RNA_2nd_ALL_2/filtered_feature_bc_matrix/")

# Load data - WT OX 3rd Experiment (leaf 5 and 6)
COL_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_3rd_ALL/filtered_feature_bc_matrix/")

# Load data - WT OX 7th Experiment (leaf 6 and 7)
COL_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_5th_ALL_2/filtered_feature_bc_matrix/")


###
# RCOg A. thaliana
###

# Data from 3rd experiment
RCO_data_3E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_3rd_ALL/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_5E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_5th_ALL_2/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_6E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_6th_ALL/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_7E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_7th_ALL_2/filtered_feature_bc_matrix/")

# Data from 3rd experiment
RCO_data_8E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_8th_New/filtered_feature_bc_matrix/", gene.column = 1)

# Data from 3rd experiment
RCO_data_10E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_RCOg_RNA_10th_New/filtered_feature_bc_matrix/", gene.column = 1)


# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(thaliana_genes, PP_genes)

##### Remove the protoplasting induced genes
COL_data_1E <- COL_data_1E[genes_to_keep, ]
COL_data_2E <- COL_data_2E[genes_to_keep, ]
COL_data_3E <- COL_data_3E[genes_to_keep, ]
COL_data_5E <- COL_data_5E[genes_to_keep, ]

RCO_data_3E <- RCO_data_3E[genes_to_keep, ]
RCO_data_5E <- RCO_data_5E[genes_to_keep, ]
RCO_data_6E <- RCO_data_6E[genes_to_keep, ]
RCO_data_7E <- RCO_data_7E[genes_to_keep, ]
RCO_data_8E <- RCO_data_8E[genes_to_keep, ]
RCO_data_10E <- RCO_data_10E[genes_to_keep, ]


# Create seurat object and perform initial filtering - 
# 1. Remove genes if their expression was not detected in at least one cell out of every 500 ("min.cells"),
# 2. Remove cells if at least 200 genes were not detected to be expressed (min.features = 200),
# 3. Remove cells with a total count of more than 110000 (nCount_RNA > 110000).
# 4. Remove cells if 5% or more of the total count of a cell belongs to the mitochondiral genes.
# 5. Remove cells if 10% or more of the total count of a cell belongs to the chloroplast genes.

###
# COL - 1 E
###

# First replicate - COL 1E - total cells 4850;
COL_1E <- CreateSeuratObject(counts = COL_data_1E, project = "COL_1E", min.features = 200)

# Add metadata information to the seurat object
COL_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_1E <- subset(COL_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_1E[["percent.mt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_1E[["percent.pt"]] <- PercentageFeatureSet(COL_1E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_1E <- subset(COL_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_1E <- NormalizeData(COL_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_1E <- FindVariableFeatures(COL_1E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W1 <- GetAssayData(COL_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W1) <- paste("C1", colnames(COL_W1), sep = "_")

###
# COL - 2 E
###

# First replicate - COL 2E - total cells 10760; 
COL_2E <- CreateSeuratObject(counts = COL_data_2E, project = "COL_2E", min.features = 200)

# Add metadata information to the seurat object
COL_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_2E <- subset(COL_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_2E[["percent.mt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_2E[["percent.pt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_2E <- subset(COL_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_2E <- NormalizeData(COL_2E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_2E <- FindVariableFeatures(COL_2E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W2 <- GetAssayData(COL_2E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W2) <- paste("C2", colnames(COL_W2), sep = "_")

###
# COL - 3 E
###

# First replicate - COL 3E - total cells 4100;
COL_3E <- CreateSeuratObject(counts = COL_data_3E, project = "COL_3E", min.features = 200)

# Add metadata information to the seurat object
COL_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_3E <- subset(COL_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_3E[["percent.mt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_3E[["percent.pt"]] <- PercentageFeatureSet(COL_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_3E <- subset(COL_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_3E <- NormalizeData(COL_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_3E <- FindVariableFeatures(COL_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W3 <- GetAssayData(COL_3E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W3) <- paste("C3", colnames(COL_W3), sep = "_")

###
# COL - 5 E
###

# First replicate - COL 5E - total cells 8420;
COL_5E <- CreateSeuratObject(counts = COL_data_5E, project = "COL_5E", min.features = 200)

# Add metadata information to the seurat object
COL_5E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-5", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_5E <- subset(COL_5E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_5E[["percent.mt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_5E[["percent.pt"]] <- PercentageFeatureSet(COL_5E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_5E <- subset(COL_5E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
COL_5E <- NormalizeData(COL_5E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
COL_5E <- FindVariableFeatures(COL_5E, selection.method = "vst", nfeatures = 3000)

# Extract the count table from the seurat object
COL_W5 <- GetAssayData(COL_5E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W5) <- paste("C5", colnames(COL_W5), sep = "_")


###
# Intersection of highly variable genes among replicates - WT Arabidopsis thaliana
###

# Lets find common variable features for the replicates of A. thaliana
wt_hvgs = Reduce(intersect, list(COL_1E@assays$RNA@var.features, COL_2E@assays$RNA@var.features, COL_3E@assays$RNA@var.features, COL_5E@assays$RNA@var.features))

fileGenerator(wt_hvgs, fileName = "Shared_highly_variable_genes_between_replicates_WT_AT.txt")


###
# RCOg - 3 E
###

# Total cells 4550;
RCO_3E <- CreateSeuratObject(counts = RCO_data_3E, project = "RCO_3E", min.features = 200)

# Add metadata information to the seurat object
RCO_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-3", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_3E <- subset(RCO_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_3E[["percent.mt"]] <- PercentageFeatureSet(RCO_3E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_3E[["percent.pt"]] <- PercentageFeatureSet(RCO_3E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_3E <- subset(RCO_3E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_3E <- NormalizeData(RCO_3E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_3E <- FindVariableFeatures(RCO_3E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G3 <- GetAssayData(RCO_3E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G3) <- paste("R3", colnames(RCO_G3), sep = "_")


###
# RCOg - 5 E
###

# Total cells 7000; 
RCO_5E <- CreateSeuratObject(counts = RCO_data_5E, project = "RCO_5E", min.features = 200)

# Add metadata information to the seurat object
RCO_5E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-5", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_5E <- subset(RCO_5E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_5E[["percent.mt"]] <- PercentageFeatureSet(RCO_5E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_5E[["percent.pt"]] <- PercentageFeatureSet(RCO_5E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_5E <- subset(RCO_5E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_5E <- NormalizeData(RCO_5E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_5E <- FindVariableFeatures(RCO_5E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G5 <- GetAssayData(RCO_5E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G5) <- paste("R5", colnames(RCO_G5), sep = "_")


###
# RCOg - 6 E
###

# Total cells 4550; 
RCO_6E <- CreateSeuratObject(counts = RCO_data_6E, project = "RCO_6E", min.features = 200)

# Add metadata information to the seurat object
RCO_6E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-6", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_6E <- subset(RCO_6E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_6E[["percent.mt"]] <- PercentageFeatureSet(RCO_6E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_6E[["percent.pt"]] <- PercentageFeatureSet(RCO_6E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_6E <- subset(RCO_6E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_6E <- NormalizeData(RCO_6E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_6E <- FindVariableFeatures(RCO_6E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G6 <- GetAssayData(RCO_6E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G6) <- paste("R6", colnames(RCO_G6), sep = "_")


###
# RCOg - 7 E
###

# Total cells 4550;
RCO_7E <- CreateSeuratObject(counts = RCO_data_7E, project = "RCO_7E", min.features = 200)

# Add metadata information to the seurat object
RCO_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-7", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_7E <- subset(RCO_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_7E[["percent.mt"]] <- PercentageFeatureSet(RCO_7E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_7E[["percent.pt"]] <- PercentageFeatureSet(RCO_7E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_7E <- subset(RCO_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_7E <- NormalizeData(RCO_7E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_7E <- FindVariableFeatures(RCO_7E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G7 <- GetAssayData(RCO_7E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G7) <- paste("R7", colnames(RCO_G7), sep = "_")


###
# RCOg - 8 E
###

# Total cells 4550;
RCO_8E <- CreateSeuratObject(counts = RCO_data_8E, project = "RCO_8E", min.features = 200)

# Add metadata information to the seurat object
RCO_8E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-8", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_8E <- subset(RCO_8E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_8E[["percent.mt"]] <- PercentageFeatureSet(RCO_8E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_8E[["percent.pt"]] <- PercentageFeatureSet(RCO_8E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_8E <- subset(RCO_8E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_8E <- NormalizeData(RCO_8E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_8E <- FindVariableFeatures(RCO_8E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G8 <- GetAssayData(RCO_8E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G8) <- paste("R8", colnames(RCO_G8), sep = "_")


###
# RCOg - 10 E
###

# Total cells 1000;
RCO_10E <- CreateSeuratObject(counts = RCO_data_10E, project = "RCO_10E", min.features = 200)

# Add metadata information to the seurat object
RCO_10E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "RCOg-10", "RCO", "Leaf")

# Remove cells with a total count more than 110000
RCO_10E <- subset(RCO_10E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
RCO_10E[["percent.mt"]] <- PercentageFeatureSet(RCO_10E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
RCO_10E[["percent.pt"]] <- PercentageFeatureSet(RCO_10E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
RCO_10E <- subset(RCO_10E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
RCO_10E <- NormalizeData(RCO_10E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
RCO_10E <- FindVariableFeatures(RCO_10E, selection.method = "vst", nfeatures = 3000)

# Extract the count matrix
RCO_G10 <- GetAssayData(RCO_10E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(RCO_G10) <- paste("R10", colnames(RCO_G10), sep = "_")


###
# Intersection of highly variable genes among replicates - WT Arabidopsis thaliana
###

# Lets find common variable features for the replicates of A. thaliana
RCOg_hvgs = Reduce(intersect, list(RCO_3E@assays$RNA@var.features, RCO_5E@assays$RNA@var.features, RCO_6E@assays$RNA@var.features, RCO_7E@assays$RNA@var.features, RCO_8E@assays$RNA@var.features, RCO_10E@assays$RNA@var.features))

fileGenerator(RCOg_hvgs, fileName = "Shared_highly_variable_genes_between_replicates_RCOg_AT.txt")


###
# combine the highly variable genes between the genotypes
###

HVGs_combined = union(wt_hvgs, RCOg_hvgs) # Total = 2294

fileGenerator(HVGs_combined, "Seurat_3000_HVG_intersect_reps_union_genotypes_without_mincells.txt")


# Integration of the replicates - find anchors
anchFeatures <- SelectIntegrationFeatures(object.list = list(COL_1E, COL_2E, COL_3E, COL_5E, RCO_3E, RCO_5E, RCO_6E, RCO_7E, RCO_8E, RCO_10E))

fileGenerator(anchFeatures, "Seurat_HVG_standard.txt")


###
# Creating a Liger object, pre-processing, and performing integration
###

# Lets create the liger object
AT_genotypes <- createLiger(list(WC1 = COL_W1, WC2 = COL_W2, WC3 = COL_W3, WC5 = COL_W5, RG3 = RCO_G3, RG5 = RCO_G5, RG6 = RCO_G6, RG7 = RCO_G7, RG8 = RCO_G8, RG10 = RCO_G10), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
AT_genotypes@cell.data$Replicates <- AT_genotypes@cell.data$dataset
AT_genotypes@cell.data$Replicates <- factor(AT_genotypes@cell.data$Replicates, levels = c("WC1", "WC2", "WC3", "WC5", "RG3", "RG5", "RG6", "RG7", "RG8", "RG10"), labels = c("WT-COL-1", "WT-COL-2", "WT-COL-3", "WT-COL-5", "RCOg-3", "RCOg-5", "RCOg-6", "RCOg-7", "RCOg-8", "RCOg-10"))

# Add species information
AT_genotypes@cell.data$Genotype <- str_sub(AT_genotypes@cell.data$dataset, 1, nchar(AT_genotypes@cell.data$dataset) - 1)
AT_genotypes@cell.data$Genotype <- factor(AT_genotypes@cell.data$Genotype, levels = c("WC", "RG"), labels = c("WT", "RCOg"))

# Add tissue information
AT_genotypes@cell.data$Tissue <- "Leaf"
AT_genotypes@cell.data$Tissue <- factor(AT_genotypes@cell.data$Tissue)

AT_genotypes@norm.data <- list(WC1 = COL_1E@assays$RNA@data, 
                               WC2 = COL_2E@assays$RNA@data, 
                               WC3 = COL_3E@assays$RNA@data, 
                               WC5 = COL_5E@assays$RNA@data,
                               RG3 = RCO_3E@assays$RNA@data,
                               RG5 = RCO_5E@assays$RNA@data,
                               RG6 = RCO_6E@assays$RNA@data,
                               RG7 = RCO_7E@assays$RNA@data,
                               RG8 = RCO_8E@assays$RNA@data,
                               RG10 = RCO_10E@assays$RNA@data)

# Selecting a set of highly variable genes
# AT_genotypes <- selectGenes(AT_genotypes, num.genes = 3000, do.plot = FALSE)

AT_genotypes@var.genes <- anchFeatures

# Scale the feature count
AT_genotypes <- scaleNotCenter(AT_genotypes)

# Check which datasets are we integrating
table(AT_genotypes@cell.data$dataset)

# Run liger integration - factorization of the matrices
AT_genotypes <- optimizeALS(AT_genotypes, k = 50, nrep = 10, lambda = 5)

# Quantile normalization of the data - integration in the shared space
AT_genotypes <- quantile_norm(AT_genotypes)

# Run liger implemented UMAP
AT_genotypes <- runUMAP(AT_genotypes)

Liger_object_K_50 <- AT_genotypes

# save the object
save(Liger_object_K_50, file = "integratet_wt_RCO_thaliana_liger.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_wt_RCO_thaliana_liger.txt")