# Load all the packages
library("Seurat")

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

# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_1E)




# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(thaliana_genes, PP_genes)


##### Remove the protoplasting induced genes
COL_data_1E <- COL_data_1E[genes_to_keep, ]
COL_data_2E <- COL_data_2E[genes_to_keep, ]
COL_data_3E <- COL_data_3E[genes_to_keep, ]
COL_data_5E <- COL_data_5E[genes_to_keep, ]

# Create seurat object and perform initial filtering - 
# 1. Remove genes if their expression was not detected in at least one cell out of every 500 ("min.cells"),
# 2. Remove cells if at least 200 genes were not detected to be expressed (min.features = 200),
# 3. Remove cells with a total count of more than 110000 (nCount_RNA > 110000).
# 4. Remove cells if 5% or more of the total count of a cell belongs to the mitochondiral genes.
# 5. Remove cells if 10% or more of the total count of a cell belongs to the chloroplast genes.

###
# COL - 1 E
###

# First replicate - COL 1E - total cells 4850; filter out genes that are not detected in at least 13 cells
COL_1E <- CreateSeuratObject(counts = COL_data_1E, project = "COL_1E", min.cells = 10, min.features = 200)

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

# Find a set of highly avariable genes - 2000 HVGs
COL_1E <- FindVariableFeatures(COL_1E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
COL_W1 <- GetAssayData(COL_1E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W1) <- paste("C1", colnames(COL_W1), sep = "_")

###
# COL - 2 E
###

# First replicate - COL 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
COL_2E <- CreateSeuratObject(counts = COL_data_2E, project = "COL_2E", min.cells = 27, min.features = 200)

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

# Find a set of highly avariable genes - 2000 HVGs
COL_2E <- FindVariableFeatures(COL_2E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
COL_W2 <- GetAssayData(COL_2E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W2) <- paste("C2", colnames(COL_W2), sep = "_")

###
# COL - 3 E
###

# First replicate - COL 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
COL_3E <- CreateSeuratObject(counts = COL_data_3E, project = "COL_3E", min.cells = 6, min.features = 200)

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

# Find a set of highly avariable genes - 2000 HVGs
COL_3E <- FindVariableFeatures(COL_3E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
COL_W3 <- GetAssayData(COL_3E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W3) <- paste("C3", colnames(COL_W3), sep = "_")

###
# COL - 5 E
###

# First replicate - COL 5E - total cells 8420; filter out genes that are not detected in at least 18 cells
COL_5E <- CreateSeuratObject(counts = COL_data_5E, project = "COL_5E", min.cells = 17, min.features = 200)

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

# Find a set of highly avariable genes - 2000 HVGs
COL_5E <- FindVariableFeatures(COL_5E, selection.method = "vst", nfeatures = 2000)

# Extract the count table from the seurat object
COL_W5 <- GetAssayData(COL_5E, assay = "RNA", slot = "counts")

# Prepend a distinctive string to the cell labels 
colnames(COL_W5) <- paste("C5", colnames(COL_W5), sep = "_")



## Select integration features using Seurat's integration feature selection approach

# Integration of the replicates - find features from the datasets to anchor cells from different sources
anchFeatures <- SelectIntegrationFeatures(object.list = list(COL_1E, COL_2E, COL_3E, COL_5E))

fileGenerator(anchFeatures, "Seurat_HVG_standard.txt")


###
# Creating a Liger object, pre-processing, and performing integration
###

# Lets create the liger object
wt_thaliana <- createLiger(list(WC1 = COL_W1, WC2 = COL_W2, WC3 = COL_W3, WC5 = COL_W5), remove.missing = F)

# Add metadata information to the liger object for the replicates in the same way as it was added in the seurat object

# Add replicate information
wt_thaliana@cell.data$Replicates <- wt_thaliana@cell.data$dataset
wt_thaliana@cell.data$Replicates <- factor(wt_thaliana@cell.data$Replicates, levels = c("WC1", "WC2", "WC3", "WC5"), labels = c("WT-COL-1", "WT-COL-2", "WT-COL-3", "WT-COL-5"))

# Add species information
wt_thaliana@cell.data$Genotype <- "WT"
wt_thaliana@cell.data$Genotype <- factor(wt_thaliana@cell.data$Genotype)

# Add tissue information
wt_thaliana@cell.data$Tissue <- "Leaf"
wt_thaliana@cell.data$Tissue <- factor(wt_thaliana@cell.data$Tissue)

wt_thaliana@norm.data <- list(WO1 = OX_1E@assays$RNA@data,
                               WO2 = OX_2E@assays$RNA@data, 
                               WO3 = OX_3E@assays$RNA@data, 
                               WO7 = OX_7E@assays$RNA@data)


wt_thaliana@var.genes <- anchFeatures

# Scale the feature count
wt_thaliana <- scaleNotCenter(wt_thaliana)

# Check which datasets are we integrating
table(wt_thaliana@cell.data$dataset)

# Run liger integration - factorization of the matrices
wt_thaliana <- optimizeALS(wt_thaliana, k = 50, nrep = 10, lambda = 5)

# Quantile normalization of the data - integration in the shared space
wt_thaliana <- quantile_norm(wt_thaliana)

# Run liger implemented UMAP
wt_thaliana <- runUMAP(wt_thaliana)

Liger_object_K_50 <- wt_thaliana

save(Liger_object_K_50, file = "integrated_col_wt_liger.RData")

writeLines(capture.output(sessionInfo()), "Session_info_integrated_ox_wt_liger.txt")