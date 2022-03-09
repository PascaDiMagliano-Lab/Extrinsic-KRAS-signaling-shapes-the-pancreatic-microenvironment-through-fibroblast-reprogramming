#This scipt reproduces the single cell RNAseq analysis of iKras* PanIN lesions as described in Velez-Delgago, et al. 2022.
#Data was processed in line with the Seurat workflow:
#Website: https://satijalab.org/seurat/index.html
#
#Reference: Stuart et al., Comprehensive Integration of Single-Cell Data. 
#Cell, 2019;177(7):1888-1902.e21. PMID: 31178118 PMCID: PMC6687398. 
#
#Figures 3G-H, 4A-E, 5C, and S4A-B are represented in this script. Preprocessing of data for Figures 4A-E is described; final visualizations
#were completed by importing the resulting datatables into Cytoscape V3.7.2 (4A-B) or Circos software V0.69-9 (circos.ca). 
#
#R Version 3.2.2
#Seurat Version 3.2.3

#Load Required Packages:
library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(magrittr)
library(data.table)

#### Load Required Functions--------------------------------------------------------------------------------------------------- #### 
check_genes <- function(genes, database, object_genes) {
  'Check to make sure that wanted genes are in reference you provide and in object
  
   Args:
   genes (chr vector): list of potential wanted genes
   database (data.frame): table with reference ligand/receptor pairs
   object_genes (chr vector): list of genes present in Seurat object
   
   Returns:
   not_found_list (chr list): list of genes not found in database and/or object
  '
  
  database_genes <- as.vector(unlist(database))
  not_found_database <- c()
  not_found_object <- c()
  
  for (i in 1:length(genes)) {
    if (!(genes[i] %in% database_genes)) {
      not_found_database <- c(not_found_database, genes[i])
    }
    
    if (!(genes[i] %in% object_genes)) {
      not_found_object <- c(not_found_object, genes[i])
    }
  }
  
  not_found_list <- list(database=not_found_database, object=not_found_object)
  return(not_found_list)
}

make_LR_pairs <- function(ligands, receptors, database) {
  'Make all LR pairs based on wanted ligands and/or receptors and database provided
  
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   database (data.frame): table with reference ligand/receptor pairs
   
   Returns:
   wanted_LR (data.frame): data.frame with ligand/receptor pairs from wanted ligands and receptors
  '
  
  wanted_LR <- data.frame(Ligand = character(0), Receptor = character(0))
  
  for (ligand in ligands){
    # list of corresponding receptors
    corresponding_receptors <- unique(as.character(database[,2][grep(ligand, database[,1])]))
    
    for (receptor in corresponding_receptors) {
      LR_row <- data.frame(Ligand = ligand, Receptor = receptor)
      wanted_LR <- rbind(wanted_LR, LR_row)
    }
  }
  
  # filter out unwanted receptors
  wanted_LR <- wanted_LR[which(wanted_LR$Receptor %in% wanted_receptors),]
  
  return(wanted_LR)
}

create_LR_table <- function(ligands, receptors, cell_types, LRs, avg0, avg1) {
  'Create table w/ ligand, receptor, source and target cells, average expressions, and IDs
   
   Args:
   ligands (chr vector): list of wanted ligands
   receptors (chr vector): list of wanted receptors
   cell_types (chr vector): list of common cell types between two objects
   LRs (data.frame): table with with wanted ligand/receptor pairs (from make_LR_pairs function)
   avg0 (list of num data.frame): average expression table of object 0
   avg1 (list of num data.frame): average expression table of object 1
   
   Returns:
   LR_table (data.frame): table of potential ligand/receptor pairs
    Contains ligand, receptor, source, target, average expression of ligand and receptors from 
    source and target cells. Also contains IDs (important for Cytoscape). 
  '
  
  LR_table <- data.frame()
  count <- 0
  
  for (ligand in ligands) {
    known_receptors <- LRs$Receptor[which(LRs$Ligand == ligand)]
    
    for (receptor in known_receptors) {
      for (i in c(1:length(cell_types))) {
        for (j in c(1:length(cell_types))) {
          LR_table <- rbind(LR_table,data.frame(ligands = ligand, receptors = receptor, 
                                                source = cell_types[i], target = cell_types[j], 
                                                avg_lig_0 = avg0[ligand, cell_types[i]],
                                                avg_rec_0 = avg0[receptor, cell_types[j]],
                                                avg_lig_1 = avg1[ligand, cell_types[i]],
                                                avg_rec_1 = avg1[receptor, cell_types[j]]))
        }
      }
    }
    
    cat(count, ' ')
    count <- count+1
  }
  
  # Create IDs
  source_ID <- c()
  target_ID <- c()
  
  # Find optimal cell type IDs
  ids_key <- c()
  
  for (i in 1:length(cell_types)) {
    for (j in 1:min(nchar(cell_types))) {
      name_1 <- substring(cell_types[i], 1, j)
      name_list <- c()
      name_list <- c(sapply(cell_types[-i], function(x) name_list <- c(name_list, substring(x, 1, j))))
      
      if(name_1 %in% name_list) {
        next
      }
      else {
        ids_key <- c(ids_key, name_1) 
        break
      }
    }
  }
  
  names(ids_key) <- cell_types
  
  for (i in c(1:length(LR_table[,1]))) {
    letter1 <- as.character(ids_key[LR_table$source[i]])
    letter2 <- as.character(ids_key[LR_table$target[i]])
    n1 <- which(ligands == LR_table$ligands[i])
    n2 <- which(receptors == LR_table$receptors[i])
    
    source_ID <- c(source_ID, paste(letter1, 'L', n1, sep = ''))
    target_ID <- c(target_ID, paste(letter2, 'R', n2, sep = ''))
  }
  
  LR_table <- cbind(LR_table, data.frame(source_ID = source_ID, target_ID = target_ID))
  
  return(LR_table)
}

avg_LR_filt <- function(table, threshold) {
  'Calculates states (ON/OFF --> 1/0) and filter if there is no expression in both groups
  
   Args:
   table (data.frame): table with potential ligand/receptors (from create_LR_table)
   threshold (num): average expression threshold 
   
   Returns:
   table (data.frame): filtered table based on average expression
  '
  
  # Find states of pairs in each group
  LR_states <- data.frame(lig_0 = character(0), rec_0 = character(0),
                          lig_1 = character(0), rec_1 = character(0),
                          is_on = logical(0))
  
  for (i in c(1:length(LR_table$avg_lig_0))) {
    row_states <- data.frame(lig_0 = table$avg_lig_0[i] > threshold, 
                             rec_0 = table$avg_rec_0[i] > threshold,
                             lig_1 = table$avg_lig_1[i] > threshold, 
                             rec_1 = table$avg_rec_1[i] > threshold)
    row_states$is_on <- (row_states$lig_0 & row_states$rec_0) | (row_states$lig_1 & row_states$rec_1)
    
    LR_states <- rbind(LR_states, row_states)
  }
  
  table <- cbind(table, ifelse((LR_states$lig_0 & LR_states$rec_0), 1, 0))
  table <- cbind(table, ifelse((LR_states$lig_1 & LR_states$rec_1), 1, 0))
  
  colnames(table)[11] <- 'is_on_0'
  colnames(table)[12] <- 'is_on_1'
  
  # Filter out pairs if pairs in both group are less than threshold 
  table <- table[which(LR_states$is_on == TRUE),]
  
  return(table)
}

LR_diff <- function(table, data_0, data_1, genes, label, alpha = 0.05) {
  'Calculate Wilcoxon-rank test on ligands and receptors (separately) between groups
    Order of test: data_0 --> data_1
  
    Args:
    table (data.frame): table with potential ligand/receptors (from create_LR_table)
    data_0 (data.frame): table of gene expression from object 0
    data_1 (data.frame): table of gene expression from object 1
    genes (chr vector): list of wanted genes
    label (chr): name of meta.data slot in Seurat object
    alpha (num): alpha level for Wilcoxon-rank test
    
    Returns:
    table (data.frame): filtered table based on Wilcoxon-rank test
  '
  
  table$lig_diff <- rep(0, length(table[,1]))
  table$rec_diff <- rep(0, length(table[,1]))
  
  table$lig_diff_p <- rep(0, length(table[,1]))
  table$rec_diff_p <- rep(0, length(table[,1]))
  
  for (i in 1:length(table[,1])) {
    ligand <- table$ligands[i]
    receptor <- table$receptors[i]
    source <- table$source[i]
    target <- table$target[i]
    
    lig_0_data <- data_0[which(data_0[,label] == source), ligand]
    rec_0_data <- data_0[which(data_0[,label] == target), receptor]
    
    lig_1_data <- data_1[which(data_1[,label] == source), ligand]
    rec_1_data <- data_1[which(data_1[,label] == target), receptor]
    
    lig_wilcox <- wilcox.test(lig_0_data, lig_1_data, exact = F, paired = F)
    rec_wilcox <- wilcox.test(rec_0_data, rec_1_data, exact = F, paired = F)
    
    table$lig_diff_p[i] <- lig_wilcox$p.value
    table$rec_diff_p[i] <- rec_wilcox$p.value
  }
  
  table$lig_diff_p_adj <- p.adjust(table$lig_diff_p, method = 'bonferroni')
  table$rec_diff_p_adj <- p.adjust(table$rec_diff_p, method = 'bonferroni')
  
  # If not significant, then 0
  # If significant, then 1
  for (i in 1:length(table[,1])) {
    table$lig_diff[i] <- ifelse(table$lig_diff_p_adj[i] < alpha, 1, 0)
    table$rec_diff[i] <- ifelse(table$rec_diff_p_adj[i] < alpha, 1, 0)
  }
  
  # If there is difference, then find if increase/decrease for ligands
  for (i in 1:length(table[,1])) {
    if (table$lig_diff[i] == 1) {
      table$lig_diff[i] <- ifelse(table$avg_lig_0[i] > table$avg_lig_1[i], yes = -1, no = 1)
    }
    
    if (table$rec_diff[i] == 1) {
      table$rec_diff[i] <- ifelse(table$avg_rec_0[i] > table$avg_rec_1[i], yes = -1, no = 1)
    }
  }
  
  return(table)
}

generate_supplement <- function(LR_table, cytoscape_nodes) {
  'Generate supplemental information for nodes to input back into Cytoscape
    1) ligand/receptor 2) ligand/receptor name 
  
   Args:
    LR_table (data.frame): table with ligand/receptor pairs
    cytoscape_ids (data.frame): exported cytoscape node table
    
   Returns:
    node_names_ids (data.frame): table with gene names matching IDs and if ligand/receptor
  '
  
  # List of node names
  ids <- c(as.character(LR_table$source_ID), as.character(LR_table$target_ID))
  gene_names <- c(as.character(LR_table$ligands), as.character(LR_table$receptors))
  node_names_ids <- data.frame(ID = ids, names = gene_names, LR = rep(c('L', 'R'),each=length(LR_table$source_ID)),stringsAsFactors = F)
  node_names_ids <- unique(node_names_ids)
  
  # sort supplemental names in order of what is in cytoscape
  cytoscape_ids <- gsub('""', '', cytoscape_nodes$`shared name`)
  order <- match(cytoscape_ids, node_names_ids$ID)
  node_names_ids <- node_names_ids[order,]
  node_names_ids$ID <- paste('"', node_names_ids$ID, '"', sep = '')
  
  return(node_names_ids)
}

#### Preprocessing ------------------------------------------------------------------------------------------------------------ ####

#Load in Datasets:
iKras_on_data <- Read10X("~/raw_feature_bc_matrix/")
iKras_off_data <- Read10X("~/raw_feature_bc_matrix/")
 
#Create Seurat Objects:
iKras_on <- CreateSeuratObject(iKras_on_data, min.cells = 3, min.features = 100)
iKras_off <- CreateSeuratObject(iKras_off_data, min.cells = 3, min.features = 100)

#Add Desired Metadata:
iKras_on[["KRAS"]] <- "ON"
iKras_off[["KRAS"]] <- "OFF"

#Merge Seurat Objects:
iKras <- merge(x = iKras_on, y = iKras_off, add.cell.ids = (c("on", "off")))

#Normalize Data:
iKras <- NormalizeData(object = iKras, normalization.method = "LogNormalize", scale.factor = 10000)

#Apply Unbiased QC Cutoffs:
iKras[["percent.mt"]] <- PercentageFeatureSet(object = iKras, pattern = "^mt-")
iKras <- subset(x = iKras, subset = nCount_RNA > 1000 & nCount_RNA < 60000 & percent.mt < 15)

#Find Variable Genes:
iKras <- FindVariableFeatures(object = iKras, selection.method = "vst", nfeatures = 2000)

#Scale Data:
iKras <- ScaleData(object = iKras, vars.to.regress = "nCount_RNA", features = rownames(iKras)) 

#Run PCA and Determine Dimensions for 90% Variance:
iKras <- RunPCA(iKras, verbose = FALSE)
st_dev <- iKras@reductions$pca@stdev
var <- st_dev^2
sum(var[1:28])/ sum(var) 

#Find Neighbors and Cluster Cells:
iKras <- FindNeighbors(object = iKras, dims = 1:28) 
iKras <- FindClusters(object = iKras, resolution = 1) 

#Run UMAP:  
iKras <- RunUMAP(object = iKras, dims = 1:28) 
DimPlot(object = iKras, reduction = "umap", label = TRUE, pt.size = 0.5)

#Identify Red Blood Cell Clusters:
FeaturePlot(iKras, features = c("Hbb-bt", "Ptprc", "Krt19", "Try4", "Col1a2"), cols = c("gainsboro", "firebrick1"))
DotPlot(iKras, features = rev(c("Hbb-bt", "Ptprc", "Krt19", "Try4", "Col1a2")), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#Label Red Blood Cell Clusters:
iKras <- RenameIdents(iKras,"0" = "Keep", 
                      "1" = "RBC", 
                      "2" = "RBC", 
                      "3" = "Keep", 
                      "4" = "Keep", 
                      "5" = "Keep", 
                      "6" = "Keep", 
                      "7" = "Keep", 
                      "8" = "Keep",
                      "9" = "Keep", 
                      "10" = "Keep", 
                      "11" = "Keep", 
                      "12" = "RBC", 
                      "13" = "Keep", 
                      "14" = "Keep", 
                      "15" = "Keep", 
                      "16" = "Keep", 
                      "17" = "Keep",
                      "18" = "RBC",
                      "19" = "Keep",
                      "20" = "Keep",
                      "21" = "Keep",
                      "22" = "Keep",
                      "23" = "Keep",
                      "24" = "Keep",
                      "25" = "Keep",
                      "26" = "RBC",
                      "27" = "Keep")

#Remove Red Blood Cells:
iKras <- subset(iKras, idents = "Keep")

#### Creating Annotated Global Object ----------------------------------------------------------------------------------------- ####

#Find Variable Genes:
iKras <- FindVariableFeatures(object = iKras, selection.method = "vst", nfeatures = 2000)

#Scale Data:
iKras <- ScaleData(object = iKras, vars.to.regress = "nCount_RNA", features = rownames(iKras)) 

#Run PCA and Determine Dimensions for 90% Variance:
iKras <- RunPCA(iKras, verbose = FALSE)
st_dev <- iKras@reductions$pca@stdev
var <- st_dev^2
sum(var[1:29])/ sum(var) 

#Find Neighbors and Cluster Cells:
iKras <- FindNeighbors(object = iKras, dims = 1:29) 
iKras <- FindClusters(object = iKras, resolution = 3) 

#Run UMAP:
iKras <- RunUMAP(object = iKras, dims = 1:29) 
DimPlot(object = iKras, reduction = "umap", label = TRUE, pt.size = 0.5)

#Identify Clusters Based On Marker Expression:
DotPlot(iKras, features = c("Krt8", "Krt18","Try4", "Col1a2", "Acta2", "Clec3b", "Ptprc", "Cd3e", "Trdc","Cd4", "Foxp3", "Ccr7", "Cd8a", "Nkg7","Il2ra","Il1rl1","Arg1", "Retnla","Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                "Itgae","Clec9a","Batf3", "Ly6c2", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Ccna2", "Pecam1", "Cdh5"), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#Annotate Clusters:
iKras <- RenameIdents(iKras,"0" = "Fibroblasts", 
                      "1" = "Fibroblasts", 
                      "2" = "Fibroblasts", 
                      "3" = "CD8+ T Cells", 
                      "4" = "B Cells", 
                      "5" = "Fibroblasts", 
                      "6" = "Fibroblasts", 
                      "7" = "CD4+ T Cells", 
                      "8" = "Tregs",
                      "9" = "Tregs", 
                      "10" = "Fibroblasts", 
                      "11" = "NK T Cells", 
                      "12" = "Fibroblasts", 
                      "13" = "Plasma Cells", 
                      "14" = "MDSCs", 
                      "15" = "Memory T Cells", 
                      "16" = "CD8+ T Cells", 
                      "17" = "Acinar Cells", 
                      "18" = "Epithelial Cells", 
                      "19" = "Fibroblasts", 
                      "20" = "NK Cells",
                      "21" = "ILC2s", 
                      "22" = "Fibroblasts", 
                      "23" = "Plasma Cells", 
                      "24" = "Macrophages",
                      "25" = "Plasma Cells",
                      "26" = "Macrophages",
                      "27" = "yd T Cells",
                      "28" = "Proliferating Immune Cells",
                      "29" = "yd T Cells",
                      "30" = "Dendritic Cells",
                      "31" = "CD4+ T Cells",
                      "32" = "Monocytes",
                      "33" = "Fibroblasts",
                      "34" = "Mast Cells",
                      "35" = "Endothelial Cells")
iKras[["manual_clusters"]] <- iKras@active.ident
new_order <- c("Epithelial Cells","Acinar Cells", "Endothelial Cells", "Fibroblasts", "Macrophages", "MDSCs", 
               "Monocytes","Dendritic Cells", "CD4+ T Cells","Tregs", "Memory T Cells","CD8+ T Cells",
               "yd T Cells", "NK T Cells", "NK Cells", "ILC2s", "Plasma Cells","B Cells", "Mast Cells", 
               "Proliferating Immune Cells") 
iKras@active.ident <- factor(iKras@active.ident, levels = new_order)
iKras[["manual_clusters"]] <- iKras@active.ident

#Identity Switches (run to switch the identity class of your object): 
Idents(object = iKras) <- "seurat_clusters" 
Idents(object = iKras) <- "KRAS"
Idents(object = iKras) <- "manual_clusters"

# save Seurat object
save(iKras, file = 'iKras.RData')

#### Figures: ------------------------------------------------------------------------------------------------------------------ ####

#### Figure 3G, Global UMAP ####
DimPlot(object = iKras, reduction = "umap", label = F, pt.size = .8, cols = c("Epithelial Cells" = "gold1",
                                                                             "Fibroblasts" = "steelblue1", 
                                                                             "Acinar Cells" = "turquoise4", 
                                                                             "MDSCs" = "mediumblue", 
                                                                             "Macrophages" = "orange1", 
                                                                             "Dendritic Cells" = "gray0", 
                                                                             "Monocytes" = "thistle4", 
                                                                             "CD4+ T Cells" = "plum1",
                                                                             "Tregs" = "maroon4", 
                                                                             "CD8+ T Cells" = "turquoise", 
                                                                             "ILC2s" = "yellow", 
                                                                             "Memory T Cells" = "hotpink",
                                                                             "yd T Cells" = "gold4", 
                                                                             "NK Cells" = "mediumpurple1", 
                                                                             "NK T Cells" = "purple4",
                                                                             "B Cells" = "yellowgreen", 
                                                                             "Plasma Cells" = "wheat1",
                                                                             "Mast Cells" = "firebrick2", 
                                                                             "Endothelial Cells" = "palegreen2", 
                                                                             "Proliferating Immune Cells" = "honeydew3"))

#### Figure 3H, Average Expression Macrophage Heatmap ####
iKras_mac <- subset(iKras, idents = "Macrophages")

mac_genes <- c("Mrc1", "Il1a", "Il1b", "Il6", "Tnf", "Igf1", "Vegfa", "Tgfb1", "Pdcd1lg2", "Arg1", "Il13",
               "Timp1", "Cxcl1", "Pdgfa", "Pdgfb", "Il2", "Ifng", "Irf4", "Irf5", "Mmp2", "Mmp9", "Tgfb2", "Il10",
               "Retnla", "H2-Eb1", "Il4", "Ccr2", "Il18", "Ccl2", "Tnfsf9", "Ccr1", "Apoe", "Trem2", "Chil3", "C1qa",
               "C1qb", "C1qc")

mac_data <- FetchData(iKras_mac, vars = c(mac_genes, "KRAS"))

df_avg <- data.frame()

for (id in levels(factor(mac_data$KRAS))) {
  data_subset <- mac_data %>% filter(KRAS == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-1)], 2, mean)
  df_avg <- rbind(df_avg, data_subset_avg)
}

off_score <- c()
on_score <- c()

for (i in 1:(ncol(df_avg))) {
  if (df_avg[1,i] > df_avg[2,i]){
    ratio_1 <- df_avg[2,i]/df_avg[1,i]
    off_score <- append(off_score, ratio_1)
    on_score <- append(on_score, 1)
  } else {
    ratio_2 <- df_avg[1,i]/df_avg[2,i]
    on_score <- append(on_score, ratio_2)
    off_score <- append(off_score, 1)
  }
  df_avg_ratios <- rbind(df_avg, on_score, off_score)
  df_avg_ratios <- df_avg_ratios[3:4,]
}

colnames(df_avg_ratios) <- colnames(mac_data)[1:(ncol(mac_data)-1)]
rownames(df_avg_ratios) <- levels(factor(mac_data$KRAS))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "YlOrRd"), scale = 'none', cluster_rows = F, cluster_cols = T)

#### Figure 4A/B/C Ligand/Receptor Interactomes ####
#Remove Proliferating Immune Cells:
Idents(object = iKras) <- "manual_clusters"
iKras_interactome <- subset(iKras, idents = c("Epithelial Cells",
                                              "Acinar Cells", 
                                              "Endothelial Cells", 
                                              "Fibroblasts", 
                                              "Macrophages", 
                                              "MDSCs", 
                                              "Monocytes",
                                              "Dendritic Cells", 
                                              "CD4+ T Cells",
                                              "Tregs", 
                                              "Memory T Cells",
                                              "CD8+ T Cells",
                                              "yd T Cells", 
                                              "NK T Cells", 
                                              "NK Cells", 
                                              "ILC2s", 
                                              "Plasma Cells",
                                              "B Cells", 
                                              "Mast Cells"))

#Rename Clusters to two letter abbreviations for Cytoscape compatibility:
iKras_interactome <- RenameIdents(iKras_interactome,
                                  "Epithelial Cells" = "EP",
                                  "Acinar Cells" = "AC", 
                                  "Endothelial Cells" = "EN", 
                                  "Fibroblasts" = "FB", 
                                  "Macrophages" = "MO", 
                                  "MDSCs" = "MD", 
                                  "Monocytes" = "MN",
                                  "Dendritic Cells" = "DC", 
                                  "CD4+ T Cells" = "T4",
                                  "Tregs" = "YQ", 
                                  "Memory T Cells" = "MT",
                                  "CD8+ T Cells" = "T8",
                                  "yd T Cells" = "YD", 
                                  "NK T Cells" = "NT", 
                                  "NK Cells" = "NK", 
                                  "ILC2s" = "IL", 
                                  "Plasma Cells" = "PC",
                                  "B Cells" = "BC", 
                                  "Mast Cells" = "MA")
iKras_interactome[["interactome_names"]] <- iKras_interactome@active.ident
cell_types <- levels(iKras_interactome@active.ident)

#Create Separate Objects for Comparison:
Idents(iKras_interactome) <- "KRAS"

object_0 <- subset(iKras_interactome, idents = "ON")
object_1 <- subset(iKras_interactome, idents = "OFF")

#Import Ligand/Receptor Gene Table and Generate Gene Lists:
LR_Pairs_database <- read.csv("~/Literature Supported Receptor Ligand Pairs Ramilowski with Additions.csv", header = T, stringsAsFactors = F) #Gene list from Ramilowski et al. 2015 PMID: 26198319 plus additions made by our group
LR_Pairs_database <- LR_Pairs_database[,3:4]
ligand_genes <- as.character(LR_Pairs_database$Ligand.ApprovedSymbol.1)
receptor_genes <- as.character(LR_Pairs_database$Receptor.ApprovedSymbol.1)

#Remerge Seurat Objects with Simplified Names:
meta_label <- "interactome_names"
Idents(object_0) <- meta_label
Idents(object_1) <- meta_label
merge <- merge(object_0, object_1)

#Get List of Genes in Merged Object:  
merge_genes <- rownames(merge)

#Create Experimental Group Seurat Objects:
Idents(merge) <- "KRAS"
iKras_ON <- subset(merge, idents = "ON")
iKras_OFF <- subset(merge, idents = "OFF")
Idents(iKras_ON) <- meta_label
Idents(iKras_OFF) <- meta_label

#Remove Ligands and Receptors Genes Not Found in Dataset:
not_found_ligands <- check_genes(ligand_genes, LR_Pairs_database, merge_genes)
not_in_lig <- unique(not_found_ligands$object)

not_found_receptors <- check_genes(receptor_genes, LR_Pairs_database, merge_genes)
not_in_rec <- unique(not_found_receptors$object)

wanted_ligands <- unique(ligand_genes[!ligand_genes %in% not_in_lig])
wanted_receptors <- unique(receptor_genes[!receptor_genes %in% not_in_rec])

genes <- unique(c(wanted_ligands, wanted_receptors))

#Retrieve Data for Genes of Interest:
iKras_ON_data <- FetchData(object = iKras_ON, vars = c(genes, meta_label))
iKras_OFF_data <- FetchData(object = iKras_OFF, vars = c(genes, meta_label))

#Generate Table of All Possible Ligand and Receptor (LR) Pairs:
LRs <- make_LR_pairs(wanted_ligands, wanted_receptors, LR_Pairs_database)

#Calculate Average Expression for All Cell Types in Each Group: 
avg_0 <- AverageExpression(object = iKras_ON, features = genes, verbose = F)
avg_1 <- AverageExpression(object = iKras_OFF, features = genes, verbose = F)

#Create LR Average Expression Table for Each Source/Target Cell Combination:
LR_table <- create_LR_table(wanted_ligands, wanted_receptors, cell_types, LRs, avg_0$RNA, avg_1$RNA)
LR_table[,1:4] <- data.frame(apply(LR_table[,1:4], 2, as.character), stringsAsFactors = FALSE)

#Determine Threshold for Expression (here we use median):
summary(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))
threshold <- median(c(LR_table$avg_lig_0, LR_table$avg_lig_1, LR_table$avg_rec_0, LR_table$avg_rec_1))

#Filter LR pairs based on Average Expresion and Threshold:
LR_table <- avg_LR_filt(LR_table, threshold)

#Calculate p-value between groups:
LR_table <- LR_diff(LR_table, iKras_ON_data, iKras_OFF_data, genes, meta_label)

#Sort LR_table by ligand:
LR_table <- arrange(LR_table, desc(lig_diff != 0 & rec_diff != 0), lig_diff_p_adj)

#Export LR Table for Cytoscape:
write.csv(LR_table, file = "iKras_3wk_on_off_LR_Table_all genes.csv", row.names = F)

#Cytoscape Ver. 3.7.2 Instructions:
#Open Cytoscape
#File -> Import New Network from file -> import LR_table
#souce_ID -> green circle, target_ID -> target symbol
#Export node table to csv via table panel as "all gene node table.csv"
#Import node table into R: 

#Generate Supplemental Data for Cytoscape: 
cytoscape_nodes <- fread("all gene node table.csv")
node_supplement <- generate_supplement(LR_table, cytoscape_nodes)
write.csv(node_supplement, file = "all_gene_supplement.csv", row.names = F)

#In Cytoscape:
#Import  node_supplement via import table from file
#Arrange and color clusters as desired (All ligand receptor pairs here, Figure 4A)
#Filter on Lig_diff_padj for signifigance, p < 0.05, and sort for genes higher in iKras* ON (Figure 4B)

#Filter LR_table to ligand/receptor pairs where ligands originate from fibroblasts and are significantly higher in iKras* ON:
Fibro_LR <- LR_table %>% filter(source == "FB" & lig_diff_p_adj < 0.05 & lig_diff == -1)

#Export Fibro_LR Table for Cicos analysis (Figure 4C):
write.csv(Fibro_LR, file = "iKras_up_in_3wk_on_Fibro_ligands.csv", row.names = F) 

#see additional GitHub upload in 
#https://github.com/PascaDiMagliano-Lab/Extrinsic-KRAS-signaling-shapes-the-pancreatic-microenvironment-through-fibroblast-reprogramming 
#for visualization details

#### Figure 4D, Average Expression Fibroblast Ligands Heatmap ####
Fibro_LR <- LR_table %>% filter(source == "FB") #LR_table generated in interactome analysis in Figure 4A
wanted_genes <- unique(as.character(Fibro_LR$ligands)) 

iKras_fibro <- subset(iKras, idents = "Fibroblasts")
fibro_data <- FetchData(iKras_fibro, vars = c(wanted_genes, "KRAS"))

df_avg <- data.frame()

for (id in levels(factor(fibro_data$KRAS))) {
  data_subset <- fibro_data %>% filter(KRAS == id)
  data_subset_avg <- apply(data_subset[,1:(ncol(data_subset)-1)], 2, mean)
  df_avg <- rbind(df_avg, data_subset_avg)
}

off_score <- c()
on_score <- c()

for (i in 1:(ncol(df_avg))) {
  if (df_avg[1,i] > df_avg[2,i]){
    ratio_1 <- df_avg[2,i]/df_avg[1,i]
    off_score <- append(off_score, ratio_1)
    on_score <- append(on_score, 1)
  } else {
    ratio_2 <- df_avg[1,i]/df_avg[2,i]
    on_score <- append(on_score, ratio_2)
    off_score <- append(off_score, 1)
  }
  df_avg_ratios <- rbind(df_avg, on_score, off_score)
  df_avg_ratios <- df_avg_ratios[3:4,]
}

colnames(df_avg_ratios) <- colnames(fibro_data)[1:(ncol(fibro_data)-1)]
rownames(df_avg_ratios) <- levels(factor(fibro_data$KRAS))
pheatmap(as.matrix(df_avg_ratios), fontsize = 14, color = brewer.pal(9, "YlGnBu"), scale = 'none', cluster_rows = F, cluster_cols = T)


#### Figure 4E, Ligand and Receptor Violin Plots from Global Object ####
plot_colors <- c("Epithelial Cells" = "gold1", "Fibroblasts" = "steelblue1", "Acinar Cells" = "turquoise4", 
                 "MDSCs" = "mediumblue", "Macrophages" = "orange1", "Dendritic Cells" = "gray0", "Monocytes" = "thistle4", 
                 "CD4+ T Cells" = "plum1","Tregs" = "maroon4", "CD8+ T Cells" = "turquoise", "ILC2s" = "yellow", 
                 "Memory T Cells" = "hotpink", "yd T Cells" = "gold4", "NK Cells" = "mediumpurple1", "NK T Cells" = "purple4", 
                 "B Cells" = "yellowgreen", "Plasma Cells" = "wheat1", "Mast Cells" = "firebrick2", 
                 "Endothelial Cells" = "palegreen2", "Proliferating Immune Cells" = "honeydew3")
VlnPlot(iKras, features = "Il33", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Il1rl1", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Cxcl1", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Cxcr2", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Il6", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Il6ra", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Saa3", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "P2rx7", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Saa1", pt.size = 0.1, cols = plot_colors) + NoLegend()
VlnPlot(iKras, features = "Scarb1", pt.size = 0.1, cols = plot_colors) + NoLegend()

#### Figure 5C, Average Expression Epithelial Cells and Fibroblasts ####
iKras_ef <- subset(iKras, idents = c("Epithelial Cells","Fibroblasts"))
Idents(object = iKras_ef) <- "KRAS"
iKras_ef_on <- subset(iKras_ef, idents = "ON")
Idents(object = iKras_ef_on) <- "manual_clusters"
average_iKras_ef <- AverageExpression(object = iKras_ef_on, return.seurat = T)

scale_col <- rev(brewer.pal(11,"RdYlBu"))
DoHeatmap(average_iKras_ef, features =c("Hbegf", "Il6", "Tnf", "Pdgf", "Hgf"), draw.lines = F) + scale_fill_gradientn(colors = scale_col)




#### Figure S4A, Cluster Marker Dotplot ####
DotPlot(iKras, features = rev(c("Krt19", "Try4", "Pecam1", "Cdh5", "Col1a2", "Pdgfrb", "Ptprc", "Cd68","Fcgr3","Adgre1", "Itgam","Cd14","Mrc1", "S100a8", "Cd33","H2-Eb1",
                                "Itgae","Clec9a","Batf3", "Ly6c2", "Ly6g","Cd3e", "Cd4", "Foxp3", "Ccr7", "Cd8a", "Trdc", "Nkg7","Il2ra","Il1rl1", "Cd79a", "Cd19", "Ms4a1","Ccr10","Prdm1","Kit","Mki67")), cols = c("blue", "red"), dot.scale = 8) + RotatedAxis()

#### Figure S4B, Ligand/Receptor Global Dotplot ####
DotPlot(iKras, features = rev(c("P2rx7", "Scarb1","Il1rl1", "Il6ra", "Cxcr2","Saa3","Saa1","Il33", "Il6", "Cxcl1")), cols = c("blue", "red"), dot.scale = 10) + RotatedAxis()
