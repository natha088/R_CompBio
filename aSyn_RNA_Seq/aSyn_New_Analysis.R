# Load the libraries
library(DESeq2)
library(pathview)
library(org.Hs.eg.db)
library(edgeR)
library(ggplot2)

setwd("/Users/natha088/Documents/Research_Mac/RNA_seq/CHURP_Files/Counts")
# Read the data
counts_table <- read.delim("subread_counts_gene_symbol.txt")
counts_table <- as.data.frame(counts_table)
hmc <- counts_table[,1:10]
sh <- counts_table[,c(1,11:19)]
condition_hmc <- factor(c(rep("H_PBS", 3), rep("H_Mono", 3), rep("H_PFF", 3)))
condition_sh <- factor(c(rep("S_PBS", 3), rep("S_Mono", 3), rep("S_PFF", 3)))
condition <- factor(c(rep("H_PBS", 3), rep("H_Mono", 3), rep("H_PFF", 3), rep("S_PBS", 3), rep("S_Mono", 3), rep("S_PFF", 3)))
dge_hmc <- DGEList(counts = hmc, group = condition_hmc, annotation.columns = 1)
dge_sh <- DGEList(counts = sh, group = condition_sh, annotation.columns = 1)
dge <- DGEList(counts = counts_table, group = condition, annotation.columns = 1)
# Filter lowly expressed genes
keep_hmc <- filterByExpr(dge_hmc,min.count=5)
keep_sh <- filterByExpr(dge_sh,min.count=5)
keep <- filterByExpr(dge,min.count=5)
dge <- dge[keep, ]
dge_hmc <-dge_hmc[keep_hmc, ]
dge_sh <- dge_sh[keep_sh, ]

coldata <- data.frame(row.names = colnames(dge), condition)
coldata_hmc <- data.frame(row.names = colnames(dge_hmc), condition_hmc)
coldata_sh <- data.frame(row.names = colnames(dge_sh), condition_sh)
# Normalize the data (edgeR)
dge <- calcNormFactors(dge)
dge_hmc <- calcNormFactors(dge_hmc)
dge_sh <- calcNormFactors(dge_sh)
# Convert DGEList to DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = dge, colData = coldata, design = ~ condition)
rownames(dds) <- as.matrix(dge$genes)

dds_hmc <- DESeqDataSetFromMatrix(countData = dge_hmc, colData = coldata_hmc, design = ~ condition_hmc)
rownames(dds_hmc) <- as.matrix(dge_hmc$genes)

dds_sh <- DESeqDataSetFromMatrix(countData = dge_sh, colData = coldata_sh, design = ~ condition_sh)
rownames(dds_sh) <- as.matrix(dge_sh$genes)

# Variance stabilizing transformation (VST)
vsd <- vst(dds, blind = FALSE)
vsd_hmc <- vst(dds_hmc, blind = FALSE)
vsd_sh <- vst(dds_sh, blind = FALSE)
# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

pcaData_hmc <- plotPCA(vsd_hmc, intgroup = "condition_hmc", returnData = TRUE)
percentVar_hmc <- round(100 * attr(pcaData_hmc, "percentVar"))

ggplot(pcaData_hmc, aes(PC1, PC2, color = condition_hmc)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_hmc[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_hmc[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

pcaData_sh <- plotPCA(vsd_sh, intgroup = "condition_sh", returnData = TRUE)
percentVar_sh <- round(100 * attr(pcaData_sh, "percentVar"))

ggplot(pcaData_sh, aes(PC1, PC2, color = condition_sh)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar_sh[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar_sh[2], "% variance")) +
  ggtitle("PCA Plot") +
  theme_minimal()

dds <- DESeq(dds)

# Perform contrasts
results_H_Mono_H_PBS <- results(dds, contrast = c("condition", "H_Mono", "H_PBS"))
results_H_PFF_H_PBS <- results(dds, contrast = c("condition", "H_PFF", "H_PBS"))
results_S_Mono_S_PBS <- results(dds, contrast = c("condition", "S_Mono", "S_PBS"))
results_S_PFF_S_PBS <- results(dds, contrast = c("condition", "S_PFF", "S_PBS"))
results_H_PFF_H_Mono <- results(dds, contrast = c("condition", "H_PFF", "H_Mono"))
results_S_PFF_S_Mono <- results(dds, contrast = c("condition", "S_PFF", "S_Mono"))
results_H_Mono_S_Mono <- results(dds, contrast = c("condition", "H_Mono", "S_Mono"))
results_H_PFF_S_PFF <- results(dds, contrast = c("condition", "H_PFF", "S_PFF"))

# MD plot HMC3 
ma_data <- as.data.frame(results_H_PFF_H_PBS)
ma_data$A <- log10(ma_data$baseMean + 1) # Add 1 to avoid log of zero
ma_data$M <- ma_data$log2FoldChange
ma_data$significance <- "Not Significant"
ma_data$significance[ma_data$padj < 0.05 & ma_data$log2FoldChange > 0] <- "Upregulated"
ma_data$significance[ma_data$padj < 0.05 & ma_data$log2FoldChange < 0] <- "Downregulated"
ggplot(ma_data, aes(x = A, y = M)) +
  # First plot non-significant genes in black
  geom_point(data = subset(ma_data, significance == "Not Significant"), 
             color = "black", alpha = 0.6) +
  # Overlay upregulated genes in red
  geom_point(data = subset(ma_data, significance == "Upregulated"), 
             color = "red", alpha = 0.6) +
  # Overlay downregulated genes in blue
  geom_point(data = subset(ma_data, significance == "Downregulated"), 
             color = "blue", alpha = 0.6) +
  # Customize the plot
  theme_minimal() +
  labs(title = "HMC3", 
       x = "Mean Expression (log10)", 
       y = "Log2 Fold Change") +
    ylim(-5, 7.5)
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed")
  
# MD plot SH-SY5Y
ma_data <- as.data.frame(results_S_PFF_S_PBS)
ma_data$A <- log10(ma_data$baseMean + 1) # Add 1 to avoid log of zero
ma_data$M <- ma_data$log2FoldChange
ma_data$significance <- "Not Significant"
ma_data$significance[ma_data$padj < 0.05 & ma_data$log2FoldChange > 0] <- "Upregulated"
ma_data$significance[ma_data$padj < 0.05 & ma_data$log2FoldChange < 0] <- "Downregulated"
ggplot(ma_data, aes(x = A, y = M)) +
  # First plot non-significant genes in black
  geom_point(data = subset(ma_data, significance == "Not Significant"), 
             color = "black", alpha = 0.6) +
  # Overlay upregulated genes in red
  geom_point(data = subset(ma_data, significance == "Upregulated"), 
             color = "red", alpha = 0.6) +
  # Overlay downregulated genes in blue
  geom_point(data = subset(ma_data, significance == "Downregulated"), 
             color = "blue", alpha = 0.6) +
  # Customize the plot
  theme_minimal() +
  labs(title = "SH-SY5Y", 
       x = "Mean Expression (log10)", 
       y = "Log2 Fold Change") +
    ylim(-5, 7.5)
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed")
  
  # Create a volcano plot SH-SY5Y
results_S_PFF_S_PBS$significance <- "Not Significant"
results_S_PFF_S_PBS$significance[results_S_PFF_S_PBS$padj < 0.05 & results_S_PFF_S_PBS$log2FoldChange > 0] <- "Upregulated"
results_S_PFF_S_PBS$significance[results_S_PFF_S_PBS$padj < 0.05 & results_S_PFF_S_PBS$log2FoldChange < 0] <- "Downregulated"
  
ggplot(as.data.frame(results_S_PFF_S_PBS), aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
    theme_minimal() +
    labs(x = "Log2 Fold Change", y = "-log10 Adjusted p-value", title = "SH-SY5Y") +
    theme(legend.title = element_blank())

# Create a volcano plot HMC3
results_H_PFF_H_PBS$significance <- "Not Significant"
results_H_PFF_H_PBS$significance[results_H_PFF_H_PBS$padj < 0.05 & results_H_PFF_H_PBS$log2FoldChange > 0] <- "Upregulated"
results_H_PFF_H_PBS$significance[results_H_PFF_H_PBS$padj < 0.05 & results_H_PFF_H_PBS$log2FoldChange < 0] <- "Downregulated"

ggplot(as.data.frame(results_H_PFF_H_PBS), aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10 Adjusted p-value", title = "HMC3") +
  theme(legend.title = element_blank())

library(ggrepel)

# Ensure gene names match by converting to lowercase
results_H_PFF_H_PBS$gene <- tolower(rownames(results_H_PFF_H_PBS))
genes_to_label <- tolower(c("TNF", "TNFRSF1A", "ADAM17", "RHBDF2", "FRMD8", "BAG4", "TRADD", 
                            "TRAF2", "TRAF5", "RIPK1", "BIRC2", "BIRC3", "RBCK1", "RNF31", 
                            "SHARPIN", "MAP3K5", "MAP3K7", "MAP3K14", "CYLD", "MAP3K8", 
                            "NFKB1", "IKBKB", "IKBKG", "CHUK", "MAP2K4", "MAP2K7", "NFKBIA", 
                            "MAPK8", "MAPK9", "MAPK14", "MAPK3", "MAPK1", "RPS6KA5", "RPS6KA4", 
                            "ITCH", "JUN", "CEBPB", "CREB1", "CFLAR", "RIPK3", "FADD", "CASP8", 
                            "CASP10", "CASP3", "CASP7", "XIAP", "MLKL", "PGAM5", "DNM1L", 
                            "CCL2", "CCL5", "CCL20", "CXCL1", "CXCL2", "CXCL3", "CXCL5", 
                            "CXCL10", "CXCL11", "CSF1", "CSF2", "FAS", "IL18R1", "JAG1", 
                            "IL1B", "IL6", "IL15", "LIF", "BCL3", "SOCS3", "TNFAIP3", "TRAF1", 
                            "FOS", "JUN", "JUNB", "MMP3", "MMP9", "MMP14", "EDN1", "VEGFC", 
                            "NOD2", "ICAM1", "SELE", "VCAM1", "PTGS2", "TNF", "LTA", 
                            "TNFRSF1B", "TRAF1", "TRAF2", "TRAF3", "PIK3CA", "MAP3K14", 
                            "BIRC2", "BIRC3", "WDR1", "AKT1", "RIPK1", "MAPK8", "IKBKB", 
                            "IKBKG", "CHUK", "MAP2K3", "NFKB1", "MAPK14", "JUN", "IRF1", 
                            "IFNB1", "BCL2L1", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", 
                            "FOSL2", "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", "ATF6", "BATF", 
                            "MAF", "MAFA", "MAFB", "MAFG", "MAFF", "MAFK", "NRL"))

# Add a column for gene labels (only label genes in genes_to_label list)
results_H_PFF_H_PBS$label <- ifelse(results_H_PFF_H_PBS$gene %in% genes_to_label, results_H_PFF_H_PBS$gene, NA)

# Volcano plot with labeled genes
ggplot(as.data.frame(results_H_PFF_H_PBS), aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  geom_text_repel(aes(label = label), box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50', max.overlaps = 50) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10 Adjusted p-value", title = "HMC3") +
  ylim(0, 5) + 
  xlim(-1.25, 1.25) +
  theme(legend.title = element_blank())

# Ensure gene names match by converting to lowercase
results_S_PFF_S_PBS$gene <- tolower(rownames(results_S_PFF_S_PBS))
# Add a column for gene labels (only label genes in genes_to_label list)
results_S_PFF_S_PBS$label <- ifelse(results_S_PFF_S_PBS$gene %in% genes_to_label, results_S_PFF_S_PBS$gene, NA)

# Volcano plot with labeled genes
ggplot(as.data.frame(results_S_PFF_S_PBS), aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "black")) +
  geom_text_repel(aes(label = label), box.padding = 0.35, point.padding = 0.3, segment.color = 'grey50', max.overlaps = 50) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-log10 Adjusted p-value", title = "SH-SY5Y") +
  ylim(0, 5) + 
  xlim(-1.25, 1.25) +
  theme(legend.title = element_blank())

# Extract significant genes
sig_genes_H_Mono_H_PBS <- subset(results_H_Mono_H_PBS, padj < 0.05)
sig_genes_H_PFF_H_PBS <- subset(results_H_PFF_H_PBS, padj < 0.05)
sig_genes_S_Mono_S_PBS <- subset(results_S_Mono_S_PBS, padj < 0.05)
sig_genes_S_PFF_S_PBS <- subset(results_S_PFF_S_PBS, padj < 0.05)
sig_genes_H_PFF_H_Mono <- subset(results_H_PFF_H_Mono, padj < 0.05)
sig_genes_S_PFF_S_Mono <- subset(results_S_PFF_S_Mono, padj < 0.05)
sig_genes_H_Mono_S_Mono <- subset(results_H_Mono_S_Mono, padj < 0.05)
sig_genes_H_PFF_S_PFF <- subset(results_H_PFF_S_PFF, padj < 0.05)


# Write to CSV files
write.csv(sig_genes_H_Mono_H_PBS, file = "sig_genes_H_Mono_H_PBS.csv", row.names = TRUE)
write.csv(sig_genes_H_PFF_H_PBS, file = "sig_genes_H_PFF_H_PBS.csv", row.names = TRUE)
write.csv(sig_genes_S_Mono_S_PBS, file = "sig_genes_S_Mono_S_PBS.csv", row.names = TRUE)
write.csv(sig_genes_S_PFF_S_PBS, file = "sig_genes_S_PFF_S_PBS.csv", row.names = TRUE)
write.csv(sig_genes_H_PFF_H_Mono, file = "sig_genes_H_PFF_H_Mono.csv", row.names = TRUE)
write.csv(sig_genes_S_PFF_S_Mono, file = "sig_genes_S_PFF_S_Mono.csv", row.names = TRUE)
write.csv(sig_genes_H_Mono_S_Mono, file = "sig_genes_H_Mono_S_Mono.csv", row.names = TRUE)
write.csv(sig_genes_H_PFF_S_PFF, file = "sig_genes_H_PFF_S_PFF.csv", row.names = TRUE)


# Extract log2 fold changes
log2fc_H_Mono_H_PBS <- results_H_Mono_H_PBS$log2FoldChange
names(log2fc_H_Mono_H_PBS) <- rownames(results_H_Mono_H_PBS)

log2fc_H_PFF_H_PBS <- results_H_PFF_H_PBS$log2FoldChange
names(log2fc_H_PFF_H_PBS) <- rownames(results_H_PFF_H_PBS)

log2fc_S_Mono_S_PBS <- results_S_Mono_S_PBS$log2FoldChange
names(log2fc_S_Mono_S_PBS) <- rownames(results_S_Mono_S_PBS)

log2fc_S_PFF_S_PBS <- results_S_PFF_S_PBS$log2FoldChange
names(log2fc_S_PFF_S_PBS) <- rownames(results_S_PFF_S_PBS)

setwd("/Users/natha088/Documents/Research_Mac/RNA_seq/CHURP_Files/pathview")

# Parkinsons
pathway_id <- "hsa05012"
species <- "human"
# Pathview for each contrast
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species,gene.idtype = "SYMBOL", out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species,gene.idtype = "SYMBOL", out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Neurodegeneration
pathway_id <- "hsa05022"
species <- "human"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Immune signaling
pathway_id <- "hsa04060"
species <- "human"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Dopaminergic synapse
pathway_id <- "hsa04728"
species <- "human"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# JAK/STAT
pathway_id <- "hsa04630"
species <- "human"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Endocytosis
pathway_id <- "hsa04144"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))


# Oxidative phosphorylation
pathway_id<- "hsa00190"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Mitophagy
pathway_id<- "hsa04137"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Autophagy
pathway_id <- "hsa04140"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Fatty acid degradation
pathway_id<- "hsa00071"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# ER Stress
pathway_id <- "hsa04141"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))

# Toll-like
pathway_id <- "hsa04620"
species <- "human"
pathview(gene.data = log2fc_H_Mono_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_Mono_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_H_PFF_H_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "H_PFF_vs_H_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_Mono_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_Mono_vs_S_PBS",limit=c(-2, 2))
pathview(gene.data = log2fc_S_PFF_S_PBS, pathway.id = pathway_id, species = species, gene.idtype = "SYMBOL",out.suffix = "S_PFF_vs_S_PBS",limit=c(-2, 2))


# Extract significant Kegg pathways first (unbiased approach)
convert_to_entrez <- function(res) {
  genes <- rownames(res)
  entrez_ids <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  res$entrez <- entrez_ids
  res <- res[!is.na(res$entrez), ]
  return(res)
}
results_H_Mono_H_PBS_id <- convert_to_entrez(results_H_Mono_H_PBS)
results_H_PFF_H_PBS_id <- convert_to_entrez(results_H_PFF_H_PBS)
results_S_Mono_S_PBS_id <- convert_to_entrez(results_S_Mono_S_PBS)
results_S_PFF_S_PBS_id <- convert_to_entrez(results_S_PFF_S_PBS)
sig_genes_H_Mono_H_PBS_id <- convert_to_entrez(sig_genes_H_Mono_H_PBS)
sig_genes_H_PFF_H_PBS_id <- convert_to_entrez(sig_genes_H_PFF_H_PBS)
sig_genes_S_Mono_S_PBS_id <- convert_to_entrez(sig_genes_S_Mono_S_PBS)
sig_genes_S_PFF_S_PBS_id <- convert_to_entrez(sig_genes_S_PFF_S_PBS)
sig_genes_H_PFF_H_Mono_id <- convert_to_entrez(sig_genes_H_PFF_H_Mono)
sig_genes_S_PFF_S_Mono_id <- convert_to_entrez(sig_genes_S_PFF_S_Mono)

library(GO.db)
library(clusterProfiler)
# Function to perform GO enrichment analysis
perform_go_enrichment <- function(res) {
  gene_list <- res$log2FoldChange
  names(gene_list) <- res$entrez
  gene_list <- sort(gene_list, decreasing = TRUE)
  go_results <- enrichGO(gene = names(gene_list), OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",  pAdjustMethod = "fdr", qvalueCutoff = 0.05)
  return(go_results)
}

go_H_Mono_H_PBS <- perform_go_enrichment(sig_genes_H_Mono_H_PBS_id)
go_H_PFF_H_PBS <- perform_go_enrichment(sig_genes_H_PFF_H_PBS_id)
go_S_Mono_S_PBS <- perform_go_enrichment(sig_genes_S_Mono_S_PBS_id)
go_S_PFF_S_PBS <- perform_go_enrichment(sig_genes_S_PFF_S_PBS_id)
go_H_PFF_H_Mono <- perform_go_enrichment(sig_genes_H_PFF_H_Mono_id)
go_S_PFF_S_Mono <- perform_go_enrichment(sig_genes_S_PFF_S_Mono_id)

go_H_Mono_H_PBS_up <- perform_go_enrichment(sig_genes_H_Mono_H_PBS_id[sig_genes_H_Mono_H_PBS_id$log2FoldChange>0,])
go_H_PFF_H_PBS_up <- perform_go_enrichment(sig_genes_H_PFF_H_PBS_id[sig_genes_H_PFF_H_PBS_id$log2FoldChange>0,])
go_S_Mono_S_PBS_up <- perform_go_enrichment(sig_genes_S_Mono_S_PBS_id[sig_genes_S_Mono_S_PBS_id$log2FoldChange>0,])
go_S_PFF_S_PBS_up <- perform_go_enrichment(sig_genes_S_PFF_S_PBS_id[sig_genes_S_PFF_S_PBS_id$log2FoldChange>0,])
go_H_PFF_H_Mono_up <- perform_go_enrichment(sig_genes_H_PFF_H_Mono_id[sig_genes_H_PFF_H_Mono_id$log2FoldChange>0,])
go_S_PFF_S_Mono_up <- perform_go_enrichment(sig_genes_S_PFF_S_PBS_id[sig_genes_S_PFF_S_Mono_id$log2FoldChange>0,])


go_H_Mono_H_PBS_down <- perform_go_enrichment(sig_genes_H_Mono_H_PBS_id[sig_genes_H_Mono_H_PBS_id$log2FoldChange<0,])
go_H_PFF_H_PBS_down <- perform_go_enrichment(sig_genes_H_PFF_H_PBS_id[sig_genes_H_PFF_H_PBS_id$log2FoldChange<0,])
go_S_Mono_S_PBS_down <- perform_go_enrichment(sig_genes_S_Mono_S_PBS_id[sig_genes_S_Mono_S_PBS_id$log2FoldChange<0,])
go_S_PFF_S_PBS_down <- perform_go_enrichment(sig_genes_S_PFF_S_PBS_id[sig_genes_S_PFF_S_PBS_id$log2FoldChange<0,])
go_H_PFF_H_Mono_down <- perform_go_enrichment(sig_genes_H_PFF_H_Mono_id[sig_genes_H_PFF_H_Mono_id$log2FoldChange<0,])
go_S_PFF_S_Mono_down <- perform_go_enrichment(sig_genes_S_PFF_S_PBS_id[sig_genes_S_PFF_S_Mono_id$log2FoldChange<0,])


dotplot(go_H_Mono_H_PBS, showCategory=10)
dotplot(go_H_PFF_H_PBS, showCategory=10)
dotplot(go_S_Mono_S_PBS, showCategory=10)
dotplot(go_S_PFF_S_PBS, showCategory=10)
dotplot(go_H_PFF_H_Mono, showCategory=20)
dotplot(go_S_PFF_S_Mono, showCategory=10)


dotplot(go_H_Mono_H_PBS_up, showCategory=10)
dotplot(go_H_PFF_H_PBS_up, showCategory=10)
dotplot(go_S_Mono_S_PBS_up, showCategory=10)
dotplot(go_S_PFF_S_PBS_up, showCategory=10)
dotplot(go_H_PFF_H_Mono_up, showCategory=10)
dotplot(go_S_PFF_S_Mono_up, showCategory=10)

dotplot(go_H_Mono_H_PBS_down, showCategory=10)
dotplot(go_H_PFF_H_PBS_down, showCategory=10)
dotplot(go_S_Mono_S_PBS_down, showCategory=10)
dotplot(go_S_PFF_S_PBS_down, showCategory=10)
dotplot(go_H_PFF_H_Mono_down, showCategory=10)
dotplot(go_S_PFF_S_Mono_down, showCategory=10)
