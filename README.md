# R_CompBio
Transcriptomics projects developed in RStudio

# aSyn_RNA_Seq
Generated microglial (HMC3) and differentiated neuron (SH-SY5Y + retinoic acid) samples treated with PBS, aSyn monomer and aSyn pre-formed fibrils in triplicate for RNA-sequencing. The file aSyn_New_Analysis.r contains the RStudio script used to analyze the count matrix table (subread_counts_gene_symbol.txt). Below are some of the plots that were generated for quality control of the RNA-seq data and analysis of the PFF vs. PBS contrast in HMC3's and SH-SY5Y.

1. Principal component analysis (PCA) plot of all samples
![alt text](images/PCA_all.png)

2. PCA plot of HMC3 samples
![alt text](images/PCA_HMC3.png)

3. PCA plot of SH-SY5Y samples
![alt text](images/PCA_SHSY5Y.png)

4. Mean-difference (MD) plot of PFF vs. PBS HMC3 contrast
![alt text](images/MDPlot_HMC3.png)

5. Mean-difference (MD) plot of PFF vs. PBS SH-SY5Y contrast
![alt text](images/MDPlot_SHSY5Y.png)

6. Volcano plot of downregulated and upregulated genes in PFF vs. PBS HMC3 contrast
![alt text](images/VolcanoPlot_HMC3.png)

7. Volcano plot of downregulated and upregulated genes in PFF vs. PBS SH-SY5Y contrast
![alt text](images/VolcanoPlot_SHSY5Y.png)

8. Top 10 Gene Ontology enriched terms in PFF vs. PBS HMC3 contrast
![alt text](images/top10_HMC3.png)

9. Top 10 Gene Ontology enriched terms in PFF vs. PBS SH-SY5Y contrast
![alt text](images/top10_SHSY5Y.png)


