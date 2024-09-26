#!/usr/bin/env Rscript

#- Install Packages
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("DESeq2")
# BiocManager::install("coseq")
# BiocManager::install("apeglm")
# BiocManager::install("HTSFilter")

#- Import Libraris
library(DESeq2)
library(ggplot2)
library(coseq)
library(Cairo)

#- parse arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Es müssen mindestens 2 Argumente übergeben werden: input_file und parameter.")
}

#- Import Data
countData_all <- read.csv(args[1], header = TRUE, sep = "\t")
#head(countData_all)

metaData_all <- read.csv(args[2], header = TRUE, sep = ",")

#head(metaData_all)

metaData_split  <- split(metaData_all, metaData_all$Time_point)

metaData_72h <- metaData_split[["A_72h"]]
metaData_96h <- metaData_split[["B_96h"]]
metaData_120h <- metaData_split[["C_120h"]]

#head(metaData_72h)

ss_72h <- c("X","dmel_72h_A", "dmel_72h_B", "dmel_72h_C", "dmau_72h_A", "dmau_72h_B", "dmau_72h_C")
countData_72h <- countData_all[ ,ss_72h]
ss_96h <- c("X","dmel_96h_A", "dmel_96h_B", "dmel_96h_C", "dmau_96h_A", "dmau_96h_B", "dmau_96h_C")
countData_96h <- countData_all[ ,ss_96h]
ss_120h <- c("X","dmel_120h_A", "dmel_120h_B", "dmel_120h_C", "dmau_120h_A", "dmau_120h_B", "dmau_120h_C")
countData_120h <- countData_all[ ,ss_120h]

#- generate DESeqDataSet from imported data
dds_all <- DESeqDataSetFromMatrix(countData = countData_all,
                              colData = metaData_all,
                              design = ~Organism + Time_point,
                              tidy = TRUE)
dds_72h <- DESeqDataSetFromMatrix(countData = countData_72h,
                              colData = metaData_72h,
                              design = ~Organism,
                              tidy = TRUE)
dds_96h <- DESeqDataSetFromMatrix(countData = countData_96h,
                              colData = metaData_96h,
                              design = ~Organism,
                              tidy = TRUE)
dds_120h <- DESeqDataSetFromMatrix(countData = countData_120h,
                              colData = metaData_120h,
                              design = ~Organism,
                              tidy = TRUE)

dds_all <- DESeq(dds_all)
dds_72h <- DESeq(dds_72h)
dds_96h <- DESeq(dds_96h)
dds_120h <- DESeq(dds_120h)

res_all <- results(dds_all)
res_72h <- results(dds_72h)
res_96h <- results(dds_96h)
res_120h <- results(dds_120h)

#- Table S1
# 72h
# total number of DEGs
upreg_72h_total <- subset(res_72h, padj < 0.05)
write.csv(upreg_72h_total , file = "DEG/upreg_72h_total.csv")
# upregulated in D. mauritiana
upreg_72h_dmau <- subset(res_72h, padj < 0.05 & log2FoldChange > 0)
write.csv(upreg_72h_dmau , file = "DEG/upreg_72h_dmau.csv")


# 96h
# total number of DEGs
upreg_96h_total <- subset(res_96h, padj < 0.05)
write.csv(upreg_96h_total , file = "DEG/upreg_96h_total.csv")
# upregulated in D. mauritiana
upreg_96h_dmau <- subset(res_96h, padj < 0.05 & log2FoldChange > 0)
write.csv(upreg_96h_dmau , file = "DEG/upreg_96h_dmau.csv")

# 120h
# total number of DEGs
upreg_120h_total <- subset(res_120h, padj < 0.05)
write.csv(upreg_120h_total , file = "DEG/upreg_120h_total.csv")

# upregulated in D. mauritiana
upreg_120h_dmau <- subset(res_120h, padj < 0.05 & log2FoldChange > 0)
write.csv(upreg_120h_dmau , file = "DEG/upreg_120h_dmau.csv")



#- filter count Matrix and export
#? not used in R pipeline
library(HTSFilter)

dds_all_filtered_HTS <- HTSFilter(dds_all, plot = TRUE, normalization = "DESeq")
dds_all_filtered <- dds_all_filtered_HTS$filteredData

filtered_countMatrix <- counts(dds_all_filtered)
write.csv(filtered_countMatrix, file = "DEG/filtered_countMatrix.csv")

filtered_countMatrix_collapsed <- counts(collapseReplicates(dds_all_filtered, 
                                                            dds_all_filtered$Time_point, 
                                                            dds_all_filtered$Organism), normalized = TRUE)

#! floor() -> round off collapsed values
write.csv(floor(filtered_countMatrix_collapsed), file = "DEG/filtered_countMatrix_collapsed.csv")



#- generate PCA Plot
rld <- rlog(dds_all, blind = FALSE)
PCA_plot <- plotPCA(rld, intgroup=c("Organism", "Time_point"), ntop = 500)
Cairo::Cairo(file="DEG/PCA.png",
      type="svg", 
      width=15, 
      height=10,
      units="cm",
      dpi=300)
PCA_plot
invisible(dev.off)


#* Check available coefs before continuing to next step
# resultsNames(dds_all)

#? pairwise differential expression analysis between the two species 
#? at each time point using the apeglm option as shrinkage estimator
resLFC_all <- lfcShrink(dds_all, coef = "Organism_dmel_vs_dmau",
                    type = "apeglm")
resLFC_72h <- lfcShrink(dds_72h, coef = "Organism_dmel_vs_dmau",
                    type = "apeglm")
resLFC_96h <- lfcShrink(dds_96h, coef = "Organism_dmel_vs_dmau",
                    type = "apeglm")
resLFC_120h <- lfcShrink(dds_120h, coef = "Organism_dmel_vs_dmau",
                    type = "apeglm")


#- filter for downregulated genes compared to 72h
resLFC_all_96 <- lfcShrink(dds_all, coef = "Time_point_B_96h_vs_A_72h",
                    type = "apeglm")
resLFC_all_120 <- lfcShrink(dds_all, coef = "Time_point_C_120h_vs_A_72h",
                    type = "apeglm")

sigGenes_all_96  <- subset(resLFC_all_96, log2FoldChange < 0 & padj < 0.05)
sigGenes_all_120  <- subset(resLFC_all_120, log2FoldChange < 0 & padj < 0.05)

sigGenes_ids_all_96  <- rownames(sigGenes_all_96)
sigGenes_ids_all_120  <- rownames(sigGenes_all_120)

write.csv(sigGenes_ids_all_96, file = "DEG/sigGenes_ids_all_downreg_96.csv")
write.csv(sigGenes_ids_all_120, file = "DEG/sigGenes_ids_all_downreg_120.csv")

#- combine the read counts that were significantly differentially expressed
#? (log2FC > 0 | log2FC < 0 and padj < 0.05)
#? between the two species in at least one stage
sigGenes_all  <- subset(resLFC_all, abs(log2FoldChange) > 0 & padj < 0.05)
sigGenes_72h  <- subset(resLFC_72h, abs(log2FoldChange) > 0 & padj < 0.05)
sigGenes_96h  <- subset(resLFC_96h, abs(log2FoldChange) > 0 & padj < 0.05)
sigGenes_120h <- subset(resLFC_120h, abs(log2FoldChange) > 0 & padj < 0.05)

sigGenes_ids_all  <- rownames(sigGenes_all)
sigGenes_ids_72h  <- rownames(sigGenes_72h)
sigGenes_ids_96h  <- rownames(sigGenes_96h)
sigGenes_ids_120h <- rownames(sigGenes_120h)

expGenes_all  <- subset(resLFC_all, resLFC_all$log2FoldChange > 0)
expGenes_72h  <- subset(resLFC_72h, abs(log2FoldChange) > 0)
expGenes_96h  <- subset(resLFC_96h, abs(log2FoldChange) > 0)
expGenes_120h <- subset(resLFC_120h, abs(log2FoldChange) > 0)

expGenes_ids_all  <- rownames(expGenes_all)
expGenes_ids_72h  <- rownames(expGenes_72h)
expGenes_ids_96h  <- rownames(expGenes_96h)
expGenes_ids_120h <- rownames(expGenes_120h)

write.csv(sigGenes_ids_all, file = "DEG/sigGenes_ids_all.csv")
write.csv(sigGenes_ids_72h, file = "DEG/sigGenes_ids_72h.csv")
write.csv(sigGenes_ids_96h, file = "DEG/sigGenes_ids_96h.csv")
write.csv(sigGenes_ids_120h, file = "DEG/sigGenes_ids_120h.csv")

#! 'dds_sigGenes' <= dds gets reduced to significantly differentially expressed genes
dds_sigGenes_all  <- dds_all[sigGenes_ids_all, ]
dds_sigGenes_72h  <- dds_all[sigGenes_ids_72h, ]
dds_sigGenes_96h  <- dds_all[sigGenes_ids_96h, ]
dds_sigGenes_120h <- dds_all[sigGenes_ids_120h, ]


#? clustering according to their expression dynamics
#! CAVE - even with seed, clustering results differ !
coseq_sigGenes_all <- coseq::coseq(dds_sigGenes_all, K = 2:25,
                               transformation = "arcsin",
                               norm = "TMM", model = "Normal", parallel = TRUE, seed = 2602112)

sink(file = "DEG/coseq-summary_all.txt")
summary(coseq_sigGenes_all)
sink(file = NULL)

#- Export clusters and profiles
#? clusters ? "FBgnXXXXXXX","#cluster"
clusters_all <- clusters(coseq_sigGenes_all)
write.csv(clusters_all, file = "DEG/coseq-clusters_all.csv")
# profiles ? table with avg. expression profiles
profiles_all <- profiles(coseq_sigGenes_all)
write.csv(profiles_all, file = "DEG/coseq-profiles_all.csv")


#- plotting time <3
#conds_all <- dds_all$Time_point
conds_all <- dds_sigGenes_all$Time_point
print(conds_all)
profiles_plot_all <- plot(coseq_sigGenes_all, graphs = "profiles", conds = conds_all, collapse_reps = "average", order = TRUE)
#profiles_plot_all
Cairo::Cairo(file="DEG/coseq_profile_plot.png",
      type="svg", 
      width=17, 
      height=10,
      units="cm",
      dpi=300)
profiles_plot_all
invisible(dev.off) # to close the file

# #* plot OreR profiles only
# #? zero all values of TAM
# profiles_OreR_all <- profiles_all[1:9]
# profiles_OreR_all[["dmau_120h_A"]] <- 0
# profiles_OreR_all[["dmau_120h_B"]] <- 0
# profiles_OreR_all[["dmau_120h_C"]] <- 0
# profiles_OreR_all[["dmau_72h_A"]] <- 0
# profiles_OreR_all[["dmau_72h_B"]] <- 0
# profiles_OreR_all[["dmau_72h_C"]] <- 0
# profiles_OreR_all[["dmau_96h_A"]] <- 0
# profiles_OreR_all[["dmau_96h_B"]] <- 0
# profiles_OreR_all[["dmau_96h_C"]] <- 0

# profiles_plot_OreR_all <- plot(coseq_sigGenes_all, y_profiles = profiles_OreR_all, graphs = "profiles", conds = conds_all, collapse_reps = "average", order = TRUE)
# profiles_plot_OreR_all

# #* plot TAM profiles only
# #? zero all values of OreR
# profiles_TAM_all <- profiles_all[10:18]
# profiles_TAM_all[["dmel_120h_A"]] <- 0
# profiles_TAM_all[["dmel_120h_B"]] <- 0
# profiles_TAM_all[["dmel_120h_C"]] <- 0
# profiles_TAM_all[["dmel_72h_A"]] <- 0
# profiles_TAM_all[["dmel_72h_B"]] <- 0
# profiles_TAM_all[["dmel_72h_C"]] <- 0
# profiles_TAM_all[["dmel_96h_A"]] <- 0
# profiles_TAM_all[["dmel_96h_B"]] <- 0
# profiles_TAM_all[["dmel_96h_C"]] <- 0

# profiles_plot_TAM_all <- plot(coseq_sigGenes_all, y_profiles = profiles_TAM_all, graphs = "profiles", conds = conds_all, collapse_reps = "average", order = TRUE)
# profiles_plot_TAM_all