# set working directory --------------------------------------
# setwd(getwd())
# -----------------------------------------------------------
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("WGCNA")
usePackage("openxlsx")
usePackage("readxl")
usePackage("flashClust")
usePackage("cluster")
usePackage("affy")
usePackage("stats")
usePackage("limma")
# install.packages("BiocInstaller",repos="https://bioconductor.org/packages/3.2/bioc")
# source("http://bioconductor.org/biocLite.R")
# biocLite("multtest")
# biocLite("edgeR")
usePackage("psych")

library(WGCNA)
library(openxlsx)
library(readxl)
library(flashClust)
library(cluster)
library(affy)
library(reshape2)
library(ggplot2)
library(plotly)
library(survival)
library(dplyr)
library(stats)
library(multtest)
library(psych)
library(igraph)
library(edgeR)
library(limma)
library(gridExtra)
library(grid)
library(gtable)


options(stringsAsFactors = FALSE)
# DE and HC ------------------------------------------------
gene.prog <- read.delim(file = "co-expression_2/output/genesprog_HC.txt")
gene.nonprog <- read.delim(file = "co-expression_2/output/genesnonprog_HC.txt")
all.genes <- rbind(gene.prog, gene.nonprog)
annot.clinical.all$bcr.patient.barcode[!which(annot.clinical.all$bcr.patient.barcode %in% rownames(all.genes))]

# make a local copy of original --------------------------------
gnp <- gene.nonprog 
gp <- gene.prog

# clustering ----------------------------------------------------
minModuleSize <- 30
build.clustering <- function(case, softPower, minModuleSize){
  adjacency <- adjacency(case, power = softPower, type = "unsigned") 
 
  # Turn adjacency into topological overlap ----------------------
  TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
  dissTOM <- 1 - TOM

  # Call the hierarchical clustering function ---------------------
  geneTree <- hclust(as.dist(dissTOM), method = "average")
 
  # Module identification using dynamic tree cut ------------------
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = minModuleSize)
  
  dynamicColors <- labels2colors(dynamicMods)
  return(list("geneTree" = geneTree, "dynamicMods" = dynamicMods, "dynamicColors" = dynamicColors, "adj" = adjacency))
}
cluster.prog <- build.clustering(gene.prog, 5, minModuleSize)
cluster.nonprog <- build.clustering(gene.nonprog, 6, minModuleSize)

save(cluster.nonprog, file = "co-expression_3/output/treeGeneNonProg.RData")
save(cluster.prog, file = "co-expression_3/output/treeGeneProg.RData")

# ---------------------------------------------------------------
load(file = "co-expression_3/output/treeGeneNonProg.RData")
load(file = "co-expression_3/output/treeGeneProg.RData")

# Calculate eigengenes ----------------------------------
MEDissThres <- 0.05 # height cut corresponds to 0.95 cor to merge
cnp <- "Non-Progressed"
cp <- "Progressed"

calc.eigen <- function(case, dynamicColors, MEDissThres, c){
  MEList <- moduleEigengenes(case, colors = dynamicColors)
  MEs <- MEList$eigengenes
  # Calculate dissimilarity of module eigengenes ---------
  MEDiss <- 1 - cor(MEs)
  # Cluster module eigengenes ----------------------------
  METree <- hclust(as.dist(MEDiss), method = "average")
  return(list("METree" = METree, "MEList" = MEList, "MEs" = MEs, "MEDiss" = MEDiss))
}
prog.eigen <- calc.eigen(gene.prog, cluster.prog$dynamicColors, MEDissThres, cp)
nonprog.eigen <- calc.eigen(gene.nonprog, cluster.nonprog$dynamicColors, MEDissThres, cnp)

save(prog.eigen, file = "co-expression_3/output/eigenGeneNonProg.RData")
save(nonprog.eigen, file = "co-expression_3/output/eigenGeneProg.RData")

# highly connected and cleaned --------------------------
load(file = "co-expression_3/output/eigenGeneNonProg.RData")
load(file = "co-expression_3/output/eigenGeneProg.RData")

METree.prog.eigen <- prog.eigen$METree
METree.nonprog.eigen <- nonprog.eigen$METree

# Plot the result ----------------------------------------
# high correlations will have low merging height ---------
sizeGrWindow(7, 6)
pdf(file = "co-expression_3/plots/MEtree_coexpr_unsigned.pdf",  wi = 10, he = 8)
plot(METree.prog.eigen, main = paste("Clustering of module eigengenes", "progressed", sep = " "),
     xlab = "", sub = "")
# Plot the cut line into the dendrogram ------------------
abline(h = MEDissThres, col = "red")
plot(METree.nonprog.eigen, main = paste("Clustering of module eigengenes", "Non-progressed", sep = " "),
     xlab = "", sub = "")
# Plot the cut line into the dendrogram ------------------
abline(h = MEDissThres, col = "red")
dev.off()


MEDissThres <- 0.05 # height cut corresponds to 0.95 cor to merge
merging <- function(case, colorDynam){
  # Call an automatic merging function ---------------------
  merge <- mergeCloseModules(case, colorDynam, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors --------------------------------
  mergedColors <- merge$colors 
  # Eigengenes of the new merged modules --------------------
  mergedMEs <- merge$newMEs
  return(list("merge" = merge, "mergedColors" = mergedColors, "mergedMEs" = mergedMEs))
}

# merged color and module eigengenes -------------------------
nonprog.merge <- merging(gene.nonprog, cluster.nonprog$dynamicColors)
prog.merge <- merging(gene.prog, cluster.prog$dynamicColors)


# Plot merged color and module eigengenes -------------------------
sizeGrWindow(12, 9)
pdf(file = "co-expression_3/plots/dendro_unsigned.pdf", wi = 9, he = 6)
plotDendroAndColors(cluster.nonprog$geneTree, cbind(cluster.nonprog$dynamicColors, nonprog.merge$mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    main = "Clustering of NonProgressed",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(cluster.prog$geneTree, cbind(cluster.prog$dynamicColors, prog.merge$mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    main = "Clustering of Progressed",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# module preservation ----------------------------------------------------------
setLabels <- c("Progressed", "Non_Progressed")
multiData <- list(Progressed = list(data = gene.prog), Non.Progressed = list(data = gene.nonprog))
multiColor <- list(Progressed = cluster.prog$dynamicColors, Non.Progressed = cluster.nonprog$dynamicColors)
nSets <- 2

system.time( {
  mp <- modulePreservation(multiData, 
                           multiColor,
                           referenceNetworks = c(1:2),
                           corFnc = "bicor",
                           corOptions = "use = 'p'",
                           networkType = "unsigned",
                           nPermutations = 200,
                           randomSeed = 12345,
                           maxGoldModuleSize = 1000, 
                           maxModuleSize = 1000, 
                           verbose = 5)
} )
system(command = "say job done!")
# user   system  elapsed 
# 6159.608 1079.008 7405.687 

save(mp, file = "co-expression_3/output/coexp_modulePreserv.RData")
# load(file = "output_2/coexp_modulePreserv.RData")

# --------------------------------------------------------------------------
ref <- 1
test <- 2
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality -------------------------------------------
compare.preservation.scores <- cbind(rownames(mp$preservation$observed[[ref]][[test]]), statsObs[, c("medianRank.pres", "medianRank.qual")],
                                     signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) 

openxlsx::write.xlsx(cbind(rownames(mp$preservation$observed[[ref]][[test]]), 
                           mp$preservation$observed[[ref]][[test]]$moduleSize,
                           statsObs[, c("medianRank.pres", "medianRank.qual")],
                           signif(statsZ, 2)), 
                     file = "co-expression_3/preservation_AllZscores_unsigned.xlsx")

preserv.zscores <- cbind(rownames(mp$preservation$observed[[ref]][[test]]), 
                         mp$preservation$observed[[ref]][[test]]$moduleSize,
                         statsObs[, c("medianRank.pres", "medianRank.qual")],
                         signif(statsZ, 2))

non.preserved <- filter(preserv.zscores, 
                        preserv.zscores$Zsummary.pres < 2 & 
                        preserv.zscores$Zdensity.pres < 2 & 
                        preserv.zscores$Zconnectivity.pres < 2 &
                        preserv.zscores$`rownames(mp$preservation$observed[[ref]][[test]])` != "grey")

# Extract Module labels and module sizes  ----------------------------------
modColors <- rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes <- mp$preservation$Z[[ref]][[test]][, 1]
# leave grey and gold modules out ------------------------------------------
plotMods <- !(modColors %in% c("grey", "gold"))
# Text labels for points ---------------------------------------------------
text <- modColors[plotMods]
# Extract each set to plot  ------------------------------------------------
plotData <- cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])

# Main titles for the plot -------------------------------------------------
mains <- c("Preservation Median rank", "Preservation Zsummary")
# Start the plot -----------------------------------------------------------
sizeGrWindow(10, 5) 
pdf(file  = "co-expression_3/plots/modulePreservation_plot.pdf", wi = 10, h = 5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{ 
  min = min(plotData[, p], na.rm = TRUE)
  max = max(plotData[, p], na.rm = TRUE)
  # Adjust ploting ranges ------------------------------------------------
  if (p == 2) {
    if (min > -max/10) {
      min = -max/10
      ylim = c(min - 0.1 * (max - min), max + 0.1 * (max - min)) 
    }
  } else {
    #ylim = c(max + 0.1 * (max - min), min - 0.1 * (max - min)) 
  }
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines between 2 (low) and 10 (high) ----------------------
  if (p == 2) {
    abline(h = 0)
    abline(h = 2, col = "darkblue", lty = 2)
    abline(h = 10, col = "darkgreen", lty = 2)
  }
}
dev.off()

#########################################################################
# consensus Blockwise Auto Generated 
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-auto.pdf
# A 4GB standard desktop or a laptop may handle up to 8000-10000 probes
#########################################################################
setLabels <- c("Progressed", "Non_Progressed")
multiExpr <- list(Progressed = list(data = gene.prog), Non.Progressed = list(data = gene.nonprog))

system.time( bnet <- blockwiseConsensusModules(
  multiExpr, 
  power = 6, 
  corType = "bicor",
  minModuleSize = 30, 
  deepSplit = 2,
  maxBlockSize = 10024, 
  randomSeed = 12345,
  networkType = "unsigned",
  TOMType = "signed",
  TOMDenom = "min",
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.005, #.995 cut for module merging 
  numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, 
  verbose = 5))
system(command = "say job done!")

# user      system  elapsed 
# 308.758  21.787 242.918 

#save(bnet, file = "co-expression_3/consensus/Consensus-NetworkConstruction-auto-entireNetObject.RData")
load("co-expression_3/consensus/Consensus-NetworkConstruction-auto-entireNetObject.RData")
bconsMEs <- bnet$multiMEs
bmoduleLabels <- bnet$colors
bmoduleColors <- labels2colors(bmoduleLabels)
bconsTree <- bnet$dendrograms[[1]]

#bwLabels <- matchLabels(bmoduleColors, bmoduleLabels, pThreshold = 1e-7)
#bwColors <- labels2colors(bwLabels)
table(bmoduleColors)

btb <- table(bmoduleColors)
btb <- as.data.frame(btb)

maxGene <- max(btb$Freq)
minGene <- min(btb$Freq)
counted <- nrow(btb) - 1

head_foot <- function(d, h, f){
  
  title <- textGrob(h, gp=gpar(fontsize=14))
  foot <- textGrob(f, x=0, hjust=0,
                   gp=gpar(fontsize=10, fontface="italic"))
  
  padding <- unit(0.5,"line")
  
  d  <- gtable_add_rows(d , heights = grobHeight(title) + padding,
                        pos = 0)
  d  <- gtable_add_rows(d , heights = grobHeight(foot)+ padding)
  
  d  <- gtable_add_grob(d , list(title,foot),
                        t=c(1, nrow(d)), l=c(1,2), 
                        r=ncol(d))
  return(d)
}

tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
colnames(btb) <- c("moduleColors", "Freq")
bmodule.summary <- tableGrob(t(cbind(btb[1:7,], btb[8:14,], btb[15:21,])), theme = tt)
bmodule.summary <- head_foot(bmodule.summary, "Summary of auto-blockwise consensus module detection\n",
                            paste("Table x: The Dynamic Tree Cut returned", counted,
                                  "proper modules with sizes ranging from\n",  maxGene,
                                  "to", minGene, "genes.The label grey is reserved for genes not assigned to any of \n the modules.", sep = " "))
grid.arrange(bmodule.summary)


# man consensus module scaled heatmap --------------------------------------------------
colnames(bconsMEs[[1]]$data) <- labels2colors(as.numeric(gsub("ME", "", colnames(bconsMEs[[1]]$data))))
colnames(bconsMEs[[2]]$data) <- labels2colors(as.numeric(gsub("ME", "", colnames(bconsMEs[[2]]$data))))

# scaling data to fit elipsoid of 1st PCA eigengene and why
# http://stats.stackexchange.com/questions/110508/questions-on-pca-when-are-pcs-independent-why-is-pca-sensitive-to-scaling-why
# http://stats.stackexchange.com/questions/22329/how-does-centering-the-data-get-rid-of-the-intercept-in-regression-and-pca/22331#22331

# This function produces a red and green color image of a data matrix
# using an RGB color specification. Larger entries are represented with
# reds of increasing intensity, and smaller entries are represented with greens
# of increasing intensity.
# (scale is generic function whose default method centers and/or
# scales the columns of a numeric matrix.)
pdf(file = "co-expression_3/consensus/prog_scaledheatmaps_and_ME.pdf", wi = 8, he = 7)
for(i in 1:length(unique(bmoduleColors))){
  which.module <- unique(bmoduleColors)[[i]]
  datME <- bconsMEs[[1]]$data[which.module]
  colnames(datME) <- paste("ME", which.module, sep = "")
  ME <- datME
  par(mfrow = c(2,1), mar = c(0.3, 5.5, 3, 2))
  plotMat(t(scale(gene.prog[, bmoduleColors == which.module ]) ),
          nrgcols = 30, rlabels = F,rcols = which.module,
          main = which.module, cex.main = 2)
  par(mar = c(5, 4.2, 0, 0.7))
  barplot(t(as.matrix(ME)), col = which.module, main = "", cex.main = 2,
          ylab = "eigengene expression", xlab = "array sample")
}
dev.off()


pdf(file = "co-expression_3/consensus/nonprog_scaledheatmaps_and_ME.pdf", wi = 8, he = 7)
for(i in 1:length(unique(bmoduleColors))){
  which.module <- unique(bmoduleColors)[[i]]
  datME <- bconsMEs[[2]]$data[which.module]
  colnames(datME) <- paste("ME", which.module, sep = "")
  ME <- datME
  par(mfrow = c(2,1), mar = c(0.3, 5.5, 3, 2))
  plotMat(t(scale(gene.nonprog[, bmoduleColors == which.module ]) ),
          nrgcols = 30, rlabels = F,rcols = which.module,
          main = which.module, cex.main = 2)
  par(mar = c(5, 4.2, 0, 0.7))
  barplot(t(as.matrix(ME)), col = which.module, main = "", cex.main = 2,
          ylab = "eigengene expression", xlab = "array sample")
}
dev.off()

# sample.drivers.nonprog <- rownames(gene.nonprog[which( bmoduleColors == "lightyellow")][which(bconsMEs[[2]]$data["lightyellow"] > 0.2),])
# sample.drivers.prog <- rownames(gene.prog[which( bmoduleColors == "lightyellow")][which(bconsMEs[[1]]$data["lightyellow"] > 0.2),])
# "TCGA-BB-4224" "TCGA-BA-5559"
# "TCGA-DQ-7591" "TCGA-BB-4223" "TCGA-DQ-7593"

clinical <- readxl::read_excel(path = "output/clinial_interest_reduced.xlsx", sheet = 1)
# demograph <- clinical[which(clinical$Unique_ID %in% union(sample.drivers.prog, sample.drivers.nonprog)),]

# extract modules gene names --------------------------------------
for(i in 1:length(unique(bmoduleColors))){
  the.module <- unique(bmoduleColors)[i]
  gp <- cluster.prog$adj[which(bmoduleColors %in% the.module) , which(bmoduleColors %in% the.module)]
  write.table(rownames(gp) , file = paste("co-expression_3/consens_module_genes/", the.module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

# Calculate module eigengenes -----------------------------
bunmergedMEs <- multiSetMEs(multiExpr, colors = NULL, universalColors = bmoduleColors)
# Calculate consensus dissimilarity of consensus module eigengenes
bconsMEDiss <- consensusMEDissimilarity(bunmergedMEs)
# Cluster consensus modules --------------------------------
bconsMETree <- hclust(as.dist(bconsMEDiss), method = "average")
# Plot the result ------------------------------------------
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(bconsMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h = 0.25, col = "red")

#  Merging of modules whose expression profiles are very similar ----------------
bmerge <- mergeCloseModules(multiExpr, bmoduleLabels, cutHeight = 0.01, verbose = 3)
# Numeric module labels
MEmoduleLabels <- bmerge$colors
# Convert labels to colors
MEmoduleColors <- labels2colors(MEmoduleLabels)
# Eigengenes of the new merged modules:
consMEs <- bmerge$newMEs

sizeGrWindow(9,6)


pdf(file = "co-expression_3/consensus/ConsensusDendrogram-auto.pdf", wi = 12, he = 9)
plotDendroAndColors(bconsTree, cbind(bmoduleColors, MEmoduleColors),
                    main = "Consensus Clustering of Conditions",
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

grid.arrange(bmodule.summary)

plot(bconsMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h = 0.25, col = "red")
dev.off()


# finding genes kME in modules ---------------------------------------------
bdatME.prog <- bnet$multiMEs[[1]]
bdatKME.prog <- signedKME(gene.prog, bdatME.prog, outputColumnName = "kME.")
kMEcomparisonScatterplot(multiExpr[[1]]$data, multiExpr[[2]]$data, bmoduleColors, nameA = "Progressor", nameB = "Non-pregressor")

bdatME.nonprog <- bnet$multiMEs[[2]]
bdatKME.nonprog <- signedKME(gene.nonprog, bdatME.nonprog, outputColumnName = "kME.")

bm <- melt(bdatKME.prog)
bm$variable <- labels2colors(gsub("kME.ta.", "", bm$variable))
bnm <- melt(bdatKME.nonprog)
bnm$variable <- labels2colors(gsub("kME.ta.", "", bnm$variable))

p1 <- ggplot(data = bm, aes(x = variable, y = value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p2 <- ggplot(data = bnm, aes(x = variable, y = value)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
g1 <- ggplotGrob(p1)
g1 <- gtable_add_cols(g1, unit(0,"mm"))     
g2 <- ggplotGrob(p2)
g2 <- gtable_add_cols(g2, unit(0,"mm"))
g <- rbind(g1, g2, size = "first")                          
g$widths <- unit.pmax(g1$widths, g2$widths) 
g$layout[grepl("guide", g$layout$name), c("t","b")] <- c(1,nrow(g))
grid.newpage()

pdf(file = "co-expression_3/consensus/kME_ranges.pdf", wi = 12, he = 9)
  grid.draw(g)
dev.off()

##############################################################################
# 1. How many genes are DE/DV/DW between progressors and nonprogressors?
# 2. How many DE/DV/DW genes are in progressor modules?
# 3. What is the overlap between each module and the DE/DV/DW genes?
# 4. Is the overlap of each progressor modules genes with DV/DE/DW genes 
#    greater than expected by chance? This is what the Fisher exact test does.
# ----------------------------------------------------------
# DE Which genes are differentially expressed based on intensity (mean)
# optain p-values for mean differences between consditions 
# using ebays library by conducting moderated t-statistics 
# ----------------------------------------------------------
plotmean.list <- lapply(seq_along(unique(bmoduleColors)), function(i){
  the.module <- unique(bmoduleColors)[i]
  module.prog <- gene.prog[, which(bmoduleColors %in% the.module)]
  module.prog <- (2^module.prog) - 1 # unlog and revert shift 
  module.mean.prog <- apply(module.prog, 2, mean)
  
  module.nonprog <- gene.nonprog[, which(bmoduleColors %in% the.module)]
  module.nonprog <- (2^module.nonprog) - 1 # unlog and revert shift 
  module.mean.nonprog <- apply(module.nonprog, 2, mean)
  de.melted <- melt(cbind(module.mean.prog, module.mean.nonprog))
  names(de.melted) <- c("genes", "condition", "mean")
  de.melted <- as.data.frame(de.melted, stringsAsfactors = T)
  de.melted <- de.melted[order(de.melted$mean),]
  
  g <- ggplot(de.melted, aes(x = reorder(genes, +mean), 
                             y = mean, 
                             colour = condition, 
                             group = condition)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = paste("Module", the.module , sep = " "), x = "Gene HGCN", y = "Mean")
  return(g)
})

pdf(file = "co-expression_3/consensus/meanDiff_plots.pdf",  wi = 10, he = 8)
  plotmean.list
dev.off()

plotmean.list <- lapply(seq_along(unique(bmoduleColors)), function(i){
  the.module <- unique(bmoduleColors)[i]
  module.prog <- gene.prog[, which(bmoduleColors %in% the.module)]
  #module.prog <- (2^module.prog) - 1 # unlog and revert shift 
  module.mean.prog <- apply(module.prog, 2, mean)
  
  module.nonprog <- gene.nonprog[, which(bmoduleColors %in% the.module)]
  #module.nonprog <- (2^module.nonprog) - 1 # unlog and revert shift 
  module.mean.nonprog <- apply(module.nonprog, 2, mean)
  de.melted <- melt(cbind(module.mean.prog, module.mean.nonprog))
  names(de.melted) <- c("genes", "condition", "mean")
  de.melted <- as.data.frame(de.melted, stringsAsfactors = T)
  de.melted <- de.melted[order(de.melted$mean),]
  
  g <- ggplot(de.melted, aes(x = reorder(genes, +mean), 
                             y = mean, 
                             colour = condition, 
                             group = condition)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = paste("Module", the.module , sep = " "), x = "Gene HGCN", y = "Mean")
  return(g)
})

pdf(file = "co-expression_3/consensus/meanDifflogged_plots2.pdf",  wi = 10, he = 8)
 plotmean.list
dev.off()

#  Moderated t-statistic --------------------------------------
# http://ugrad.stat.ubc.ca/R/library/limma/html/ebayes.html

prog.unlogged <- (2^gene.prog) - 1 # unlog and revert shift 
prog.means <- apply(prog.unlogged , 2, mean)
nonprog.unlogged <- (2^gene.nonprog) - 1 # unlog and revert shift 
nonprog.means <- apply(nonprog.unlogged , 2, mean)
gene.means <- cbind(prog.means, nonprog.means)
fit <- lmFit(gene.means)
fit <-  eBayes(fit)
# topTable(fit)
#  Ordinary t-statistic
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma

#  Q-Q plots of t statistics
#  Points off the line may be differentially expressed
pdf(file = "co-expression_3/consensus/unloggedmeanQQ_plots_Ebayes.pdf",  wi = 10, he = 8)
  par(mfrow=c(1,2))
  qqt(ordinary.t, df=fit$df.residual, main="Ordinary t")
  abline(0,1)
  qqt(fit$t, df=fit$df.total,main="Moderated t")
  abline(0,1)
  par(mfrow=c(1,1))
dev.off()


m.DE <- lapply(fit$p.value , function(x) x < 0.01)
names(m.DE) <- rownames(fit)
de.genes <- cbind(unlist(lapply(m.DE , function(x) sum(x))), 
                  unlist(lapply(m.DE , function(x) nrow(x))))
de.genes <- names(which(de.genes[,1] == 1))


modules.enrich.DE <- lapply(seq_along(unique(bmoduleColors)), function(i) {
  x <- factor(colnames(all.genes) %in% colnames(gene.prog[, which(bmoduleColors %in% unique(bmoduleColors)[i])]), levels = c("TRUE", "FALSE"))
  y <- factor(colnames(all.genes) %in% de.genes, levels = c("TRUE", "FALSE"))
  p.val <- fisher.test(x, y, alternative = "greater")$p.value
  return(p.val)
})

names(modules.enrich.DE) <- unique(bmoduleColors)
modules.enrich.DE <- unlist(modules.enrich.DE)
DE.modules <- names(modules.enrich.DE)[modules.enrich.DE < (0.05/length(unique(bmoduleColors)))]

# ------------------------------------------------------------
# DV genes
# for each gene in each condition performs an F test to compare the 
# variance of two samples from normal population. 
# library(stats) function var.test()
# ------------------------------------------------------------
dv.list <- lapply(seq_along(unique(bmoduleColors)), function(i){
  the.module <- unique(bmoduleColors)[i]
  module.prog <- gene.prog[, which(bmoduleColors %in% the.module)]
  module.nonprog <- gene.nonprog[, which(bmoduleColors %in% the.module)]
  module.prog <- (2^module.prog) - 1
  module.nonprog <- (2^module.nonprog) - 1
  
  pval.var <- lapply(seq_along(colnames(module.prog)), function(i){
    var.test(x = module.prog[,i], y = module.nonprog[,i])$p.value
  })
  pval.var <- as.data.frame(do.call(rbind, pval.var), stringsAsfactors = F)
  rownames(pval.var) <- colnames(module.prog)
  
  # adjust p-values for simple multiple testing procedure 
  # library(multtest) function mt.rawp2adjp()
  # use benjamini and Hochberg (1995) for strong control of the false 
  # discovery rate (FDR) 
  # FDR of >= .5 might be quite meaningful
  adjp.out <- mt.rawp2adjp(pval.var$V1, proc = "BH")
  fdr.var <- adjp.out$adjp[order(adjp.out$index), 2]
  fdr.var <- as.data.frame(fdr.var, stringsAsfactors = F)
  rownames(fdr.var) <- colnames(module.prog)
  rownames(fdr.var)[fdr.var >= 0.5]
  confidence <- 1 - fdr.var
  
  module.sd.prog <- apply(module.prog, 2, sd)
  module.sd.nonprog <- apply(module.nonprog, 2, sd)
  
  dv.module <- cbind(pval.var, fdr.var, confidence, module.sd.prog, module.sd.nonprog, module.sd.prog - module.sd.nonprog)
  names(dv.module) <- c(paste("pvalVar", the.module, sep = "_"), 
                        paste("FDR", the.module, sep = "_"), 
                        paste("Confidence", the.module, sep = "_"), 
                        paste("progSD", the.module, sep = "_"), 
                        paste("nonprogSD", the.module, sep = "_"),
                        paste("diffSD", the.module, sep = "_"))
  return(dv.module)
})
names(dv.list) <- unique(bmoduleColors)


prog.unlogged <- (2^gene.prog) - 1 # unlog and revert shift 
prog.sd <- apply(prog.unlogged , 2, sd)
nonprog.unlogged <- (2^gene.nonprog) - 1 # unlog and revert shift 
nonprog.sd <- apply(nonprog.unlogged , 2, sd)
gene.sd <- cbind(prog.means, nonprog.means)

pval.var <- lapply(seq_along(colnames(gene.prog)), function(i){
  var.test(x = prog.unlogged[,i], y = nonprog.unlogged[,i])$p.value
})
pval.var <- as.data.frame(do.call(rbind, pval.var), stringsAsfactors = F)
rownames(pval.var) <- colnames(gene.prog)

# adjust p-values for simple multiple testing procedure 
# library(multtest) function mt.rawp2adjp()
# use benjamini and Hochberg (1995) for strong control of the false 
# discovery rate (FDR) 
# FDR of < 0.05 might be quite meaningful
adjp.out <- mt.rawp2adjp(pval.var$V1, proc = "BH")
fdr.var <- adjp.out$adjp[order(adjp.out$index), 2]
fdr.var <- as.data.frame(fdr.var, stringsAsfactors = F)
rownames(fdr.var) <- colnames(gene.prog)
dv.genes <- rownames(fdr.var)[fdr.var < 0.05]
confidence <- 1 - fdr.var
# summary(fdr.var)
# fdr.var        
# Min.   :0.000000  
# 1st Qu.:0.004557  
# Median :0.187677  
# Mean   :0.312369  
# 3rd Qu.:0.585901  
# Max.   :0.999908  

modules.enrich.DV <- lapply(seq_along(unique(bmoduleColors)), function(i) {
  x <- factor(colnames(all.genes) %in% colnames(gene.prog[, which(bmoduleColors %in% unique(bmoduleColors)[i])]), levels = c("TRUE", "FALSE"))
  y <- factor(colnames(all.genes) %in% dv.genes, levels = c("TRUE", "FALSE"))
  p.val <- fisher.test(x, y, alternative = "greater")$p.value
  return(p.val)
})

names(modules.enrich.DV) <- unique(bmoduleColors)
modules.enrich.DV <- unlist(modules.enrich.DV)
DV.modules <- names(modules.enrich.DV)[modules.enrich.DV < (0.05/length(unique(bmoduleColors)))]


plotsd.list <- lapply(seq_along(unique(bmoduleColors)), function(i){
  the.module <- unique(bmoduleColors)[i]
  module.prog <- gene.prog[, which(bmoduleColors %in% the.module)]
  #module.prog <- (2^module.prog) - 1 # unlog and revert shift 
  module.sd.prog <- apply(module.prog, 2, sd)
  
  module.nonprog <- gene.nonprog[, which(bmoduleColors %in% the.module)]
  #module.nonprog <- (2^module.nonprog) - 1 # unlog and revert shift 
  module.sd.nonprog <- apply(module.nonprog, 2, sd)
  dv.melted <- melt(cbind(module.sd.prog, module.sd.nonprog))
  names(dv.melted) <- c("genes", "condition", "sd")
  dv.melted <- as.data.frame(dv.melted, stringsAsfactors = T)
  dv.melted <- dv.melted[order(dv.melted$sd),]
  
  g <- ggplot(dv.melted, aes(x = reorder(genes, +sd), 
                             y = sd, 
                             colour = condition, 
                             group = condition)) + 
    geom_point() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = paste("Module", the.module , sep = " "), x = "Gene HGCN", y = "Standard Deviation")
  return(g)
})

pdf(file = "co-expression_3/consensus/sdDifflogged_plots.pdf",  wi = 10, he = 8)
  plotsd.list
dev.off()

# extract affected modules -----------------------------------
# optimized and re-written code via: https://github.com/iancuo/cosplicingNetworks 
# filter de and dv p_values less than 0.01 and conduct F test 
# modules enriched at p < 0.05, with Bonferroni corrected for multiple module comparison 
# The Fisher test checks whether some modules have more significant genes than expected by chance.
# The expectation is that, if no module is preferentially affected, then the signif genes will be 
# distributed among all modules at random. However, if some modules contain more signif genes than 
# their "fair share", then those modules are flagged as affected.
#-------------------------------------------------------------
modules.enrich.DEDV <- lapply(seq_along(unique(bmoduleColors)), function(i) {
  x <- factor(colnames(all.genes) %in% colnames(gene.prog[, which(bmoduleColors %in% unique(bmoduleColors)[i])]), levels = c("TRUE", "FALSE"))
  y <- factor(colnames(all.genes) %in% union(dv.genes, de.genes), levels = c("TRUE", "FALSE"))
  p.val <- fisher.test(x, y, alternative = "g")$p.value
  return(p.val)
})
names(modules.enrich.DEDV) <- unique(bmoduleColors)
modules.enrich.DEDV <- unlist(modules.enrich.DEDV)
names(modules.enrich.DEDV)[modules.enrich.DEDV < (0.05/length(unique(bmoduleColors)))]

# ------------------------------------------------------------
# DW genes differentially wired
# library(psych) function r.test 
# Tests the significance of a single correlation, the difference 
# between two independent correlations, the difference between two 
# dependent correlations sharing one variable (Williams's Test), 
# or the difference between two dependent correlations with different variables (Steiger Tests).
# ------------------------------------------------------------
pThreshold <- 0.01
adjThreshold <- 0.5
rawAdj1 <- adjacency(gene.prog, power = 5, type = "unsigned") 
rawAdj2 <- adjacency(gene.nonprog, power = 6, type = "unsigned") 
n1 <- dim(gene.prog)[1]
n2 <- dim(gene.nonprog)[1]
geneNames <- rownames(rawAdj1)
nGenes <- length(geneNames)

vectorEdges <- rep(1, nGenes)
names(vectorEdges) <- geneNames

mx <- lapply(geneNames, function(x){
  vector1 <- rawAdj1[x,]
  vector2 <-  rawAdj2[x,]
  idxsVector <- abs(vector1 - vector2) > adjThreshold
  for (gene2 in geneNames[idxsVector == T]){
    vectorEdges[gene2] <- r.test(n = n1, r12 = rawAdj1[x, gene2], r34 = rawAdj2[x, gene2], n2 = n2)$p
  }
  vectorEdges < pThreshold
})

mx <- do.call(rbind, mx)
af.edges <- length(which(mx == T))
colnames(mx)[mx[which(mx == T)]]
edge.change.rate <- af.edges / (length(colnames(gene.prog)))^2


geneChangeEdgeCount <- rowSums(mx)
geneChangeEdgeCount <- as.data.frame(geneChangeEdgeCount, stringsAsfacors = F)
rownames(geneChangeEdgeCount) <- colnames(gene.prog)

pValuesEdgeChange <- rep(1, length(colnames(gene.prog)))
pValuesEdgeChange <- lapply(colnames(gene.prog), function(gene){
  binom.test(x = geneChangeEdgeCount[gene,], n = length(colnames(gene.prog)), 
             p = edge.change.rate, alternative = "g")$p.value
})
names(pValuesEdgeChange) <- colnames(gene.prog)
names(pValuesEdgeChange)[which(unlist(pValuesEdgeChange) != 1)]


adjpOut <- mt.rawp2adjp(unlist(pValuesEdgeChange), proc = "BH")
fdrEdgesChange <- adjpOut$adjp[order(adjpOut$index), 2]

genesChangedEdges <- geneNames[unlist(pValuesEdgeChange) < 0.01]
geneChangeEdgeCount[genesChangedEdges,]
mean(geneChangeEdgeCount$V1)
mean(geneChangeEdgeCount[genesChangedEdges,])

median(geneChangeEdgeCount$geneChangeEdgeCount)
median(geneChangeEdgeCount[genesChangedEdges,])



modules.enrich.DW <- lapply(seq_along(unique(bmoduleColors)), function(i) {
  x <- factor(colnames(all.genes) %in% colnames(gene.prog[, which(bmoduleColors %in% unique(bmoduleColors)[i])]), levels = c("TRUE", "FALSE"))
  y <- factor(colnames(all.genes) %in% genesChangedEdges, levels = c("TRUE", "FALSE"))
  p.val <- fisher.test(x, y, alternative = "g")$p.value
  return(p.val)
})
names(modules.enrich.DW) <- unique(bmoduleColors)
modules.enrich.DW <- unlist(modules.enrich.DW)
DW.modules <- names(modules.enrich.DW)[modules.enrich.DW < (0.05/length(unique(bmoduleColors)))]
DW.modules <- DW.modules[which(!DW.modules %in% "grey")]
DV.modules <- DV.modules[which(!DV.modules %in% "grey")]

dw.genes <- genesChangedEdges
n <- max(length(de.genes), length(dv.genes), length(dw.genes))
length(de.genes) <- n                      
length(dv.genes) <- n
length(dw.genes) <- n
diff.analysis.genes <- cbind(de.genes, dv.genes, dw.genes)
colnames(diff.analysis.genes) <- c("DE_GENES", "DV_GENES", "DW_GENES")
openxlsx::write.xlsx(diff.analysis.genes , file = "co-expression_3/consensus/netDiffAnalysis_genes.xlsx")


n <- max(length(DE.modules), length(DV.modules), length(DW.modules))
length(DE.modules) <- n                      
length(DV.modules) <- n
length(DW.modules) <- n

net.analysis.modules <- rbind(DE.modules, DV.modules, DW.modules)

tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
net.module.summary <- tableGrob(net.analysis.modules, theme = tt)

net.module.summary <- head_foot(net.module.summary, "Summary of Consensus Network analysis\n", 
                                paste("Table  : Modules with DE/DV/DW genes"))

pdf(file = "co-expression_3/consensus/netanalysis_modules.pdf",  wi = 10, he = 8)
 grid.arrange(net.module.summary)
dev.off()


#----------------------------------------------------------------------
# extract modules gene names 
# --------------------------------------------------------------------
for(i in 1:length(unique(bmoduleColors))){
  the.module <- unique(bmoduleColors)[i]
  gp <- cluster.prog$adj[which(bmoduleColors %in% the.module) , which(bmoduleColors %in% the.module)]
  write.table(rownames(gp) , file = paste("co-expression_3/consens_module_genes/", the.module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

#-------------------------------------------------------------------
# extract kME of each module 
# module eigengenes corresponding to the modules returned in colors, in multi-set format.
# -------------------------------------------------------------------
bnet$multiMEs[[1]]$data
pc <- as.numeric(gsub("ME", "", colnames(bnet$multiMEs[[1]]$data)))
colnames(bnet$multiMEs[[1]]$data) <- labels2colors(pc)

bnet$multiMEs[[2]]$data
npc <- as.numeric(gsub("ME", "", colnames(bnet$multiMEs[[2]]$data)))
colnames(bnet$multiMEs[[2]]$data) <- labels2colors(npc)

progkME <- signedKME(gene.prog, bnet$multiMEs[[1]]$data, outputColumnName = "kME.")
colnames(progkME) <- paste("kME.", colnames(bnet$multiMEs[[1]]$data), sep = "") 
nonprogkME <- signedKME(gene.nonprog, bnet$multiMEs[[2]]$data, outputColumnName = "kME.")
colnames(nonprogkME) <- paste("kME.", colnames(bnet$multiMEs[[2]]$data), sep = "") 

prog.topPoskME <- lapply(colnames(progkME), function(x){ 
   orderd.kME <- progkME[order(-progkME[x]),][x]
   return(head(orderd.kME, 20))
})
names(prog.topPoskME) <- colnames(progkME)

prog.topNegkME <- lapply(colnames(progkME), function(x){ 
  orderd.kME <- progkME[order(-progkME[x]),][x]
  return(tail(orderd.kME, 10))
})
names(prog.topNegkME) <- colnames(progkME)

nonprog.topPoskME <- lapply(colnames(nonprogkME), function(x){ 
  orderd.kME <- nonprogkME[order(-nonprogkME[x]),][x]
  return(head(orderd.kME, 10))
})
names(nonprog.topPoskME) <- colnames(nonprogkME)

nonprog.topNegkME <- lapply(colnames(nonprogkME), function(x){ 
  orderd.kME <- nonprogkME[order(-nonprogkME[x]),][x]
  return(tail(orderd.kME, 10))
})
names(nonprog.topNegkME) <- colnames(nonprogkME)

save(prog.topPoskME, prog.topNegkME, nonprog.topPoskME, nonprog.topNegkME, file = "co-expression_3/consensus/consKME.RData")
load("co-expression_3/consensus/consKME.RData")

prog.kME <- lapply(seq_along(prog.topPoskME), function(i){
  g <- as.data.frame(cbind(rownames(prog.topPoskME[[i]]), rownames(prog.topNegkME[[i]])))
  colnames(g) <- c("poskME", "negkME")
  return(g)
})
names(prog.kME) <- names(prog.topPoskME)

nonprog.kME <- lapply(seq_along(nonprog.topPoskME), function(i){
  g <- as.data.frame(cbind(rownames(nonprog.topPoskME[[i]]), rownames(nonprog.topNegkME[[i]])))
  colnames(g) <- c("poskME", "negkME")
  return(g)
})
names(nonprog.kME) <- names(nonprog.topPoskME)

save(prog.kME, nonprog.kME, file = "co-expression_3/consensus/topconsKME.RData")
load("co-expression_3/consensus/topconsKME.RData")

prog.tkME <- lapply(seq_along(prog.topPoskME[!names(prog.topPoskME) %in% "kME.grey"]), function(i){
  g <- as.data.frame(cbind(rownames(prog.topPoskME[[i]]), prog.topPoskME[[i]]))
  names(prog.topPoskME)[[i]]
  colnames(g) <- c(names(prog.topPoskME)[[i]], "value")
  return(g)
})

prog.tkME <- do.call(cbind, prog.tkME)
openxlsx::write.xlsx(prog.tkME, rowNames = F, file = "co-expression_3/consensus/topkME.xlsx")
# -----------------------------------------------------------------------------
# Comparing eigengene networks in progressor and Non-progressor conditions ----
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-EigengeneNetworks.pdf
# -----------------------------------------------------------------------------
nSets <- 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels <- c("Progressor", "Non-Progressor")
shortLabels <- c("Prog", "Non-Prog")

prog.clinical.trait <- as.data.frame(clinical[which(clinical$Unique_ID %in% rownames(gene.prog)), c(5,6)])
rownames(prog.clinical.trait) <- clinical$Unique_ID[which(clinical$Unique_ID %in% rownames(gene.prog))]
nonprog.clinical.trait <- as.data.frame(clinical[which(clinical$Unique_ID %in% rownames(gene.nonprog)), c(5,6)])
rownames(nonprog.clinical.trait) <- clinical$Unique_ID[which(clinical$Unique_ID %in% rownames(gene.nonprog))]


age.prog <- as.data.frame(annot.clinical$birth_days_to[which(annot.clinical$bcr.patient.barcode %in% rownames(gene.prog))])
rownames(age.prog) <- annot.clinical$bcr.patient.barcode[which(annot.clinical$bcr.patient.barcode %in% rownames(gene.prog))]
colnames(age.prog) <- "Age"

age.nonprog <- as.data.frame(annot.clinical$birth_days_to[which(annot.clinical$bcr.patient.barcode %in% rownames(gene.nonprog))])
rownames(age.nonprog) <- annot.clinical$bcr.patient.barcode[which(annot.clinical$bcr.patient.barcode %in% rownames(gene.nonprog))]
colnames(age.nonprog) <- "Age"

prog.clinical.trait <- cbind(prog.clinical.trait, age.prog$Age)
nonprog.clinical.trait <- cbind(nonprog.clinical.trait, age.nonprog$Age)

Trait <- vector(mode = "list", length = nSets)
for (set in 1:1)
{
  Trait[[set]]$data <- prog.clinical.trait
  Trait[[set + 1]]$data <- nonprog.clinical.trait
  names(Trait[[set]]$data) <- c("Smoking", "Alcohol")
  names(Trait[[set + 1]]$data) <- c("Smoking", "Alcohol")
  
}
Trait <- WGCNA::fixDataStructure(Trait)
Trait[[1]]$data <- transform(Trait[[1]]$data, 
          Smoking = as.numeric(Smoking),
          Alcohol = as.numeric(Alcohol))

Trait[[2]]$data <- transform(Trait[[2]]$data,
          Smoking = as.numeric(Smoking),
          Alcohol = as.numeric(Alcohol))

# Recalculate consMEs to give them color names ----------------
consMEsC <- multiSetMEs(multiExpr, universalColors = bmoduleColors)
# We add the trait to the eigengenes and order ---------------
# them by consesus hierarchical clustering 
MET <- consensusOrderMEs(addTraitToMEs(consMEsC, Trait))

sizeGrWindow(8, 10)
pdf(file = "co-expression_3/consensus/EigengeneNetworks2_wage.pdf", width = 8, height = 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()

# textMatrix <- paste(signif(Trait[[1]]$data, 2), "\n(",
#                    signif(Trait[[1]]$data, 1), ")", sep = "")
# dim(textMatrix) <- dim(Trait[[1]]$data)
# par(mar = c(6, 8.8, 3, 2.2))
# labeledHeatmap(Matrix = Trait[[1]]$data ,
#                xLabels = names(Trait[[1]]$data ),
#                yLabels = bmoduleColors,
#                ySymbols = bmoduleColors,
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text = 0.5,
#                zlim = c(-1,1),
#                main = "Module--trait relationships in progressor")

# -----------------------------------------------------------------------
# connectivity and GS 
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-07-Membership.pdf
# -----------------------------------------------------------------------
datME <- bnet$multiMEs[[1]]$data
p.values.smoking <- corPvalueStudent(cor(prog.clinical.trait$tobacco_smoking_pack_years_smoked, datME, use = "p"), nSamples = 68)
GS1.smoking <- as.numeric(cor(prog.clinical.trait$tobacco_smoking_pack_years_smoked, gene.prog, use = "p"))
GeneSignificance.smoking <- abs(GS1.smoking)
# Module significance is defined as average gene significance ------------
ModuleSignificance.smoking <- tapply(GeneSignificance.smoking, bmoduleColors, mean, na.rm = T)

datME <- bnet$multiMEs[[1]]$data
p.values.alcohol <- corPvalueStudent(cor(prog.clinical.trait$alcohol_consumption_per_day, datME, use = "p"), nSamples = 68)
GS1.alcohol <- as.numeric(cor(prog.clinical.trait$alcohol_consumption_per_day, gene.prog, use = "p"))
GeneSignificance.alcohol <- abs(GS1.alcohol)
# Module significance is defined as average gene significance ------------
ModuleSignificance.alcohol <- tapply(GeneSignificance.alcohol, bmoduleColors, mean, na.rm = T)
t1 <- tableGrob(as.data.frame(round(ModuleSignificance.smoking,4)), theme = tt)
grid.arrange(t1)
t1 <- tableGrob(as.data.frame(round(ModuleSignificance.alcohol,4)), theme = tt)
grid.arrange(t1)

round(ModuleSignificance.smoking,4) 
round(ModuleSignificance.alcohol,4) 
sizeGrWindow(8,7)
par(mfrow = c(1,1))
pdf(file = "co-expression_3/consensus/consModuleGS_smokingANDalcohol.pdf", wi = 8, he = 7)
  plotModuleSignificance(GeneSignificance.smoking, bmoduleColors, boxplot = T)
  plotModuleSignificance(GeneSignificance.smoking, bmoduleColors)
  plotModuleSignificance(GeneSignificance.alcohol, bmoduleColors, boxplot = T)
  plotModuleSignificance(GeneSignificance.alcohol, bmoduleColors)
dev.off()

# -------------------------------------------------------------------------------------------
# The function intramodularConnectivity computes the whole network connectivity/degree kTotal,
# the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal
# -------------------------------------------------------------------------------------------
# ADJ1 <- abs(cor(gene.prog, use = "p"))^6
Alldegrees1 <- intramodularConnectivity(cluster.prog$adj, bmoduleColors)
Alldegrees2 <- intramodularConnectivity(cluster.nonprog$adj, bmoduleColors)
openxlsx::write.xlsx(Alldegrees1, rowNames = T, file = "co-expression_3/consensus/prog_networkconnectivity.xlsx")
colorlevels <- unique(bmoduleColors)

# alcohol ------------------------------------------------------------------------------------
pdf(file = "co-expression_3/consensus/consensus_progressorGSalcohol_kwithinConnectivity.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance.alcohol[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kWithin Connectivity", ylab = "Gene Significance drink per day", abline = TRUE)
}

for(i in c(11:length(colorlevels))){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance.alcohol[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kWithin Connectivity", ylab = "Gene Significance drink per day", abline = TRUE)
}
dev.off()

# smoking ------------------------------------------------------------------------------------
pdf(file = "co-expression_3/consensus/consensus_progressorGSsmoking_kwithinConnectivity.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance.smoking[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kWithin Connectivity", ylab = "Gene Significance pack year", abline = TRUE)
}
for(i in c(11:length(colorlevels))){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance.smoking[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kWithin Connectivity", ylab = "Gene Significance pack year", abline = TRUE)
}
dev.off()

# --------------------------------------------------------------------------------------
# Relationship between the module membership measures and intramodular connectivity
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Simulated-07-Membership.pdf
# --------------------------------------------------------------------------------------
names(bdatKME.prog) <- labels2colors(as.numeric(gsub("kME.ta.ME", "", names(bdatKME.prog))))
names(bdatKME.nonprog) <- labels2colors(as.numeric(gsub("kME.ta.ME", "", names(bdatKME.nonprog))))

# progressors -------------------------------------------------------------------------------------
pdf(file = "co-expression_3/consensus/consensus_progKME_over_kwithinConnectivity.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
which.color = colorlevels[[i]]
restrictGenes = bmoduleColors == which.color
verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
                   (bdatKME.prog[restrictGenes, which.color])^6,
                   col = which.color,
                   xlab = "kWithin Intramodular Connec",
                   ylab = "kME (Module Membership)^6")
}

for(i in c(11:length(colorlevels))){
  which.color = colorlevels[[i]]
  restrictGenes = bmoduleColors == which.color
  verboseScatterplot(Alldegrees1$kWithin[restrictGenes],
                     (bdatKME.prog[restrictGenes, which.color])^6,
                     col = which.color,
                     xlab = "kWithin Intramodular Connec",
                     ylab = "kME (Module Membership)^6")
}
dev.off()

# non-progressors -----------------------------------------------------------------------------------
pdf(file = "co-expression_3/consensus/consensus_nonprogKME_over_kwithinConnectivity.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
  which.color = colorlevels[[i]]
  restrictGenes = bmoduleColors == which.color
  verboseScatterplot(Alldegrees2$kWithin[restrictGenes],
                     (bdatKME.nonprog[restrictGenes, which.color])^6,
                     col = which.color,
                     xlab = "kWithin Intramodular Connec",
                     ylab = "kME (Module Membership)^6")
}

for(i in c(11:length(colorlevels))){
  which.color = colorlevels[[i]]
  restrictGenes = bmoduleColors == which.color
  verboseScatterplot(Alldegrees2$kWithin[restrictGenes],
                     (bdatKME.nonprog[restrictGenes, which.color])^6,
                     col = which.color,
                     xlab = "kWithin Intramodular Connec",
                     ylab = "kME (Module Membership)^6")
}
dev.off()

#------kme and GS 
pdf(file = "co-expression_3/consensus/consensus_GSalcohol_KME.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(abs(bdatKME.prog[restrict1, whichmodule]),
                     GeneSignificance.alcohol[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kME Connectivity", ylab = "Gene Significance alcohol per day", abline = TRUE)
}

for(i in c(11:length(colorlevels))){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(abs(bdatKME.prog[restrict1, whichmodule]),
                     GeneSignificance.alcohol[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kME Connectivity", ylab = "Gene Significance alcohol per day", abline = TRUE)
}
dev.off()

pdf(file = "co-expression_3/consensus/consensus_GSsmoking_KME.pdf", wi = 13, he = 9)
par(mfrow = c(2, as.integer(0.5 + 10 / 2)))
par(mar = c(4,5,3,1))
for(i in c(1:10)){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(abs(bdatKME.prog[restrict1, whichmodule]),
                     GeneSignificance.smoking[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kME Connectivity", ylab = "Gene Significance pack years", abline = TRUE)
}

for(i in c(11:length(colorlevels))){
  whichmodule = colorlevels[[i]]
  restrict1 = (bmoduleColors == whichmodule)
  verboseScatterplot(abs(bdatKME.prog[restrict1, whichmodule]),
                     GeneSignificance.smoking[restrict1], col = bmoduleColors[restrict1],
                     main = whichmodule,
                     xlab = "kME Connectivity", ylab = "Gene Significance pack years", abline = TRUE)
}
dev.off()
# -------------------------------------------------------------------------------------------
# Finding genes with high gene significance and high intramodular connectivity 
# -------------------------------------------------------------------------------------------
proggenes.highGS.highkME <- lapply(names(bdatKME.prog), function(x){
  the.module <- x
  filter <- abs(GeneSignificance.alcohol) > 0.2 & abs(bdatKME.prog[the.module]) > 0.8
  table(filter)
  filter <- as.data.frame(dimnames(multiExpr[[1]]$data)[[2]][filter])
  names(filter) <- the.module
  return(filter)
})
names(proggenes.highGS.highkME) <- names(bdatKME.prog)

max_length <- max(sapply(proggenes.highGS.highkME,nrow))
proggenes.highGS.highkME <- sapply(seq_along(proggenes.highGS.highkME), function(i){
  c(proggenes.highGS.highkME[[i]][[1]], rep(NA, max_length - length(proggenes.highGS.highkME[[i]][[1]])))
})
proggenes.highGS.highkME <- as.data.frame(proggenes.highGS.highkME)
names(proggenes.highGS.highkME) <- names(bdatKME.prog)

openxlsx::write.xlsx(proggenes.highGS.highkME, rowNames = F, file = "co-expression_3/consensus/progGenes_highGCandkME_alcohol.xlsx")

for(i in 1:length(names(proggenes.highGS.highkwithin))){
  the.module <- names(proggenes.highGS.highkwithin)[i]
  the.names <- as.data.frame(proggenes.highGS.highkME[,i][!is.na(proggenes.highGS.highkwithin[,i])])
  write.table(the.names, file = paste("co-expression_3/proggenes_highGSandkME_alcohol/", the.module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

#---------------------------------------------------------------------------------
# high KME and smoking 
#---------------------------------------------------------------------------------
proggenes.highGS.highkwithin <- lapply(names(bdatKME.prog), function(x){
  the.module <- x
  filter <- abs(GeneSignificance.smoking) > 0.2 & abs(bdatKME.prog[the.module]) > 0.8
  table(filter)
  filter <- as.data.frame(dimnames(multiExpr[[1]]$data)[[2]][filter])
  names(filter) <- the.module
  return(filter)
})
names(proggenes.highGS.highkwithin) <- names(bdatKME.prog)

max_length <- max(sapply(proggenes.highGS.highkwithin,nrow))
proggenes.highGS.highkwithin <- sapply(seq_along(proggenes.highGS.highkwithin), function(i){
  c(proggenes.highGS.highkwithin[[i]][[1]], rep(NA, max_length - length(proggenes.highGS.highkwithin[[i]][[1]])))
})
proggenes.highGS.highkwithin <- as.data.frame(proggenes.highGS.highkwithin)
names(proggenes.highGS.highkwithin) <- names(bdatKME.prog)

openxlsx::write.xlsx(proggenes.highGS.highkwithin, rowNames = F, file = "co-expression_3/consensus/progGenes_highGCandkME_smoking.xlsx")

for(i in 1:length(names(proggenes.highGS.highkwithin))){
  the.module <- names(proggenes.highGS.highkwithin)[i]
  the.names <- as.data.frame(proggenes.highGS.highkwithin[,i][!is.na(proggenes.highGS.highkwithin[,i])])
  write.table(the.names, file = paste("co-expression_3/proggenes_highGSandkME_smoking/", the.module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

#---------------------------------------------------------------------------------
# alcohol and smoking 
#---------------------------------------------------------------------------------
proggenes.highGS.highkwithin <- lapply(names(bdatKME.prog), function(x){
  the.module <- x
  filter <- abs(GeneSignificance.alcohol) > 0.2 & abs(GeneSignificance.smoking) > 0.2 & abs(bdatKME.prog[the.module]) > 0.8
  table(filter)
  filter <- as.data.frame(dimnames(datExpr[[1]]$data)[[2]][filter])
  names(filter) <- the.module
  return(filter)
})
names(proggenes.highGS.highkwithin) <- names(bdatKME.prog)

max_length <- max(sapply(proggenes.highGS.highkwithin,nrow))
proggenes.highGS.highkwithin <- sapply(seq_along(proggenes.highGS.highkwithin), function(i){
  c(proggenes.highGS.highkwithin[[i]][[1]], rep(NA, max_length - length(proggenes.highGS.highkwithin[[i]][[1]])))
})
proggenes.highGS.highkwithin <- as.data.frame(proggenes.highGS.highkwithin)
names(proggenes.highGS.highkwithin) <- names(bdatKME.prog)

openxlsx::write.xlsx(proggenes.highGS.highkwithin, rowNames = T, file = "co-expression_3/consensus/progGenes_highGCandKwithin_both.xlsx")

for(i in 1:length(names(proggenes.highGS.highkwithin))){
  the.module <- names(proggenes.highGS.highkwithin)[i]
  the.names <- as.data.frame(proggenes.highGS.highkwithin[,i][!is.na(proggenes.highGS.highkwithin[,i])])
  write.table(the.names, file = paste("co-expression_3/proggenes_highGSandkME_both/", the.module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

# -------------------------------------------------------------------------------------------
# Overlapping samples with anatomical region, HPV, plates and center 
# -------------------------------------------------------------------------------------------
annot.clinical <- read_excel("Data/HNSCC_Table1_Data.xlsx", col_names = TRUE, col_types = NULL, na = "", skip = 0)
annot.clinical[[1]] <- gsub(".", '-', annot.clinical[[1]], fixed = T)
clinical.trait <- annot.clinical[, c(1, 7, 8)]
raw.clinial.dat <- read.delim(file = "Data/nationwidechildrens.org_clinical_patient_hnsc.txt", header = T, sep = "\t")
raw.clinial.dat <- raw.clinial.dat[-c(1:2), ]
raw.clinial.nonprog <- raw.clinial.dat[which(raw.clinial.dat$bcr_patient_barcode %in% rownames(gene.nonprog)),]
raw.clinial.prog <- raw.clinial.dat[which(raw.clinial.dat$bcr_patient_barcode %in% rownames(gene.prog)),]

maps <- read.delim(file = "Data/mappingcombo.txt")
maps <- cbind(maps, substr(maps$barcode.s., 22, 28), substr(maps$barcode.s., 6, 7)) # extract plate numbers 
names(maps)[5:6] <- c("plate-center", "TSS")

submaps <- maps[, c("simple_barcode", "plate-center", "TSS")]
names(clinical.trait)[1] <- "simple_barcode"
submaps <- unique(submaps)
clinical.trait <- merge(clinical.trait, submaps, by = "simple_barcode")

prog.clinical.overlap <- clinical.trait[which(clinical.trait$simple_barcode %in% rownames(gene.prog)), ]
nonprog.clinical.overlap <- clinical.trait[which(clinical.trait$simple_barcode %in% rownames(gene.nonprog)), ]

rownames(prog.clinical.overlap) <- prog.clinical.overlap$simple_barcode
prog.clinical.overlap <- as.data.frame(prog.clinical.overlap)

rownames(nonprog.clinical.overlap) <- nonprog.clinical.overlap$simple_barcode
nonprog.clinical.overlap <- as.data.frame(nonprog.clinical.overlap)
dim(prog.clinical.overlap)
for (i in 2:5) {
  prog.clinical.overlap[[i]] <- as.numeric(factor(prog.clinical.overlap[[i]]))
  nonprog.clinical.overlap[[i]] <- as.numeric(factor(nonprog.clinical.overlap[[i]]))
}
# -------------------------------------------------------------------------------------------
# https://wikis.utexas.edu/display/bioiteam/Clustering+using+WGCNA
# -------------------------------------------------------------------------------------------
sampleTree.prog <- flashClust(dist(gene.prog), method = "average") 
traitColors.prog <- numbers2colors(prog.clinical.overlap[,2:5], signed = FALSE)
sampleTree.nonprog <- flashClust(dist(gene.nonprog), method = "average") 
traitColors.nonprog <- numbers2colors(nonprog.clinical.overlap[,2:5], signed = FALSE)
names(prog.clinical.overlap)[2:5] <- c("hpv_p16_ish", "anatomy_region", "plates-center", "TSS")
names(nonprog.clinical.overlap)[2:5] <- c("hpv_p16_ish", "anatomy_region", "plates-center", "TSS")

pdf(file = "co-expression_3/consensus/sample_pheno_overlap.pdf", wi = 13, he = 9)
  plotDendroAndColors(sampleTree.prog, traitColors.prog, groupLabels = names(prog.clinical.overlap)[2:5], main = "Progressor Sample Average Dissimilarity Dendrogram and phynotype category overlap heatmap")
  plotDendroAndColors(sampleTree.nonprog, traitColors.nonprog, groupLabels = names(nonprog.clinical.overlap)[2:5], main = "None-Progressor Sample Average Dissimilarity Dendrogram and phynotype category overlap heatmap")
dev.off()

# -------------------------------------------------------------------------------------------
#  extract samples driving ME's in modules  
# -------------------------------------------------------------------------------------------
datME <- bnet$multiMEs[[1]]$data
rownames(datME) <- rownames(gene.prog)
progsample.highkME <- lapply(seq_along(names(bdatKME.prog)), function(i){
  the.module <- names(bdatKME.prog)[[i]]
  filter <- abs(bnet$multiMEs[[1]]$data[the.module]) > 0.2
  filter <- as.data.frame(rownames(datME)[filter])
  names(filter) <- the.module
  return(filter)
})
names(progsample.highkME) <- names(bdatKME.prog)

datME <- bnet$multiMEs[[2]]$data
rownames(datME) <- rownames(gene.nonprog)
nonprogsample.highkME <- lapply(seq_along(names(bdatKME.nonprog)), function(i){
  the.module <- names(bdatKME.nonprog)[[i]]
  filter <- abs(bnet$multiMEs[[2]]$data[the.module]) > 0.2
  filter <- as.data.frame(rownames(datME)[filter])
  names(filter) <- the.module
  return(filter)
})
names(nonprogsample.highkME) <- names(bdatKME.nonprog)

pdf(file = "co-expression_3/consensus/highMEprogSample_pheno_overlap.pdf", wi = 13, he = 9)
for(i in 1:length(names(progsample.highkME))){
  if(length(progsample.highkME[[i]][[1]]) >= 2){
  subset.dat <- gene.prog[which(rownames(gene.prog) %in% progsample.highkME[[i]][[1]]),]
  sampleTree.subset.dat <- flashClust(dist(subset.dat), method = "average") 
  traitColors.subset.dat <- numbers2colors(prog.clinical.overlap[rownames(subset.dat),2:5], signed = FALSE)
  plotDendroAndColors(sampleTree.subset.dat, traitColors.subset.dat, 
                      groupLabels = names(prog.clinical.overlap)[2:5], 
                      main = paste("Progressor ", names(progsample.highkME)[[i]], 
                                   " Consensus Module Average Dissimilarity Dendrogram and phynotype category overlap heatmap", 
                                   sep = ""))
  }
}
dev.off()

# -----------------------------------------------------------------------
# Code to export modules to cytoscape 
# Our analysis data is too large for visualization purposes
# -----------------------------------------------------------------------
setLabels <- c("Progressed", "Non_Progressed")
datExpr <- list(Progressed = list(data = gene.prog), Non.Progressed = list(data = gene.nonprog))

TOM <- TOMsimilarityFromExpr(gene.prog, power = 5)
TOM.nonprog <- TOMsimilarityFromExpr(gene.nonprog, power = 6)

# Co-expression DE/DV/DW modules from net-analysis 
#modules <- c("brown","turquoise","pink","purple", "black", "salmon","grey60", "blue", "lightyellow")
modules <- c("turquoise")
probes <- names(gene.prog)
inModule <- is.finite(match(bmoduleColors, modules))
modProbes <- gene.prog[inModule]
modGenes <- colnames(modProbes)
# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule]

dimnames(modTOM) <- list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("co-expression_3/cyto_export/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("co-expression_3/cyto_export/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,
                               nodeNames = modGenes,
                               altNodeNames = modGenes,
                               nodeAttr = bmoduleColors[inModule])

# ----------------------------------------------------------------
# simple heatmap of traits 
# ----------------------------------------------------------------
sizeGrWindow(10,6)
# Will display correlations and their p-values
dat.trait <- as.data.frame(cbind(prog.clinical.trait$tobacco_smoking_pack_years_smoked,
      prog.clinical.trait$alcohol_consumption_per_day))
names(dat.trait) <- c("pack years", "drink per day")
moduleTraitCor <- cor(dat.trait, bconsMEs[[1]]$data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = 68)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               yLabels = names(dat.trait),
               xLabels = names(bconsMEs[[1]]$data),
               xSymbols = names(bconsMEs[[1]]$data),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste(""))

# -----------------------------------------------------------------
# KM curve of prog vs unprog 
library(survival)
# ----------------------------------------------------------------
annot.clinical.all <- read_excel("Data/HNSCC_Table1_Data.xlsx", col_names = TRUE, col_types = NULL, na = "", skip = 0)
annot.clinical.all[[1]] <- gsub(".", '-', annot.clinical.all[[1]], fixed = T)
annot.clinical.all <- annot.clinical.all[which(annot.clinical.all$bcr.patient.barcode %in% rownames(all.genes)),]

annot.clinical.all$last_contact_days_min[annot.clinical.all$last_contact_days_min == "NA"] <- NA
annot.clinical.all$death_days_final[annot.clinical.all$death_days_final == "NA"] <- NA
deaths.days <- as.numeric(do.call(rbind, lapply(seq_along(annot.clinical.all$vital_status_final), function(i) ifelse(!is.na(as.numeric(annot.clinical.all$death_days_final[[i]])), annot.clinical.all$death_days_final[[i]], annot.clinical.all$last_contact_days_min[[i]]))))

annot.clinical.all$vital_status_final[which(annot.clinical.all$vital_status_final == "Alive")] <- 0
annot.clinical.all$vital_status_final[which(annot.clinical.all$vital_status_final == "Dead")] <- 1

formula <- Surv(deaths.days, as.numeric(annot.clinical.all$vital_status_final)) ~ annot.clinical.all$Progression_FINAL
fit <- survfit(formula,  type = "kaplan-meier", conf.type = "log")
p <- survdiff(formula)
p.val <- 1 - pchisq(p$chisq, length(p$n) - 1)

quantile(fit, probs = c(0.25, 0.5, 0.75), conf.int = FALSE)
summary(fit, times = c(50, 100))
plot(fit, col = c(1:4), xlab = "time (days)", ylab = "survival probability", main = paste("Kaplan-Meier-estimate of tumor progression; p-value:", round(p.val, 4)))
d <- summary(fit, times = c(50, 100))

tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
KM.quantile <- tableGrob(quantile(fit, probs = c(0.25, 0.5, 0.75), conf.int = FALSE), theme = tt)
KM.summary <- tableGrob(d$table, theme = tt)
grid.arrange(KM.summary)

length(annot.clinical.all$anatomic_organ_subdivision[which(annot.clinical.all$anatomic_organ_subdivision == "Larynx")])
length(annot.clinical.all$anatomic_organ_subdivision[which(annot.clinical.all$anatomic_organ_subdivision == "Oral Cavity")])
length(annot.clinical.all$anatomic_organ_subdivision[which(annot.clinical.all$anatomic_organ_subdivision == "Oropharynx")])

length(annot.clinical.all$hpv_status_p16_or_ish[which(annot.clinical.all$hpv_status_p16_or_ish == "NA")])
length(annot.clinical.all$hpv_status_p16_or_ish[which(annot.clinical.all$hpv_status_p16_or_ish == "Negative")])
length(annot.clinical.all$hpv_status_p16_or_ish[which(annot.clinical.all$hpv_status_p16_or_ish == "Positive")])

# -----------------------------------------------------------------
# smoking and alcohol use barplot 
# ----------------------------------------------------------------
raw.clinial.dat <- read.delim(file = "Data/nationwidechildrens.org_clinical_patient_hnsc.txt", header = T, sep = "\t")
raw.clinial.dat <- raw.clinial.dat[which(raw.clinial.dat$bcr_patient_barcode %in% rownames(all.genes)),]
length(raw.clinial.dat$tobacco_smoking_history_indicator[which(raw.clinial.dat$tobacco_smoking_history_indicator == "Current reformed smoker for < or = 15 years")])
length(raw.clinial.dat$tobacco_smoking_history_indicator[which(raw.clinial.dat$tobacco_smoking_history_indicator == "Current smoker")])
length(raw.clinial.dat$tobacco_smoking_history_indicator[which(raw.clinial.dat$tobacco_smoking_history_indicator == "Lifelong Non-smoker")])
length(raw.clinial.dat$tobacco_smoking_history_indicator[which(raw.clinial.dat$tobacco_smoking_history_indicator == "Current reformed smoker for > 15 years")])

raw.clinial.dat[raw.clinial.dat == "[Not Available]"] <- NA
length(which(!is.na(raw.clinial.dat$tobacco_smoking_pack_years_smoked)))
length(which(!is.na(raw.clinial.dat$alcohol_consumption_per_day)))

write.table(names(all.genes)[which(names(all.genes) %in% map.to.gene$hgnc)] , file = "co-expression_3/cosplice_coexpr_commongenes.txt", quote = F, sep = "\n", append = F, row.names = F, col.names = F)

p <- clinical$hpv_status_p16_or_ish[which(clinical$Unique_ID %in% rownames(gene.nonprog))] %in% "Positive"
length(which(p == T))

dat <- as.data.frame(cbind(raw.clinial.dat$tobacco_smoking_pack_years_smoked, 
                           raw.clinial.dat$alcohol_consumption_per_day))
names(dat) <- c("pack_years", "drink_per_day")
dat <- as.data.frame(dat, stringsAsfactors = T)
dat[!is.na(dat)] <- "complete"
dat[is.na(dat)] <- "NA"
dat <- tidyr::gather(dat, "pack_years", "drink_per_day", 1:2)
names(dat) <- c("variable", "documentation")
ggplot(dat, aes(documentation, fill=documentation)) + geom_bar() + facet_grid(. ~ variable) + xlab("") + ylab("Frequency")

# ---------------------------------------------------------
# DE , DV, DW heatmaps 
# ---------------------------------------------------------
unique(raw.clinial.dat$anatomic_organ_subdivision)
unique(annot.clinical$anatomic_organ_subdivision)

p.o <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Oral Cavity" & annot.clinical$bcr.patient.barcode %in% rownames(gene.prog))]
p.l <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Larynx" & annot.clinical$bcr.patient.barcode %in% rownames(gene.prog))]
p.op <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Oropharynx" & annot.clinical$bcr.patient.barcode %in% rownames(gene.prog))]


np.o <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Oral Cavity" & annot.clinical$bcr.patient.barcode %in% rownames(gene.nonprog))]
np.l <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Larynx" & annot.clinical$bcr.patient.barcode %in% rownames(gene.nonprog))]
np.op <- annot.clinical$bcr.patient.barcode[which(annot.clinical$anatomic_organ_subdivision == "Oropharynx" & annot.clinical$bcr.patient.barcode %in% rownames(gene.nonprog))]



prog.unlogged <- (2^gene.prog) - 1 # unlog and revert shift 
prog.means <- apply(prog.unlogged , 2, mean)
nonprog.unlogged <- (2^gene.nonprog) - 1 # unlog and revert shift 
nonprog.means <- apply(nonprog.unlogged , 2, mean)
gene.means <- cbind(prog.means, nonprog.means)



de.dv.dw <- read_excel("co-expression_3/consensus/netDiffAnalysis_genes.xlsx", col_names = TRUE, col_types = NULL, na = "", skip = 0)
dw <- de.dv.dw$DW_GENES[!is.na(de.dv.dw$DW_GENES)]
de <- de.dv.dw$DE_GENES[!is.na(de.dv.dw$DE_GENES)]
dv <- de.dv.dw$DV_GENES[!is.na(de.dv.dw$DV_GENES)]


dw.nonprog <- as.matrix(gene.nonprog[rownames(gene.nonprog) %in% np.o, colnames(gene.nonprog) %in% dw])
dw.prog <- as.matrix(gene.prog[rownames(gene.prog) %in% p.o, colnames(gene.prog) %in% dw])
turq <- prog.unlogged[rownames(prog.unlogged) %in% p.o, rownames(cluster.prog$adj)[unique(bmoduleColors) == "turquoise"]]
turq <- melt(turq)
dw.prog <- melt(dw.prog)
dw.nonprog <- melt(dw.nonprog)

ggplot(turq, aes(Var2, Var1)) + geom_tile(aes(fill = value)) + xlab("") + ylab("")
ggplot(turq, aes(Var2, Var1)) + geom_tile(aes(fill = value)) + xlab("") + ylab("") 
ggplot(dw.prog, aes(Var2, Var1)) + geom_tile(aes(fill = value)) + xlab("") + ylab("")             
ggplot(dw.nonprog, aes(Var2, Var1)) + geom_tile(aes(fill = value)) + xlab("") + ylab("")

heatmap.2(dw.prog, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))

#write.table(as.data.frame(dw) , file = "co-expression_3/netanalysis/dw.txt", quote = F, sep = "\n", append = F, row.names = F, col.names = F)
#write.table(as.data.frame(de) , file = "co-expression_3/netanalysis/de.txt", quote = F, sep = "\n", append = F, row.names = F, col.names = F)
#write.table(as.data.frame(dv) , file = "co-expression_3/netanalysis/dv.txt", quote = F, sep = "\n", append = F, row.names = F, col.names = F)

library("pheatmap")
library("RColorBrewer")
all.dw <- cbind(t(dw.prog), t(dw.nonprog))

pheatmap(all.dw, cluster_cols = F, cluster_rows = F)
pheatmap(dw.prog, cluster_cols = F, cluster_rows = F)
pheatmap(dw.nonprog, cluster_cols = F, cluster_rows = F)

annot.clinical <- dplyr::filter(annot.clinical.all, annot.clinical.all$bcr.patient.barcode %in% rownames(all.genes))
annot <- as.data.frame(cbind(annot.clinical$bcr.patient.barcode, annot.clinical$anatomic_organ_subdivision))
names(annot) <- c("id", "anatomy")

the.test <- all.genes[match(annot.clinical$bcr.patient.barcode, rownames(all.genes)),]
test <- merge(annot, the.test, by.x = "id", by.y = 0)
rownames(test) <- rownames(the.test)
conditions <- c(rownames(gene.prog) ,rownames(gene.nonprog))
the.test <- test[match(conditions, rownames(test)), ]
the.test <- cbind(c(rep("progressor", 68) , rep("nonprogressor", 161)), the.test)
names(the.test)[1] <- "condition"
the.test[,2] <- NULL

library(reshape)
testing <- reshape::melt(the.test[,c("condition",dw)], id.var = "condition")
rm(testing)

testing <- the.test[,dw]
annotation <- data.frame(testing = factor(0:68 == 1, labels = c("progressor", "nonprogressor")))
pheatmap(t(the.test[which(the.test$condition %in%"progressor" & rownames(the.test) %in% p.o),dw]), cluster_rows = T, cluster_cols = F)
pheatmap(t(the.test[which(the.test$condition %in%"progressor" & rownames(the.test) %in% p.l),dw]), cluster_rows = T, cluster_cols = F)
pheatmap(t(the.test[which(the.test$condition %in%"progressor" & rownames(the.test) %in% p.op),dw]), cluster_rows = T, cluster_cols = F)



pheatmap(t(the.test[which(the.test$condition %in%"progressor"), dw]), cluster_rows = T, cluster_cols = T)
pheatmap(t(the.test[which(the.test$condition %in%"nonprogressor"),dw]), cluster_rows = T, cluster_cols = T)

pheatmap(t(the.test[which(the.test$condition %in% "progressor"), de[1000:1050]]), cluster_rows = T, cluster_cols = F)
pheatmap(t(the.test[which(the.test$condition %in% "nonprogressor"), de[200:250]]), cluster_rows = T, cluster_cols = F)

pheatmap(t(the.test[which(the.test$condition %in% "progressor"), dv[1000:1050]]), cluster_rows = T, cluster_cols = F)
pheatmap(t(the.test[which(the.test$condition %in% "nonprogressor"), dv[1000:1050]]), cluster_rows = T, cluster_cols = F)


annotation <- as.data.frame(the.test$condition, stringsAsFactors = T)
pheatmap(the.test[,dw], cluster_rows = F, cluster_cols = F, annotation = annotation, gp = gpar(fill = "green"))

tttest <- melt(the.test)
DEgenes.paper <- as.data.frame(DEgenes.paper)
# write.table(DEgenes.paper , file = "co-expression_3/netanalysis/paper_degenes.txt", quote = F, sep = "\n", append = F, row.names = F, col.names = F)
DEgenes.paper[grep("IT" , DEgenes.paper),]
#the.test[rownames(the.test) %in% c("TCGA-CR-5248", "TCGA-CN-5358", "TCGA-BA-5559", "TCGA-BA-4074") ,1:2]
summary(the.test[which(the.test$condition %in% "progressor"), "TNFRSF1B"])

low.median <- c("TCGA-BA-7269", "TCGA-CV-A45Q", "TCGA-D6-A6EM", "TCGA-CV-5973", "TCGA-CV-6436", "TCGA-D6-A6ES", "TCGA-D6-A6EK")
submaps[which(submaps$simple_barcode %in% low.median),]

batch.effect <- tableGrob(submaps[which(submaps$simple_barcode %in% low.median),], theme = tt)
grid.arrange(batch.effect)
