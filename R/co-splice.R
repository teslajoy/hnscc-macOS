usePackage <- function(p)
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}
usePackage("WGCNA")
usePackage("multtest")
usePackage("openxlsx")
usePackage("readxl")
usePackage("limma")
usePackage("plotly")
usePackage("survival")
usePackage("dplyr")
usePackage("stats")
usePackage("igraph")
usePackage("edgeR")
usePackage("gridExtra")
usePackage("grid")
usePackage("gtable")


library(WGCNA)
library(multtest)
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
library(vegan)

# ----------------------------------------------------------------------
# load function definition and similarity matrices of All, prog, non-prog 
# ----------------------------------------------------------------------
options(stringsAsFactors = FALSE)
# ----------------------------------------------------------------------
load("co-spliceEx/cormantel-matrices/corrMantelNonProg_WOlfactory_cleaned.RData")
load("co-spliceEx/cormantel-matrices/corrMantelProg_WOlfactory_cleaned.RData")
load("final_cospliceEx/corrMantelAll_WOlfactory_cleaned.RData")


# ----------------------------------------------------------------------
# calculate sft threshhold 
# ----------------------------------------------------------------------
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
diag(mantelCorMatrix) <- 1
dim(mantelCorMatrix) # 10614 10614

sft.all <- pickSoftThreshold.fromSimilarity(similarity = mantelCorMatrix, 
                                            RsquaredCut = 0.85, 
                                            powerVector = powers, 
                                            verbose = 5, 
                                            moreNetworkConcepts = T, 
                                            blockSize = 10614)

diag(mantelCorMatrixProg) <- 1
sft.prog <- pickSoftThreshold.fromSimilarity(similarity = mantelCorMatrixProg, 
                                             RsquaredCut = 0.85, 
                                             powerVector = powers, 
                                             verbose = 5, 
                                             moreNetworkConcepts = T, 
                                             blockSize = 10614)

diag(mantelCorMatrixNonProg) <- 1
sft.nonprog <- pickSoftThreshold.fromSimilarity(similarity = mantelCorMatrixNonProg, 
                                                RsquaredCut = 0.85, 
                                                powerVector = powers, 
                                                verbose = 5, 
                                                moreNetworkConcepts = T, 
                                                blockSize = 10614)


# ----------------------------------------------------------------
# plot sftThresh calculations for all, progressor, and non progressor
# ----------------------------------------------------------------
sizeGrWindow(9, 5)
plot.soft.thr <- function(sft, case){
  par(mfrow = c(1,2))
  cex1 = 0.9
  
  # Scale-free topology fit index as a function of the soft-thresholding power -----------
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n",
       main = paste("Scale independence", case, sep = " "))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red")
  
  # red line corresponds to using an R^2 cut-off of h -----------------------------------
  abline(h = 0.90, col = "red")
  
  # Mean connectivity as a function of the soft-thresholding power ----------------------
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity", 
       type = "n",
       main = paste("Mean connectivity", case, sep = " "))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  par(mfrow = c(1,1))
}

pdf(file = "co-spliceEx/plots/sftthresh.pdf", wi = 12, he = 9)
  plot.soft.thr(sft.all, "co-spliceEx All Data")      #5
  plot.soft.thr(sft.prog, "co-spliceEx Progressor")   #7
  plot.soft.thr(sft.nonprog, "co-spliceEx NonProgressor") #7
dev.off()

#----------------------------------------------------------------
# clustering
#----------------------------------------------------------------
minModuleSize <- 30
build.clustering <- function(case, softPower, minModuleSize){
  adjacency <- case^softPower
  adjacency[is.na(adjacency)] = 0
  adjacency <- abs(adjacency)
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
cluster.all <- build.clustering(mantelCorMatrix, 5, minModuleSize)
cluster.prog <- build.clustering(mantelCorMatrixProg, 7, minModuleSize)
cluster.nonprog <- build.clustering(mantelCorMatrixNonProg, 7, minModuleSize)


table(colorsCoSplicEx)
length(table(colorsCoSplicEx)) # 113 modules 

tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
the.modules.size <- table(colorsCoSplicEx)
module.summary <- tableGrob(t(cbind(as.data.frame(the.modules.size)[1:15,], 
                                    as.data.frame(the.modules.size)[16:30,],
                                    as.data.frame(the.modules.size)[31:45,],
                                    as.data.frame(the.modules.size)[46:60,],
                                    as.data.frame(the.modules.size)[61:75,],
                                    as.data.frame(the.modules.size)[76:90,], 
                                    as.data.frame(the.modules.size)[91:105,], 
                                    as.data.frame(the.modules.size)[106:120,])), theme = tt)
grid.arrange(module.summary)

modulesCoSplicEx = names(table(colorsCoSplicEx))
sum(colorsCoSplicEx == "grey")
modulesCoSplicEx = names(table(colorsCoSplicEx))

# save(adjCoSplicEx, colorsCoSplicEx, modulesCoSplicEx, gene.names, file = "co-spliceEx/output_data/adjModulesCoSplicEx.RData")
# load("co-spliceEx/output_data/adjModulesCoSplicEx.RData")

names(colorsCoSplicEx) = gene.names
modulesCoSplicEx = names(table(colorsCoSplicEx))
try(dir.create("co-spliceEx/moduleGeneList"), silent = T)

CoSplicExConn = intramodularConnectivity(adjCoSplicEx, colorsCoSplicEx, scaleByMax=T)
totalScaledConnectivity = CoSplicExConn[,"kTotal"]/max(CoSplicExConn[,"kTotal"])
CoSplicExConn = cbind(CoSplicExConn, totalScaledConnectivity)

# -----------------------------------------------------------------------------------
# For each node (gene), the function sums adjacency entries (excluding the diagonal) 
# to other nodes within the same module. Optionally, the connectivities can be scaled 
# by the maximum connectivy in each module.
# -----------------------------------------------------------------------------------
CoSplicExConnProg <- intramodularConnectivity(mantelCorMatrixProg,  colorsCoSplicEx, scaleByMax = T)
CoSplicExConnNonProg <- intramodularConnectivity(mantelCorMatrixNonProg,  colorsCoSplicEx, scaleByMax = T)

for (module in modulesCoSplicEx){
  # print(module)
  currModuleInfo = cbind(rownames(CoSplicExConn)[colorsCoSplicEx == module] , CoSplicExConn[colorsCoSplicEx == module,"kWithin"])
  write.csv(currModuleInfo, file = paste("co-spliceEx/moduleGeneList/module_", module, ".csv", sep = ""), row.names = F, col.names = F)  
}

try(dir.create("co-spliceEx/moduleGeneList_txt"), silent = T)
for(module in modulesCoSplicEx){
  write.table(rownames(CoSplicExConn)[colorsCoSplicEx == module] , file = paste("co-spliceEx/moduleGeneList_txt/", module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}

# with out 7 samples -----------------------------------------------------
save(adjCoSplicEx_Prog, adjCoSplicEx_NonProg, colorsCoSplicEx_Prog, colorsCoSplicEx_NonProg, file = "co-spliceEx/output_data/condition_clustering_reduced_testcanberra.RData")
save(adjCoSplicEx, colorsCoSplicEx, file = "co-spliceEx/output_data/alldata_similarityAdjacency_reduced.RData")
save(hierADJCoSplicEx_Prog, hierADJCoSplicEx_NonProg, hybridCoSplicEx_Prog, hybridCoSplicEx_NonProg, file = "co-spliceEx/output_data/conditions_dynamicClustering_reduced.RData")
save(hierADJCoSplicEx, hybridCoSplicEx, file = "co-spliceEx/output_data/alldata_dynamicClusteringreduced.RData")

load("co-spliceEx/output_data/condition_clustering_reduced.RData")
load("co-spliceEx/output_data/alldata_similarityAdjacency_reduced.RData")
load("co-spliceEx/output_data/condition_clustering_reduced_testcanberra.RData")
# with 7 samples ---------------------------------------------------------
# save(adjCoSplicEx_Prog, adjCoSplicEx_NonProg, colorsCoSplicEx_Prog, colorsCoSplicEx_NonProg, file = "co-spliceEx/output_data/condition_clustering.RData")
# save(adjCoSplicEx_Prog, adjCoSplicEx_NonProg, colorsCoSplicEx_Prog, colorsCoSplicEx_NonProg, file = "co-spliceEx/output_data/condition_similarityAdjacency.RData")
# save(adjCoSplicEx, colorsCoSplicEx, file = "co-spliceEx/output_data/alldata_similarityAdjacency.RData")
# save(hierADJCoSplicEx_Prog, hierADJCoSplicEx_NonProg, hybridCoSplicEx_Prog, hybridCoSplicEx_NonProg, file = "co-spliceEx/output_data/conditions_dynamicClustering.RData")
# save(hierADJCoSplicEx, hybridCoSplicEx, file = "co-spliceEx/output_data/alldata_dynamicClustering.RData")

load("co-spliceEx/output_data/condition_similarityAdjacency.RData")
load("co-spliceEx/output_data/conditions_dynamicClustering.RData")
load("co-spliceEx/output_data/alldata_similarityAdjacency.RData")
load("co-spliceEx/output_data/alldata_dynamicClustering.RData")
# ---------------------------------------------------------------------------
# plot dendrograms and color bands 
# ---------------------------------------------------------------------------
# save(cluster.all, file = "co-spliceEx/output_data/clustering_reduced_all.RData")
load("co-spliceEx/output_data/clustering_reduced_all.RData")
hierADJCoSplicEx <- cluster.all$geneTree
colorsCoSplicEx  <- cluster.all$dynamicColors
adjCoSplicEx <- cluster.all$adj

hierADJCoSplicEx_Prog <- cluster.prog$geneTree
colorsCoSplicEx_Prog  <- cluster.prog$dynamicColors
adjCoSplicEx_Prog <- cluster.prog$adj

hierADJCoSplicEx_NonProg <- cluster.nonprog$geneTree
colorsCoSplicEx_NonProg  <- cluster.nonprog$dynamicColors
adjCoSplicEx_NonProg <- cluster.nonprog$adj

save(adjCoSplicEx_Prog, adjCoSplicEx_NonProg, colorsCoSplicEx_Prog, colorsCoSplicEx_NonProg, file = "co-spliceEx/output_data/condition_updated_clustering_reduced.RData")
load("co-spliceEx/condition_updated_clustering_reduced.RData")

sizeGrWindow(12, 9)
pdf(file = "co-spliceEx/plots/dendro_cospliceEx_reduced_TOMdynamicClustering.pdf", wi = 12, he = 9)
plotDendroAndColors(hierADJCoSplicEx, colorsCoSplicEx,
                    "Dynamic Tree Cut",
                    main = "Clustering of cospliceEx of both conditions",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(hierADJCoSplicEx_Prog, colorsCoSplicEx_Prog,
                    "Dynamic Tree Cut",
                    main = "Clustering of cospliceEx of Progressors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(hierADJCoSplicEx_NonProg, colorsCoSplicEx_NonProg,
                    "Dynamic Tree Cut",
                    main = "Clustering of cospliceEx of Non-Progressors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# -----------------------------------------------------------------------------------
# co-splice module preservation of conditions on exacloud 
# -----------------------------------------------------------------------------------
multiData <- vector("list", 2)

multiData[[1]] <- list(data = adjCoSplicEx_Prog)
multiData[[2]] <- list(data = adjCoSplicEx_NonProg)

names(multiData) <- c("Progressor", "nonProgressor", "All")
checkSets(multiData, checkStructure = FALSE, useSets = NULL)

multiColor <- vector("list", 2)

multiColor[[1]] <- as.vector(colorsCoSplicEx_Prog)
multiColor[[2]] <- as.vector(colorsCoSplicEx_NonProg)

names(multiColor) <- c("Progressor", "nonProgressor")

mp <- modulePreservation(
  multiData,
  multiColor,
  referenceNetworks = 1,
  dataIsExpr = F,
  networkType = "unsigned",
  corFnc = "bicor",
  corOptions = "use = 'p'",
  nPermutations = 200,
  includekMEallInSummary = FALSE,
  restrictSummaryForGeneralNetworks = FALSE,
  calculateQvalue = FALSE,
  randomSeed = 45,
  maxGoldModuleSize = 1000,
  maxModuleSize = 1000,
  quickCor = 1,
  ccTupletSize = 2,
  calculateCor.kIMall = TRUE,
  useInterpolation = FALSE,
  checkData = F,
  greyName = "grey",
  savePermutedStatistics = FALSE,
  loadPermutedStatistics = FALSE,
  plotInterpolation = FALSE,
  discardInvalidOutput = TRUE,
  verbose = 3, indent = 0)

#save(mp, file = "co-spliceEx/co-spliceExModulePreservation_reduced_bicor_updated2.RData")
load("co-spliceEx/co-spliceExModulePreservation_reduced_bicor_updated2.RData")

# --------------------------------------------------------------------------
ref <- 1
test <- 2
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])

# Compare preservation to quality -------------------------------------------
compare.preservation.scores <- cbind(rownames(mp$preservation$observed[[ref]][[test]]), 
                                     statsObs[, c("medianRank.pres", "medianRank.qual")],
                                     signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) 

openxlsx::write.xlsx(cbind(rownames(mp$preservation$observed[[ref]][[test]]), 
                           mp$preservation$observed[[ref]][[test]]$moduleSize,
                           statsObs[, c("medianRank.pres", "medianRank.qual")],
                           signif(statsZ, 2)), 
                     file = "co-spliceEx/preservation_AllZscores_unsigned_reduced_bicor_testcanberra_pwr7.xlsx")

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
pdf(file  = "co-spliceEx/plots/modulePreservation_reduced_bicor.pdf", wi = 10, h = 5)
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
       xlim = c(10, 1000), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08)
  # For Zsummary, add threshold lines between 2 (low) and 10 (high) ----------------------
  if (p == 2) {
    abline(h = 0)
    abline(h = 2, col = "darkblue", lty = 2)
    abline(h = 10, col = "darkgreen", lty = 2)
  }
}
dev.off()

# extract preservation modules results gene name and scaled heatmaps ----------------
try(dir.create("co-spliceEx/preservationModuleGeneList_txt_reduced_bicor_PWR7"), silent = T)
for(module in unique(modColors)){
  write.table(rownames(adjCoSplicEx_Prog)[colorsCoSplicEx_Prog == module] , file = paste("co-spliceEx/preservationModuleGeneList_txt_reduced_bicor_PWR7/", module, ".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
}