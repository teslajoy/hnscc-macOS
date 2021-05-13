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
# Data link on box: https://app.box.com/folder/15828236365
gene.prog <- read.delim(file = "co-expression_2/output/genesprog_HC.txt")
gene.nonprog <- read.delim(file = "co-expression_2/output/genesnonprog_HC.txt")
all.genes <- rbind(gene.prog, gene.nonprog)



# upload updated HPV status -----------------------------------------
# Data link on box: https://app.box.com/file/180612907516
hpv.dat <- read.csv(file = "Data/Oncotarget_merge.csv", header = T, sep = ",")
all.pos <- hpv.dat[which( hpv.dat$TCGA_Nature_Any_Positive == "Positive"),]
hpv.dat <- hpv.dat[which(hpv.dat$T_PATIENT_ID %in% rownames(all.genes)),]

# make 
hpv.dat$N_ISH_BCR.DCC[hpv.dat$N_ISH_BCR.DCC == 0] = "Negative"
hpv.dat$N_ISH_BCR.DCC[hpv.dat$N_ISH_BCR.DCC == 1] = "Postive"

hpv.dat$N_p16_BCR.DCC[hpv.dat$hpv.dat$N_p16_BCR.DCC == 0] = "Negative"
hpv.dat$N_p16_BCR.DCC[hpv.dat$hpv.dat$N_p16_BCR.DCC == 1] = "Postive"

hpv.dat$N_p16_Adel_Update[hpv.dat$N_p16_Adel_Update == 0] = "Negative"
hpv.dat$N_p16_Adel_Update[hpv.dat$N_p16_Adel_Update == 1] = "Postive"

my.pos <- hpv.dat$T_PATIENT_ID[which(hpv.dat$T_HPV_STATUS_ISH  == "Positive" | hpv.dat$T_HPV_STATUS_P16 == "Positive" | hpv.dat$N_ISH_BCR.DCC == "Positive" | hpv.dat$N_p16_BCR.DCC == "Positive" | hpv.dat$N_p16_Adel_Update == "Positive" | hpv.dat$N_Final_HPV_Status == "Positive" | hpv.dat$T_Final_HPV_Status == "Positive")]
my.neg <- hpv.dat$T_PATIENT_ID[which(hpv.dat$T_HPV_STATUS_ISH  == "Negative" | hpv.dat$T_HPV_STATUS_P16 == "Negative" | hpv.dat$N_ISH_BCR.DCC == "Negative" | hpv.dat$N_p16_BCR.DCC == "Negative" | hpv.dat$N_p16_Adel_Update == "Negative" | hpv.dat$N_Final_HPV_Status == "Negative" | hpv.dat$T_Final_HPV_Status == "Negative")]
my.neg <- my.neg[which(!my.neg %in% my.pos)]


# Data link on box: https://app.box.com/file/118550648057
clinical <- readxl::read_excel(path = "output/clinical_interest_reduced.xlsx", sheet = 1)
clinical$hpv_status_p16_or_ish[which(clinical$Unique_ID %in% my.neg)] <- "Positive"
clinical$hpv_status_p16_or_ish[which(clinical$Unique_ID %in% my.pos)] <- "Negative"
clinical$hpv_status_p16_or_ish[is.na(clinical$hpv_status_p16_or_ish)] <- "unknown"
# 149 positive # 35 negative # 45 unknown

hpv.pos <- clinical$Unique_ID[which(clinical$hpv_status_p16_or_ish == "Positive")]
hpv.neg <- clinical$Unique_ID[which(clinical$hpv_status_p16_or_ish == "Negative")]
hpv.unknown <- clinical$Unique_ID[which(clinical$hpv_status_p16_or_ish == "unknown")]


pos.genes <- all.genes[hpv.pos, ]
neg.genes <- all.genes[hpv.neg, ]
unknown.genes <- all.genes[hpv.unknown, ]

setLabels <- c("Positive", "Negative", "Unknown")
multiExpr <- list(Positive = list(data = pos.genes), 
                  Negative = list(data = neg.genes), 
                  Unknown  = list(data = unknown.genes))

# data link to hpvnet: https://app.box.com/file/180613643310
# system.time( hpvnet <- blockwiseConsensusModules(
#   multiExpr, 
#   power = 6, 
#   corType = "bicor",
#   minModuleSize = 30, 
#   deepSplit = 2,
#   maxBlockSize = 10024, 
#   randomSeed = 12345,
#   networkType = "unsigned",
#   TOMType = "signed",
#   TOMDenom = "min",
#   pamRespectsDendro = FALSE,
#   mergeCutHeight = 0.005, #.995 cut for module merging 
#   numericLabels = TRUE,
#   minKMEtoStay = 0,
#   saveTOMs = TRUE, 
#   verbose = 5))
# system(command = "say job done!")
# save(hpvnet, file = "co-expression_3/consensus/hpv_status_consensus.RData")


hpvMEs <- hpvnet$multiMEs[[1]]
hpvnetmoduleLabels <- hpvnet$colors
hpvnetmoduleColors <- labels2colors(hpvnetmoduleLabels)
hpvposTree <- hpvnet$dendrograms[[1]]


#######################################################################################
pc <- as.numeric(gsub("ME", "", colnames(hpvnet$multiMEs[[1]]$data)))
colnames(hpvnet$multiMEs[[1]]$data) <- labels2colors(pc)
npc <- as.numeric(gsub("ME", "", colnames(hpvnet$multiMEs[[2]]$data)))
colnames(hpvnet$multiMEs[[2]]$data) <- labels2colors(npc)
upc <- as.numeric(gsub("ME", "", colnames(hpvnet$multiMEs[[3]]$data)))
colnames(hpvnet$multiMEs[[3]]$data) <- labels2colors(upc)

progkME <- abs(signedKME(gene.prog, hpvnet$multiMEs[[1]]$data, outputColumnName = "kME."))
colnames(progkME) <- paste("kME.", colnames(hpvnet$multiMEs[[1]]$data), sep = "") 
nonprogkME <- abs(signedKME(gene.nonprog, hpvnet$multiMEs[[2]]$data, outputColumnName = "kME."))
colnames(nonprogkME) <- paste("kME.", colnames(hpvnet$multiMEs[[2]]$data), sep = "") 
unknownkME <- abs(signedKME(unknown.genes, hpvnet$multiMEs[[3]]$data, outputColumnName = "kME."))
colnames(unknownkME) <- paste("kME.", colnames(hpvnet$multiMEs[[3]]$data), sep = "") 

# extract top kme -----------------------------------------------
prog.topPoskME <- lapply(colnames(progkME), function(x){ 
  orderd.kME <- progkME[order(-progkME[x]),][x]
  return(head(orderd.kME, 20))
})
names(prog.topPoskME) <- colnames(progkME)


nonprog.topPoskME <- lapply(colnames(nonprogkME ), function(x){ 
  orderd.kME <- nonprogkME[order(-nonprogkME[x]),][x]
  return(head(orderd.kME, 20))
})
names(nonprog.topPoskME) <- colnames(nonprogkME)

unknown.topPoskME <- lapply(colnames(unknownkME), function(x){ 
  unknown.kME <- unknownkME[order(-unknownkME[x]),][x]
  return(head(unknown.kME, 20))
})
names(unknown.topPoskME) <- colnames(unknownkME)

#######################################################################################
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
rawAdj1 <- adjacency(pos.genes, power = 6, type = "unsigned") 
rawAdj2 <- adjacency(neg.genes, power = 6, type = "unsigned") 
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
  
  for (gene2 in geneNames[which(idxsVector == T)]){
    vectorEdges[gene2] <- r.test(n = n1, r12 = rawAdj1[x, gene2], r34 = rawAdj2[x, gene2], n2 = n2)$p
  }
  vectorEdges < pThreshold
})

mx <- do.call(rbind, mx)
af.edges <- length(which(mx == T))
colnames(mx)[mx[which(mx == T)]]
edge.change.rate <- af.edges / (length(colnames(pos.genes)))^2


geneChangeEdgeCount <- rowSums(mx)
geneChangeEdgeCount <- as.data.frame(geneChangeEdgeCount, stringsAsfacors = F)
rownames(geneChangeEdgeCount) <- colnames(pos.genes)

pValuesEdgeChange <- rep(1, length(colnames(pos.genes)))
pValuesEdgeChange <- lapply(colnames(pos.genes), function(gene){
  binom.test(x = geneChangeEdgeCount[gene,], n = length(colnames(pos.genes)), 
             p = edge.change.rate, alternative = "g")$p.value
})
names(pValuesEdgeChange) <- colnames(pos.genes)
names(pValuesEdgeChange)[which(unlist(pValuesEdgeChange) != 1)]
binomial.pvals <- pValuesEdgeChange[which(unlist(pValuesEdgeChange) != 1  )]

adjpOut <- mt.rawp2adjp(unlist(pValuesEdgeChange), proc = "BH")
fdrEdgesChange <- adjpOut$adjp[order(adjpOut$index), 2]

genesChangedEdges <- geneNames[unlist(pValuesEdgeChange) < 0.01]
geneChangeEdgeCount[genesChangedEdges,]
mean(geneChangeEdgeCount[genesChangedEdges,])

median(geneChangeEdgeCount$geneChangeEdgeCount)
median(geneChangeEdgeCount[genesChangedEdges,])

# all.genes <- rbind(pos.genes, neg.genes)

modules.enrich.DW <- lapply(seq_along(unique(hpvnetmoduleColors)), function(i) {
  x <- factor(colnames(all.genes) %in% colnames(pos.genes[, which(hpvnetmoduleColors %in% unique(hpvnetmoduleColors)[i])]), levels = c("TRUE", "FALSE"))
  y <- factor(colnames(all.genes) %in% genesChangedEdges, levels = c("TRUE", "FALSE"))
  p.val <- fisher.test(x, y, alternative = "g")$p.value
  return(p.val)
})

names(modules.enrich.DW) <- unique(hpvnetmoduleColors)
modules.enrich.DW <- unlist(modules.enrich.DW)
DW.modules <- names(modules.enrich.DW)[modules.enrich.DW < (0.05/length(unique(hpvnetmoduleColors)))]
DW.modules <- DW.modules[which(!DW.modules %in% "grey")]
# "turquoise"    "midnightblue" "grey60" 

#######################################################################################
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
cluster.pos <- build.clustering(pos.genes, 6, minModuleSize)
cluster.neg <- build.clustering(neg.genes, 6, minModuleSize)
cluster.unknown <- build.clustering(unknown.genes, 6, minModuleSize)

#######################################################################################
# here on I change the variable names in code to plot the network, Its not coded systematically
# due to different cluster of interest and time 

neghpv.topPoskME <- nonprog.topPoskME
poshpv.topPoskME <- nonprog.topPoskME
unknown.topPoskME

neghpv.dw <- neghpv.topPoskME[c("kME.turquoise", "kME.midnightblue", "kME.grey60")]
neg.turq <- cluster.neg$adj[rownames(neghpv.dw[["kME.turquoise"]]) ,rownames(neghpv.dw[["kME.turquoise"]])]
neg.midblue <- cluster.neg$adj[rownames(neghpv.dw[["kME.midnightblue"]]) ,rownames(neghpv.dw[["kME.midnightblue"]])]
neg.gsix <- cluster.neg$adj[rownames(neghpv.dw[["kME.grey60"]]) ,rownames(neghpv.dw[["kME.grey60"]])]

poshpv.dw <- poshpv.topPoskME[c("kME.turquoise", "kME.midnightblue", "kME.grey60")]
pos.turq <- cluster.pos$adj[rownames(poshpv.dw[["kME.turquoise"]]) ,rownames(poshpv.dw[["kME.turquoise"]])]
pos.midblue <- cluster.pos$adj[rownames(poshpv.dw[["kME.midnightblue"]]) ,rownames(poshpv.dw[["kME.midnightblue"]])]
pos.gsix <- cluster.pos$adj[rownames(poshpv.dw[["kME.grey60"]]) ,rownames(poshpv.dw[["kME.grey60"]])]

unknown.dw <- unknown.topPoskME[c("kME.turquoise", "kME.midnightblue", "kME.grey60")]
unknown.turq <- cluster.unknown$adj[rownames(unknown.dw[["kME.turquoise"]]) ,rownames(unknown.dw[["kME.turquoise"]])]
unknown.midblue <- cluster.unknown$adj[rownames(unknown.dw[["kME.midnightblue"]]) ,rownames(unknown.dw[["kME.midnightblue"]])]
unknown.gsix <- cluster.unknown$adj[rownames(unknown.dw[["kME.grey60"]]) ,rownames(unknown.dw[["kME.grey60"]])]

# openxlsx::write.xlsx(poshpv.dw, rowNames = T, file = "HPVpos_topkME.xlsx")
# openxlsx::write.xlsx(neghpv.dw, rowNames = T, file = "HPVneg_topkME.xlsx")
# openxlsx::write.xlsx(unknown.dw , rowNames = T, file = "HPVunknown_topkME.xlsx")

binomial.pvals <- as.data.frame(unlist(binomial.pvals))
binomial.pvals <- as.data.frame(cbind(rownames(binomial.pvals)[which(rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.turquoise"]]) | 
                                                                       rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.midnightblue"]]) | 
                                                                       rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.grey60"]]))], 
                                      binomial.pvals[which(rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.turquoise"]]) | 
                                                             rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.midnightblue"]]) | 
                                                             rownames(binomial.pvals) %in% rownames(neghpv.dw[["kME.grey60"]])),]))
# openxlsx::write.xlsx(binomial.pvals , rowNames = F, colNames = F, file = "binomialpval_dw_hpv.xlsx")


# plot the net -------------------------------------------------------
net.obj <- function(x){
  temp <- x
  temp[temp < 0.6] <- 0
  
  library(igraph)
  net <- graph.adjacency(temp, mode = "undirected", weighted = TRUE)
  net <- simplify(net, remove.loops = TRUE)
  E(net)$weight <- edge.betweenness(net)
  
  V(net)$size <- 4*sqrt(graph.strength(net))
  return(net)
}

neg.net <- net.obj(neg.midblue)
pos.net <- net.obj(pos.midblue)
unknown.net <- net.obj(unknown.midblue)

lay.neg <- layout.circle(neg.net)
lay.pos  <- layout.circle(pos.net)
lay.unknown  <- layout.circle(unknown.net)

par(mfrow = c(1,3))
plot(pos.net, edge.width = E(pos.net)$weight, vertex.size = V(pos.net)$size, vertex.color = "lightblue", layout = lay.pos, vertex.label.color = "black", main = "HPV-Pos")
plot(neg.net, edge.width = E(neg.net)$weight, vertex.size = V(neg.net)$size, vertex.color = "lightblue", layout = lay.neg, vertex.label.color = "black", main = "HPV-Neg")
plot(unknown.net, edge.width = E(unknown.net)$weight, vertex.size = V(unknown.net)$size, vertex.color = "lightblue", layout = lay.unknown, vertex.label.color = "black", main = "HPV-Unknown")

dev.off()

length(clinical$Unique_ID[which(clinical$Progression_FINAL == "Progressor" & clinical$hpv_status_p16_or_ish == "Positive")])
length(clinical$Unique_ID[which(clinical$Progression_FINAL == "Progressor" & clinical$hpv_status_p16_or_ish == "Negative")])
length(clinical$Unique_ID[which(clinical$Progression_FINAL == "Progressor" & clinical$hpv_status_p16_or_ish == "unknown")])
# out of 68, 47 pos, 6 negative, 15 unknown

length(clinical$Unique_ID[which(clinical$Progression_FINAL == "NonProgressor" & clinical$hpv_status_p16_or_ish == "Positive")])
length(clinical$Unique_ID[which(clinical$Progression_FINAL == "NonProgressor" & clinical$hpv_status_p16_or_ish == "Negative")])
length(clinical$Unique_ID[which(clinical$Progression_FINAL == "NonProgressor" & clinical$hpv_status_p16_or_ish == "unknown")])
# 102 pos, 29 negative, 30 unknown
