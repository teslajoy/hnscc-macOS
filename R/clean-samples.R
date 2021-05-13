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

library(WGCNA)
library(openxlsx)
library(readxl)
library(flashClust)
library(cluster)
library(affy)
library(reshape2)
library(ggplot2)
library(geneplotter) 
library(RColorBrewer)
library(genefilter)

#---------------------------------------------------------------
options(stringsAsFactors = FALSE)
#---------------------------------------------------------------

# Read normalized_counts of gene expression --------------------
gene.prog <- read.delim(file = "Data/geneprog_matrix.txt")
gene.nonprog <- read.delim(file = "Data/genenonprog_matrix.txt")
# Read DE genes from published paper ---------------------------
# [Available At]: http://www.ncbi.nlm.nih.gov/pubmed/26747525 
DEgenes.paper <- readxl::read_excel(path = "Data/DEgenes_paper.xlsx", sheet = 1)

# make a local copy of original --------------------------------
gnp <- gene.nonprog 
gp <- gene.prog

# EDA; combine all data  ----------------------------------------
all.genes <- rbind(gene.prog, gene.nonprog)

#----------------------------------------------------------------
# log transform data based on advised suggestion and common practice 
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
# https://www.biostars.org/p/106127/
# "One could also start with normalized counts (or RPKM/FPKM data) 
# and log-transform them using log2(x+1)."
#----------------------------------------------------------------
all.genes <- log2(all.genes + 1)

# EDA; remove sum of zero for a gene's expression value --------
adds <- apply(all.genes, 2, sum) # 2 indicates columns
adds <- as.data.frame(adds, stringsAsfactors = F)
colnames(adds) <- "sum"

zero.reads.genenames <- rownames(adds)[which(adds == 0)]

adds <- as.data.frame(cbind(rownames(adds)[which(!rownames(adds) %in% zero.reads.genenames)], 
      adds[which(!rownames(adds) %in% zero.reads.genenames), ]), stringsAsfactors = F)

olfactoryrecept <- colnames(all.genes)[grep("OR", substr(colnames(all.genes), 1, 2))]
length(olfactoryrecept) # 395
all.genes <- all.genes[ ,which(!colnames(all.genes) %in% union(zero.reads.genenames, olfactoryrecept))]

dim(all.genes)
# Export: which genes had zero reads?  --------------------------
write.table(union(zero.reads.genenames, olfactoryrecept) , file = "output_2/zeroreadsAndOR.txt", quote = F, sep = "\t", append = F, row.names = F, col.names = F)

#----------------------------------------------------------------
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/HumanBrainTranscriptome/Identification%20and%20Removal%20of%20Outlier%20Samples%20-%20Illumina.pdf
#----------------------------------------------------------------
pdf(file = "output_2/IACandBOXPLOT.pdf", wi = 8, he = 7)
IAC <- cor(t(all.genes), use = "p")
hist(IAC, col = "lightblue", 
     main = paste("Histogram of inter-array correlation (IAC)  mean = ", 
                  format(mean(IAC[upper.tri(IAC)]), digits = 3) , 
                  "\n All HNSCC gene-expression HNSC IlluminaHiSeq RNASeqV2"))
boxplot(t(all.genes), las = 2, cex = 0.7, 
        main = "Overall distribution of All HNSCC \n gene-expression data", col = "lightblue") 
dev.off()  

# -------------------------------------------------------------
# Choose soft threshold 
# -------------------------------------------------------------
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# if corFunc = "bicor",then corOptions = "use = 'p', maxPOutliers = 0.1"
# blockSize = n is dependent on machines operating system 
gene.all.sft <- pickSoftThreshold(all.genes, powerVector = powers, verbose = 5, 
                                  networkType = "unsigned", corFnc = "bicor", blockSize = 10024,
                                  moreNetworkConcepts = TRUE)
system(command = "say job done!")
gene.all.sft$powerEstimate #5
#------------------------------------------------------------------------
# extract hubs based on Sum of node-connectivity quantile 
#------------------------------------------------------------------------
extract.hub.genes <- function(gene.dat, sft){
  # gene.dat format: sample X genes ---------------------------
  connCounts <- softConnectivity(gene.dat, power = sft,
                                 corFnc = "bicor", corOptions = "use = 'p'", 
                                 type = "unsigned",
                                 blockSize = 10024, 
                                 minNSamples = NULL, 
                                 verbose = 2, indent = 0)
  # extract quantiles  ----------------------------------------
  quantileConn <- quantile(connCounts, seq(0, 1, 0.1))  
  # 50% or higher quantile ------------------------------------
  geneReadsHighConn <- gene.dat[ ,connCounts > quantileConn[6]]
  geneNamesHighConn <- colnames(geneReadsHighConn)
  
  return(list("node.connectivity" = connCounts, 
              "quantile.connectivity" = quantileConn ,
              "gene.highconn.names" = geneNamesHighConn, 
              "gene.highconn.reads" = geneReadsHighConn))
}
gene.all.hc <- extract.hub.genes(all.genes, 5)

dim(gene.all.hc$gene.highconn.reads)

ag <- all.genes
all.genes <- gene.all.hc$gene.highconn.reads
missed.genes <- DEgenes.paper$DEgenes[which(!DEgenes.paper$DEgenes %in% colnames(all.genes))]
paper.gene <- ag[,which(colnames(ag) %in% missed.genes)]
all.genes <- cbind(all.genes, paper.gene)

# length(unique(colnames(all.genes)))
# length(colnames(all.genes))
# -------------------------------------------------------------
# All genes 
# Check existance of low reads in samples
# WGCNA functions won't run if gsg$allOK returns FALSE 
# -------------------------------------------------------------
gsg <- goodSamplesGenes(all.genes, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Print the gene and sample names that were removed ----------
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(all.genes)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(genes.prog)[!gsg$goodSamples], collapse = ", ")))
  # Remove the outlier genes and samples from the data ---------
  all.genes <- all.genes[gsg$goodSamples, gsg$goodGenes]
}

# -------------------------------------------------------------
# Segment Data/choose HC genes by condition 
# -------------------------------------------------------------
genes.prog <- log2(gene.prog[, colnames(all.genes)] + 1)
genes.nonprog <- log2(gene.nonprog[, colnames(all.genes)] + 1)
# -------------------------------------------------------------
# Progressor
# Check existance of low reads in samples
# WGCNA functions won't run if gsg$allOK returns FALSE 
# -------------------------------------------------------------
gsg <- goodSamplesGenes(genes.prog, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Print the gene and sample names that were removed ----------
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(genes.prog)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(genes.prog)[!gsg$goodSamples], collapse = ", ")))
  # Remove the outlier genes and samples from the data ---------
  genes.prog <- genes.prog[gsg$goodSamples, gsg$goodGenes]
}

# re-check 
gsg <- goodSamplesGenes(genes.prog, verbose = 3)
gsg$allOK

# -------------------------------------------------------------
# Non-progressed 
# Check existance of low reads in samples
# WGCNA functions won't run if gsg$allOK returns FALSE 
# -------------------------------------------------------------
gsg <- goodSamplesGenes(genes.nonprog, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Print the gene and sample names that were removed ----------
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(genes.nonprog)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(genes.nonprog)[!gsg$goodSamples], collapse = ", ")));
  # Remove the outlier genes and samples from the data ---------
  genes.nonprog <- genes.nonprog[gsg$goodSamples, gsg$goodGenes]
}
# re-check 
gsg <- goodSamplesGenes(genes.nonprog, verbose = 3)
gsg$allOK

genes.prog <- genes.prog[,intersect(colnames(genes.prog), colnames(genes.nonprog))]
genes.nonprog <- genes.nonprog[,intersect(colnames(genes.prog), colnames(genes.nonprog))]
all.genes <- all.genes[,which(colnames(all.genes) %in% intersect(colnames(genes.prog), colnames(genes.nonprog)))]

# all paper genes made it thought ------------------------------
DEgenes.paper$DEgenes %in% colnames(genes.prog) 
DEgenes.paper$DEgenes %in% colnames(genes.nonprog) 

write.table(genes.prog , file = "output_2/genesprog_HC.txt", sep = "\t", append = F)
write.table(genes.nonprog , file = "output_2/genesnonprog_HC.txt", sep = "\t", append = F)
write.table(all.genes , file = "output_2/allgenes_HC.txt", sep = "\t", append = F)

dim(all.genes) # 229 sample X 10024 final gene count selected
# ----------------------------------------------------------------
# Calculate sftThresh for conditions 
# ----------------------------------------------------------------
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
# if corFunc = "bicor",then corOptions = "use = 'p', maxPOutliers = 0.1"
gene.prog.sft <- pickSoftThreshold(gene.prog, powerVector = powers, 
                                   verbose = 5, networkType = "unsigned", blockSize = 10024,
                                   corFnc = "cor", moreNetworkConcepts = TRUE)

gene.nonprog.sft <- pickSoftThreshold(gene.nonprog, powerVector = powers, 
                                      verbose = 5, networkType = "unsigned", blockSize = 10024,
                                      corFnc = "cor", moreNetworkConcepts = TRUE)
system(command = "say job done!")
gene.prog.sft$powerEstimate    #5
gene.nonprog.sft$powerEstimate #6

# ----------------------------------------------------------------
# plot Sample clustering to detect outliers dendrogram 
# ----------------------------------------------------------------
sampleTree.nonprog <- hclust(dist(gene.nonprog), method = "average")
sampleTree.prog <- hclust(dist(gene.prog), method = "average")
sampleTree.all <- hclust(dist(all.genes), method = "average")

sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))

pdf(file = "output_2/sample_clustering.pdf",  wi = 12, he = 9)
  plot(sampleTree.all, 
       main = "Sample clustering to detect sample outliers in all HNSCC samples", 
       sub = "", xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  plot(sampleTree.prog, 
       main = "Sample clustering to detect sample outliers in HNSCC Progressors", 
       sub = "", xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  plot(sampleTree.nonprog, 
       main = "Sample clustering to detect sample outliers in HNSCC Non-Progressors", 
       sub = "", xlab = "", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
dev.off()

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

# pdf(file = "output_2/coexpr_softthr.pdf", wi = 9, he = 6)
#   plot.soft.thr(gene.prog.sft, "Progressed HC")
#   plot.soft.thr(gene.nonprog.sft, "Non-Progressed HC")
#   plot.soft.thr(gene.all.sft, "All samples HC")
# dev.off()

pdf(file = "co-expression_3/plots/sftThresh.pdf", wi = 12, he = 9)
plot.soft.thr(gene.all.sft, "co-expression All Data")
plot.soft.thr(gene.prog.sft, "co-expression Progressor")
plot.soft.thr(gene.nonprog.sft, "co-expression NonProgressor" )
dev.off()
# -------------------------------------------------------------------------------------
# IAC of final results 
# -------------------------------------------------------------------------------------
pdf(file = "output_2/final_IAC.pdf", wi = 9, he = 6)
IAC <- cor(t(all.genes), use = "p")
hist(IAC, col = "lightblue", 
     main = paste("Histogram of inter-array correlation (IAC)  mean = ", 
                  format(mean(IAC[upper.tri(IAC)]), digits = 3) , 
                  "\n All HNSCC gene-expression HNSC IlluminaHiSeq RNASeqV2"))
boxplot(t(all.genes), las = 2, cex = 0.7, 
        main = "Overall distribution of All HNSCC \n gene-expression data", col = "lightblue") 

IAC <- cor(t(genes.prog), use = "p")
hist(IAC, col = "lightblue", 
     main = paste("Histogram of inter-array correlation (IAC)  mean = ", 
                  format(mean(IAC[upper.tri(IAC)]), digits = 3) , 
                  "\n Progressor HNSCC gene-expression HNSC IlluminaHiSeq RNASeqV2"))
boxplot(t(genes.prog), las = 2, cex = 0.7, 
        main = "Overall distribution of All HNSCC \n gene-expression data", col = "lightblue") 

IAC <- cor(t(genes.nonprog), use = "p")
hist(IAC, col = "lightblue", 
     main = paste("Histogram of inter-array correlation (IAC)  mean = ", 
                  format(mean(IAC[upper.tri(IAC)]), digits = 3) , 
                  "\n None-Progressor HNSCC gene-expression HNSC IlluminaHiSeq RNASeqV2"))
boxplot(t(genes.nonprog), las = 2, cex = 0.7, 
        main = "Overall distribution of All HNSCC \n gene-expression data", col = "lightblue") 
dev.off()
