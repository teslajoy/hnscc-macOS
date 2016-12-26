#setwd("/Users/nasim/Documents/thesis_prelim/HNSC/R")
#http://statisticalrecipes.blogspot.com/2012/08/biomart-find-gene-name-using-chromosome.html
library(biomaRt) # map exon to gene by chromosome location 

# load data --------------------------------------------------
load("RNA_normalized_samples/gene_prog_files.RData")
load("RNA_normalized_samples/gene_nonprog_files.RData")
load("RNA_normalized_samples/gene_files.RData")
load("RNA_normalized_samples/exon_prog_files.RData")
load("RNA_normalized_samples/exon_nonprog_files.RData")
load("RNA_normalized_samples/exon_files.RData")


# 20531 genes and 239322 exons -----------------------------
unique(lapply(exon.prog.files, function(x) nrow(x)))
unique(lapply(gene.prog.files, function(x) nrow(x)))

# patient X gene  matrix ------------------------------------
# progressor ------------------------------------------------
gene.prog <- lapply(gene.prog.files ,function(x) cbind(x$normalized_count))
gene.prog <- do.call(cbind, gene.prog)
dim(gene.prog) # 20531   68
genenames <- gsub("\\|.*", "", gene.prog.files[[1]][[1]])
rownames(gene.prog) <- genenames
colnames(gene.prog) <- names(gene.prog.files)
gene.prog <- t(gene.prog)
gene.prog <- gene.prog[ ,which(colnames(gene.prog) != "?")]
#write.table(gene.prog, file = "geneprog_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)
#test <- read.delim(file = "geneprog_matrix.txt")


# nonprogressor --------------------------------------------------------
gene.nonprog <- lapply(gene.nonprog.files ,function(x) cbind(x$normalized_count))
gene.nonprog <- do.call(cbind, gene.nonprog)
dim(gene.nonprog) #20531   161
genenames <- gsub("\\|.*", "", gene.nonprog.files[[1]][[1]])
rownames(gene.nonprog) <- genenames
colnames(gene.nonprog) <- names(gene.nonprog.files)
gene.nonprog <- t(gene.nonprog)
gene.nonprog <- gene.nonprog[ ,which(colnames(gene.nonprog) != "?")]
#write.table(gene.nonprog, file = "genenonprog_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)
#test <- read.delim(file = "genenonprog_matrix.txt")

# all ----------------------------------------------------------------
gene.all <- lapply(gene.files ,function(x) cbind(x$normalized_count))
gene.all <- do.call(cbind, gene.all)
dim(gene.all) #20531   229
genenames <- gsub("\\|.*", "", gene.files[[1]][[1]])
rownames(gene.all) <- genenames
colnames(gene.all) <- names(gene.files)
gene.all <- t(gene.all)
gene.all <- gene.all[ ,which(colnames(gene.all) != "?")]
write.table(gene.all, file = "geneall_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)

# patient X exon matrix ------------------------------------
# progressor -----------------------------------------------
exon.prog <- lapply(exon.prog.files ,function(x) cbind(x$RPKM))
exon.prog <- do.call(cbind, exon.prog)
rownames(exon.prog) <- exon.prog.files[[1]][[1]]
colnames(exon.prog) <- names(exon.prog.files)
exon.prog <- t(exon.prog)
write.table(exon.prog, file = "exonprog_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)
#test <- read.delim(file = "exonprog_matrix.txt")
#dim(exon.prog) # 239322   68
#test <- do.call(cbind, list(exon.prog.files[[1]][[1]], exon.prog))
#rownames(exon.prog)[which(!rownames(exon.prog) %in% rownames(gene.prog))] # none 

# non progressor --------------------------------------------
exon.nonprog <- lapply(exon.nonprog.files ,function(x) cbind(x$RPKM))
exon.nonprog <- do.call(cbind, exon.nonprog)
dim(exon.nonprog) # 239322   162
rownames(exon.nonprog) <- exon.nonprog.files[[1]][[1]]
colnames(exon.nonprog) <- names(exon.nonprog.files)
exon.nonprog <- t(exon.nonprog)
rownames(exon.nonprog)[which(!rownames(exon.nonprog) %in% rownames(gene.nonprog))] #"TCGA-CQ-A4CG" 
exon.nonprog <- exon.nonprog[-which(!rownames(exon.nonprog) %in% rownames(gene.nonprog)),]
write.table(exon.nonprog, file = "exonnonprog_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)
#test <- read.delim(file = "exonnonprog_matrix.txt")

# all patients  --------------------------------------------
exon.all <- lapply(exon.files ,function(x) cbind(x$RPKM))
exon.all <- do.call(cbind, exon.all)
dim(exon.all) # 239322   230
rownames(exon.all) <- exon.files[[1]][[1]]
colnames(exon.all) <- names(exon.files)
exon.all <- t(exon.all)
rownames(exon.all)[which(!rownames(exon.all) %in% rownames(gene.all))] #"TCGA-CQ-A4CG" 
exon.all <- exon.all[-which(!rownames(exon.all) %in% rownames(gene.all)),]
write.table(exon.all, file = "exonall_matrix.txt", col.names = T, row.names = T, sep = "\t", append = F)
#test <- read.delim(file = "exonall_matrix.txt")
#gene.prog <- read.delim(file = "/home/users/sanati/thesis/Data/matrix/geneprog_matrix.txt")
