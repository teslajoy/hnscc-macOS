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
library(plotly)
library(survival)
library(dplyr)

options(stringsAsFactors = FALSE)
all.exon <- read.csv("Data/exonall_matrix.txt", sep = "\t", header = T, check.names = FALSE)

# gene annotation file from sample and data relationship file SDRF
gaf <- read.delim("/Users/nasim/Downloads/gabby_hnscc/hnsc_august_rnaseq/RNASeqV2/TCGAhg19June2011gaf.txt")
all.genes <- read.delim(file = "output/allgenes_HC_DE.txt")
unique(gaf$FeatureType)
gaf.genes <- gaf[gaf$FeatureType == "gene",]
gaf.genes$FeatureID <- gsub("\\|.*", "", gaf.genes$FeatureID)

#gaf.genes <- gaf.genes[gaf.genes$FeatureID %in% colnames(all.genes),]

chr <- unlist(lapply(gaf.genes$CompositeCoordinates, function(x) gsub(":|_", "", substr(x, 4 , 5))))
start <- unlist(lapply(gaf.genes$CompositeCoordinates, function(x) gsub("^.*\\:", "", gsub("\\-.*", "", x))))
end <- unlist(lapply(gaf.genes$CompositeCoordinates, function(x) gsub("\\:.*", "", gsub("^.*\\-", "", substr(x, 1, nchar(x) - 1)))))
upordown <- unlist(lapply(gaf.genes$CompositeCoordinates, function(x) substr(x, nchar(x) - 1 + 1, nchar(x))))
genecoord <- as.data.frame(cbind(gaf.genes$FeatureID, chr, start, end, upordown), stringsAsfactors = F)
names(genecoord)[1] <- "hgcn"
gaf.genes$CompositeCoordinates[[1]]

cords <- lapply(gaf.genes$CompositeCoordinates, function(x){
  read.table(text = substr(gsub(":", "" , gsub("^.*?:", "", x)), 1, nchar(gsub(":", "" , gsub("^.*?:", "", x))) - 1), 
             sep = ",", strip.white = T, header = F, fill = F, flush = T, stringsAsFactors = F, 
             colClasses = "character")})
names(cords) <- gaf.genes$FeatureID


load("/Users/nasim/Documents/thesis_prelim/HNSC/R/data/exon_prog_files.RData")

chr <- unlist(lapply(exon.prog.files[[1]][,1], function(x) gsub(":|_", "", substr(x, 4 , 5))))
start <- unlist(lapply(exon.prog.files[[1]][,1], function(x) gsub("^.*\\:", "", gsub("\\-.*", "", x))))
end <- unlist(lapply(exon.prog.files[[1]][,1], function(x) gsub("\\:.*", "", gsub("^.*\\-", "", substr(x, 1, nchar(x) - 1)))))
upordown <- unlist(lapply(exon.prog.files[[1]][,1], function(x) substr(x, nchar(x) - 1 + 1, nchar(x))))
sample.coord <- cbind(chr, start, end, upordown)
sample.coord <- as.data.frame(sample.coord, stringsAsfactors = F)
sample.coord <- sample.coord[order(sample.coord$chr), ]
sample.coord$hgnc <- NA 


hgnc <- lapply(seq_along(exon.prog.files[[1]][,1]), function(i){
  
  g <- genecoord$hgcn[which(genecoord$upordown == sample.coord$upordown[i] & genecoord$chr == sample.coord$chr[i] 
                            & as.numeric(sample.coord$start[i]) >= as.numeric(genecoord$start) & as.numeric(sample.coord$end[i]) <= as.numeric(genecoord$end))]
})

hgnc[lapply(hgnc,length) == 0] <- "TEMPORARY"
f <- do.call(cbind, hgnc)
# length(f[1, ])
mapped.names <- cbind(f[1, ], exon.prog.files[[1]][,1])
colnames(mapped.names) <- c("hgnc", "chrlocus")

write.table(mapped.names, file = "output/hgnc_mapped_chrlocus.txt", col.names = T, row.names = F, sep = "\t", append = F, quote = F)

# length(colnames(all.exon))
# dim(mapped.names)

mapped.names <- as.data.frame(mapped.names, stringsAsfactors = F)
all.exon.locus <- as.data.frame(colnames(all.exon), stringsAsfactors = F)
names(all.exon.locus) <- "chrlocus"
merged.names <- merge(mapped.names, all.exon.locus, by = "chrlocus")

library("dplyr")
merged.names <- dplyr::left_join(mapped.names, all.exon.locus, by = "chrlocus")
head(merged.names)
head(all.exon.locus)
all.exon <- t(all.exon) 
all.exon <- cbind(merged.names, all.exon)
mp.names <- all.exon[, 1:2]
mp.names <- as.data.frame(mp.names, stringsAsfactors = F)
write.table(mp.names, file = "output/final_hgnc_mapped_chrlocus.txt", col.names = T, row.names = F, sep = "\t", append = F, quote = F)

## extract genes in pipeline 
exon.in.gene <- mp.names[ mp.names$hgnc %in% colnames(all.genes), ] 
dim(exon.in.gene)
length(unique(exon.in.gene$hgnc))
length(colnames(all.genes))

unique(exon.in.gene$hgnc)[[1]]
all.exon$chrlocus <- NULL
# do this on data on all samples 
# reluting list will be one job per 
each.genes.exon <- lapply(seq_along(unique(exon.in.gene$hgnc)), function(i) {
dplyr::filter(all.exon, all.exon$hgnc == unique(exon.in.gene$hgnc)[[i]])
})
exon.sizes <- lapply(each.genes.exon, nrow)
#summary(unlist(exon.sizes))
logged.exons <- log2( all.exon[,3:length(colnames(all.exon))] + 1)
head(each.genes.exon)
save(each.genes.exon, file = "data/listexontogene.RData")
save(all.exon, file = "data/listexontogeneWithChrLocation.RData")
save(logged.exons, file = "data/listexontogeneWithChrLocationandLogged.RData")

adds <- apply(logged.exons, 1, sum)
adds <- as.data.frame(adds, stringsAsfactors = F)
colnames(adds) <- "sum"
zero.reads.exonnames <- rownames(adds)[which(adds == 0)]
length(zero.reads.exonnames) #1881
