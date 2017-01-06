# set working directory ------------------------------------------------------
# setwd(getwd())
# ----------------------------------------------------------------------------
usePackage <- function(p)
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

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

options(stringsAsFactors = FALSE)

cosplice.lightyellow <- read.delim(file = "co-spliceEx/preservationModuleGeneList_txt_reduced_bicor_PWR7/lightyellow.txt", header = F, sep = "\n")
cosplice.purple <- read.delim(file = "co-spliceEx/preservationModuleGeneList_txt_reduced_bicor_PWR7/purple.txt", header = F, sep = "\n")


coexpr.black <- read.delim(file = "co-expression_3/consens_module_genes/black.txt", header = F, sep = "\n")
coexpr.cyan <- read.delim(file = "co-expression_3/consens_module_genes/cyan.txt", header = F, sep = "\n")
coexpr.blue <- read.delim(file = "co-expression_3/consens_module_genes/blue.txt", header = F, sep = "\n")
coexpr.tan <- read.delim(file = "co-expression_3/consens_module_genes/tan.txt", header = F, sep = "\n")
coexpr.yellow <- read.delim(file = "co-expression_3/consens_module_genes/yellow.txt", header = F, sep = "\n")
coexpr.grey60 <-  read.delim(file = "co-expression_3/consens_module_genes/grey60.txt", header = F, sep = "\n")
coexpr.lightgreen <- read.delim(file = "co-expression_3/consens_module_genes/lightgreen.txt", header = F, sep = "\n")
coexpr.turquoise <- read.delim(file = "co-expression_3/consens_module_genes/turquoise.txt", header = F, sep = "\n")
coexpr.purple <- read.delim(file = "co-expression_3/consens_module_genes/purple.txt", header = F, sep = "\n")
coexpr.pink <- read.delim(file = "co-expression_3/consens_module_genes/pink.txt", header = F, sep = "\n")
coexpr.lightcyan <- read.delim(file = "co-expression_3/consens_module_genes/lightcyan.txt", header = F, sep = "\n")
coexpr.salmon <- read.delim(file = "co-expression_3/consens_module_genes/salmon.txt", header = F, sep = "\n")

test <- list(coexpr.black$V1, coexpr.cyan$V1, coexpr.blue$V1, coexpr.tan$V1, coexpr.yellow$V1, coexpr.grey60$V1, 
     coexpr.lightgreen$V1, coexpr.turquoise$V1, coexpr.purple$V1, coexpr.pink$V1, coexpr.lightcyan$V1, coexpr.salmon$V1)
names(test) <- c("black", "cyan", "blue", "tan", "yellow", "grey60", "lightgreen", "turquoise", "purple", "pink", "lightcyan", "salmon")
lapply(seq_along(test), function(j)lapply(seq_along(test), function(i) if(length(which(test[[j]] %in% test[[i]])) > 1){print(cbind(i,j))}))

test2 <- list(cosplice.lightyellow, cosplice.purple)
names(test2) <- c("lightyellow", "cosplice.purple")
lapply(seq_along(test), function(j)lapply(seq_along(test2), function(i) if(length(which(test[[j]] %in% test2[[i]])) > 1){print(cbind(i,j))}))
lapply(seq_along(test2), function(j)lapply(seq_along(test2), function(i) if(length(which(test2[[j]] %in% test2[[i]])) > 1){print(cbind(i,j))}))

# -----------------------------------------------------------------------
# splicing sig 
# -----------------------------------------------------------------------
#load("co-spliceEx/output_data/splicingSignifResults_WOlfactory_reduced.RData")
#load("co-spliceEx/splicingSignif_NonProg.RData")
 load("co-spliceEx/splicingSignif_Prog.RData")
splicing.alc.purple <- geneSplicingSignifAlcohol[which(rownames(geneSplicingSignifAlcohol) %in% cosplice.purple$V1),]
splicing.smoke.purple <- geneSplicingSignifSmoking[which(rownames(geneSplicingSignifSmoking) %in% cosplice.purple$V1),]

splicing.alc.lightyellow <- geneSplicingSignifAlcohol[which(rownames(geneSplicingSignifAlcohol) %in% cosplice.lightyellow$V1),]
splicing.smoke.lightyellow <- geneSplicingSignifSmoking[which(rownames(geneSplicingSignifSmoking) %in% cosplice.lightyellow$V1),]

plot.ssig <- function(dat, xtt, ytt, tt){
melted.dat <- melt(dat)
ggplot(melted.dat, aes(Var2, value, colour = Var2)) + 
  geom_line() + 
  labs(colour = "stats") +
  xlab(xtt) + 
  ylab(ytt) + 
  ggtitle(tt)
}

tt <- ttheme_default(colhead = list(fg_params = list(parse = TRUE)))
purple.ssig <- tableGrob(splicing.alc.purple[which(splicing.alc.purple[,1] > 0.2),], theme = tt)
lightyellow.ssig <- tableGrob(splicing.alc.lightyellow[which(splicing.alc.lightyellow[,1] > 0.2),], theme = tt)
grid.arrange(purple.ssig)
grid.arrange(lightyellow.ssig)


pdf(file = "co-spliceEx/plots/ssig_nonprog.pdf",  wi = 10, he = 8)
plot.ssig(splicing.alc.purple, "", "value", "purple module alcohol per day splicing sig")
plot.ssig(splicing.smoke.purple, "", "value", "purple module pack years splicing sig")
grid.arrange(purple.ssig)
plot.ssig(splicing.alc.lightyellow, "", "value", "purple module alcohol per day splicing sig")
plot.ssig(splicing.smoke.lightyellow, "", "value", "purple module pack years splicing sig")
grid.arrange(lightyellow.ssig)
dev.off()

psa <- splicing.alc.purple[which(splicing.alc.purple[,1] > 0.2),]
pss <- splicing.smoke.purple[which(splicing.smoke.purple[,1] > 0.2),]
lsa <- splicing.alc.lightyellow[which(splicing.alc.lightyellow[,1] > 0.2),]
lss <- splicing.smoke.lightyellow[which(splicing.smoke.lightyellow[,1] > 0.2),]

the.list.purple <- list(splicing.alc.purple, splicing.smoke.purple, psa, pss)
names(the.list.purple) <- c("drink_per_day", "pack_years", "drink_per_day_sigOver0.2", "pack_years_sigOver0.2")
openxlsx::write.xlsx(the.list.purple , rowNames = T, file = "co-spliceEx/ssig_nonpreserved/purpleGenes_prog.xlsx")
the.list.lightyellow <- list(splicing.alc.lightyellow, splicing.smoke.lightyellow, lsa, lss)
names(the.list.lightyellow) <- c("drink_per_day", "pack_years", "drink_per_day_sigOver0.2", "pack_years_sigOver0.2")
openxlsx::write.xlsx(the.list.lightyellow , rowNames = T, file = "co-spliceEx/ssig_nonpreserved/lightyellowGenes_prog.xlsx")

lapply(seq_along(the.list.purple), function(i){
  write.table(rownames(the.list.purple[[i]]) , file = paste("co-spliceEx/ssig_nonpreserved/purple_nonprog/",names(the.list.purple)[[i]] ,".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
})

lapply(seq_along(the.list.lightyellow), function(i){
  write.table(rownames(the.list.lightyellow[[i]]) , file = paste("co-spliceEx/ssig_nonpreserved/lightyellow_nonprog/",names(the.list.lightyellow)[[i]] ,".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
})

lapply(seq_along(the.list.purple), function(i){
  write.table(rownames(the.list.purple[[i]]) , file = paste("co-spliceEx/ssig_nonpreserved/purple_prog/",names(the.list.purple)[[i]] ,".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
})

lapply(seq_along(the.list.lightyellow), function(i){
  write.table(rownames(the.list.lightyellow[[i]]) , file = paste("co-spliceEx/ssig_nonpreserved/lightyellow_prog/",names(the.list.lightyellow)[[i]] ,".txt", sep = ""), quote = F, sep = "\n", append = F, row.names = F, col.names = F)
})

