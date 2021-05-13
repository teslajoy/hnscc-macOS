# set working directory --------------------------------------
# setwd("/Users/nasim/Documents/thesis_prelim/HNSC/R/")
# setwd(getwd())

# install packages as required---------------------------------
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}

# library calls----------------------------------------------
library(dplyr)
library(tidyr)
library(plyr)
library(reshape2)
library(ggplot2)
library(plotly)
library(stringr)
library(gridExtra)
library(grid)
library(gtable)
library(psych)
library(readxl)
library(RCurl)
library(rjson)


# Upload data ------------------------------------------------------------------
clinic.path <- "/Users/nasim/Documents/thesis_prelim/HNSC/Clinical/Biotab/nationwidechildrens.org_clinical_patient_hnsc.txt"
clinical.dat <- read.delim(clinic.path, header = T, sep = "\t", stringsAsFactors = F)

uuidandbar <- cbind(clinical.dat$bcr_patient_uuid, clinical.dat$bcr_patient_barcode)
uuidandbar <- uuidandbar[-c(1,2), ]

# EDA --------------------------------------------------------------------------
summary(clinical.dat)
head(clinical.dat) 
tail(clinical.dat)

# remove bcr_patient_uuid and CDE_ID fields -------------------------------------
clinical.dat <- clinical.dat[-c(1, 2), -c(1)] # remove 
dim(clinical.dat) 
names(clinical.dat)

annotations <- read.delim("annotation.txt", header = T, sep = "\t", stringsAsFactors = F)
clinical.dat <- clinical.dat[which(clinical.dat[[1]] %in% annotations[[1]]),]
all.clinical.count <- nrow(clinical.dat)

uuidandbar <- uuidandbar[which(uuidandbar[,2] %in% clinical.dat[[1]]),]

annot.clinical <- read_excel("HNSCC_Table1_Data.xlsx", col_names = TRUE, col_types = NULL, na = "", skip = 0)
annot.clinical[[1]] <- gsub(".", '-', annot.clinical[[1]], fixed = T)

# omics data ------------------------------------------------------------
#rnaseqv2.path <- "/Users/nasim/Documents/thesis_prelim/HNSC/omics/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
#rnaseqv2.path <- "/Users/nasim/Downloads/hnsc_august_rnaseq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/"
rnaseqv2.gene.path <- "/Users/nasim/Documents/thesis/gabby_hnscc/hnsc_august_rnaseq/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/"
rnaseqv2.exon.path <- "/Users/nasim/Documents/thesis/gabby_hnscc/hnsc_august_rnaseq/RNASeqV2/Level_3"

# map uuid to barcodes via TCGA site ------------------------------------
# alternative would be read mappings file -------------------------------
get.mappings <- function(rnaseqv2.path){
  rnaseqv2 <- list.files(path = rnaseqv2.path)
  
  rnaseqv2.uid <- lapply(rnaseqv2, function(x)  gsub("\\..*", "", substr(x[1], 9 , 50)))
  rnaseqv2.uid <- unlist(rnaseqv2.uid)
  rnaseqv2.uid <- unique(rnaseqv2.uid)
  
  #https://gist.github.com/agrueneberg/3842413 ----------------------------
  # Read sample UUIDs -----------------------------------------------------
  uuids <- rnaseqv2.uid 
  
  # Convert to character vector -------------------------------------------
  uuids <- as.vector(t(uuids))
  
  # Query TCGA's UUID to barcode Web Service -------------------------------
  resp <- getURL("https://tcga-data.nci.nih.gov/uuid/uuidws/mapping/json/uuid/batch", 
                 customrequest = "POST", httpheader = c("Content-Type: text/plain"), 
                 postfields = paste(uuids, collapse = ","))
  
  # Extract mappings from response -----------------------------------------
  mappings <- fromJSON(resp)$uuidMapping
  
  # Extract patient barcode from sample barcode-----------------------------
  mappings <- lapply(mappings, function(mapping) {
    mapping$barcode <- substr(mapping$barcode, 0, 12)
    return(mapping)
  })
     
  # convert list into data frame -------------------------------------------
  options(stringsAsFactors = FALSE)
  mappings <- ldply(mappings, data.frame)
  
  # which are in annotation list -------------------------------------------
  mappings <- mappings[which(mappings[,1] %in% annot.clinical[[1]]), ]
  
  return(mappings)
}

# extract full list of all patients normalized gene and exon mappings -----------------
gene.mappings <- get.mappings(rnaseqv2.gene.path)
exon.mappings <- get.mappings(rnaseqv2.exon.path)

# extract progressor and non-progressor barcodes ---------------------------------------
nonprog <- annot.clinical$bcr.patient.barcode[which(annot.clinical$Progression_FINAL == "NonProgressor")]
prog <- annot.clinical$bcr.patient.barcode[which(annot.clinical$Progression_FINAL == "Progressor")]

gene.mappings.nonprog <- gene.mappings[which(gene.mappings$barcode %in% nonprog), ]
gene.mappings.prog <- gene.mappings[which(gene.mappings$barcode %in% prog), ]

exon.mappings.nonprog <- exon.mappings[which(exon.mappings$barcode %in% nonprog), ]
exon.mappings.prog <- exon.mappings[which(exon.mappings$barcode %in% prog), ]

# Alternative map from file check -----------------------------------------------------
# there are 230 patients. 162 nonprogress and 68 progressor cases ---------------------
rnaseqv2.filename.map <- read.delim(file = "/Users/nasim/Downloads/gabby_hnscc/hnsc_august_rnaseq/FILE_SAMPLE_MAP.txt")
length(unique(rnaseqv2.filename.map$barcode.s.))
rnaseqv2.filename.map[ ,3] <- substr(rnaseqv2.filename.map[,2], 0, 12)
names(rnaseqv2.filename.map)[3] <- "simple_barcode"

r <- lapply(rnaseqv2.filename.map[ ,1] , function(x)  gsub("\\..*", "", substr(x[1], 9 , 50)))
rnaseqv2.filename.map[ ,4] <- unlist(r)
names(rnaseqv2.filename.map)[4] <- "simple_uuid"
#unique(rnaseqv2.filename.map[which(rnaseqv2.filename.map[,3] %in% annot.clinical[[1]]), 3]) #230

# save txt file for future ref ---------------------------------------------------
#m <- as.data.frame(rnaseqv2.filename.map, stringsAsfactors = F)
#write.table(m, file = "mappingcombo.txt", sep = "\t", append = F)
#m <- read.delim(file = "mappingcombo.txt")

##################################################################################
# read gene files -----------------------------------------------------------------
# all patients case ---------------------------------------------------------------
gene.rnaseqv2 <- list.files(path = rnaseqv2.gene.path)
gene.filenames <- lapply(gene.mappings[,2], function(x) grep(x, gene.rnaseqv2, value = TRUE))
gene.filenames <- lapply(gene.filenames, function(x) grep("genes.normalized", x, value = TRUE))
# one file is stored at the wrong directory by previous researcher -----------------
gene.filenames <- gene.filenames[-which(lapply(gene.filenames, function(x) length(x) < 1) == T)]
# one file is wrong and contains isoform information 
# gene.filenames[208] %in% gsub(".txt", "", rnaseqv2.filename.map[[1]])
# "unc.edu.eee78dec-0b18-48ea-a793-56dd12dc3d3f.2203467.rsem.isoforms.normalized_results"
# rnaseqv2.filename.map[which("eee78dec-0b18-48ea-a793-56dd12dc3d3f" %in% rnaseqv2.filename.map$simple_uuid), ] #TCGA-CQ-7065 # 002e6d9c-3085-41da-ae05-e8bac069c9ce
# rnaseqv2.filename.map[which(rnaseqv2.filename.map$simple_barcode == "TCGA-CQ-7065"), ]
# names(gene.files[[which(lapply(gene.files, function(x) nrow(x)) == 73599)]])
gene.files <- lapply(gene.filenames, function(x) read.delim(paste(rnaseqv2.gene.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(gene.files) <- lapply(gene.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(gene.files, file = "data/gene_files.RData")

existing.gene.file.id <- rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% unlist(unique(gene.filenames))), ] # 229 
existing.proggene.file.id <- existing.gene.file.id[existing.gene.file.id$simple_barcode %in% prog, ] # 68
existing.nonproggene.file.id <- existing.gene.file.id[existing.gene.file.id$simple_barcode %in% nonprog, ] # 161 (missing 6 non progressor)
missing.nonprog <- nonprog[which(!nonprog %in% existing.nonproggene.file.id$simple_barcode)]
# [1] "TCGA-CN-A63V" "TCGA-CN-A63Y" "TCGA-CN-A640" "TCGA-CQ-7064" "TCGA-CQ-A4CG" "TCGA-CR-6480"
# clinical data available but not omics data 

# nonprogressor patients case ------------------------------------------------------
gene.nonprog.filenames <- lapply(gene.mappings.nonprog[,2], function(x) grep(x, gene.rnaseqv2, value = TRUE))
gene.nonprog.filenames <- lapply(gene.nonprog.filenames, function(x) grep("genes.normalized", x, value = TRUE))
# one file is stored at the wrong directory by previous researcher -----------------
gene.nonprog.filenames <- gene.nonprog.filenames[-which(lapply(gene.nonprog.filenames, function(x) length(x) < 1) == T)]
gene.nonprog.files <- lapply(gene.nonprog.filenames, function(x) read.delim(paste(rnaseqv2.gene.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(gene.nonprog.files) <- lapply(gene.nonprog.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(gene.nonprog.files, file = "data/gene_nonprog_files.RData")

# test ------------------------------
# names(gene.nonprog.files) %in% nonprog
# annot.clinical$Progression_FINAL[annot.clinical$bcr.patient.barcode == "TCGA-CQ-7065"]

# progressor patients case ------------------------------------------------------
gene.prog.filenames <- lapply(gene.mappings.prog[,2], function(x) grep(x, gene.rnaseqv2, value = TRUE))
gene.prog.filenames <- lapply(gene.prog.filenames, function(x) grep("genes.normalized", x, value = TRUE))
gene.prog.files <- lapply(gene.prog.filenames, function(x) read.delim(paste(rnaseqv2.gene.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(gene.prog.files) <- lapply(gene.prog.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(gene.prog.files, file = "data/gene_prog_files.RData")

# test -----------------------------
# names(gene.prog.files) %in% prog
# read exon files -----------------------------------------------------------------
# all patients case ---------------------------------------------------------------
exon.rnaseqv2 <- list.files(path = rnaseqv2.exon.path)
exon.filenames <- lapply(exon.mappings[,2], function(x) grep(x, exon.rnaseqv2, value = TRUE))
exon.filenames <- lapply(exon.filenames, function(x) grep("exon", x, value = TRUE))
exon.files <- lapply(exon.filenames, function(x) read.delim(paste(rnaseqv2.exon.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(exon.files) <- lapply(exon.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(exon.files, file = "data/exon_files.RData")

# nonprogressor patients case ---------------------------------------------------------------
exon.rnaseqv2 <- list.files(path = rnaseqv2.exon.path)
exon.nonprog.filenames <- lapply(exon.mappings.nonprog[,2], function(x) grep(x, exon.rnaseqv2, value = TRUE))
exon.nonprog.filenames <- lapply(exon.nonprog.filenames, function(x) grep("exon", x, value = TRUE))
exon.nonprog.files <- lapply(exon.nonprog.filenames, function(x) read.delim(paste(rnaseqv2.exon.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(exon.nonprog.files) <- lapply(exon.nonprog.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(exon.nonprog.files, file = "data/exon_nonprog_files.RData")

# nonprogressor patients case ---------------------------------------------------------------
exon.rnaseqv2 <- list.files(path = rnaseqv2.exon.path)
exon.prog.filenames <- lapply(exon.mappings.prog[,2], function(x) grep(x, exon.rnaseqv2, value = TRUE))
exon.prog.filenames <- lapply(exon.prog.filenames, function(x) grep("exon", x, value = TRUE))
exon.prog.files <- lapply(exon.prog.filenames, function(x) read.delim(paste(rnaseqv2.exon.path, x[1], sep = "/"), header = T, sep = "\t", stringsAsFactors = F))
names(exon.prog.files) <- lapply(exon.prog.filenames, function(x) rnaseqv2.filename.map[which(rnaseqv2.filename.map$filename %in% x), 3])
save(exon.prog.files, file = "data/exon_prog_files.RData")

#---------------------------------------------------------------------------------------------
rna.not.vailable <- c("TCGA-CN-A63V", "TCGA-CN-A63Y", "TCGA-CN-A640", "TCGA-CQ-7064", "TCGA-CQ-A4CG", "TCGA-CR-6480")
annot.clinical <- annot.clinical[which(!annot.clinical$bcr.patient.barcode %in% rna.not.vailable),]
# filter data uid by features to count ----------------------------------------
# progressed ------------------------------------------------------------------
nonprogressor <- annot.clinical$bcr.patient.barcode[which(annot.clinical$Progression_FINAL == "NonProgressor")]
progressor <- annot.clinical$bcr.patient.barcode[which(annot.clinical$Progression_FINAL == "Progressor")]

# gender -----------------------------------------------------------------------
female <- filter(annot.clinical, gender == "FEMALE")[[1]]
male <- filter(annot.clinical, gender == "MALE")[[1]]

# race --------------------------------------------------------------------------
white <- filter(annot.clinical, race == "WHITE")[[1]]
black <- filter(annot.clinical, race == "BLACK OR AFRICAN AMERICAN")[[1]]
asian <- filter(annot.clinical, race == "ASIAN")[[1]]
native <- filter(annot.clinical, race == "AMERICAN INDIAN OR ALASKA NATIVE")[[1]]

# vital status -------------------------------------------------------------------
alive <- filter(annot.clinical, vital_status_final == "Alive")[[1]]
dead <- filter(annot.clinical, vital_status_final == "Dead")[[1]]

# pathology TNM stage -------------------------------------------------------------
t1.t2 <- clinical.dat$bcr_patient_barcode[grep("T1", annot.clinical$clinical_T)] 
t3 <- clinical.dat$bcr_patient_barcode[grep("T3", annot.clinical$clinical_T)] 
t4 <- clinical.dat$bcr_patient_barcode[grep("T4", annot.clinical$clinical_T)] 
tx <- clinical.dat$bcr_patient_barcode[grep("TX", annot.clinical$clinical_T)] 

n0 <- clinical.dat$bcr_patient_barcode[grep("N0", annot.clinical$clinical_N)] 
nx <- clinical.dat$bcr_patient_barcode[!clinical.dat$bcr_patient_barcode %in% n0] 

# alcohol history ------------------------------------------------------------------
alcohol.yes <- filter(annot.clinical, alcohol_history_documented  == "YES")[[1]]
alcohol.no <- filter(annot.clinical, alcohol_history_documented  == "NO")[[1]]
alcohol.na <- filter(annot.clinical, alcohol_history_documented  == "NA")[[1]]

# tobacco history ------------------------------------------------------------------
tobacco.yes <- filter(annot.clinical, tobacco_smoking_history_indicator  == "Yes")[[1]]
tobacco.no <- filter(annot.clinical, tobacco_smoking_history_indicator  == "No")[[1]]
tobacco.na <- filter(annot.clinical, tobacco_smoking_history_indicator  == "NA")[[1]]

# hpv status -----------------------------------------------------------------------
hpv.yes <- filter(annot.clinical, hpv_status_p16_or_ish  == "Positive")[[1]]
hpv.no <- filter(annot.clinical, hpv_status_p16_or_ish == "Negative")[[1]]
hpv.na <- filter(annot.clinical, hpv_status_p16_or_ish == "NA")[[1]]

# anatomical region ----------------------------------------------------------------
oralcavity <- filter(annot.clinical, anatomic_organ_subdivision  == "Oral Cavity")[[1]]
oropharynx <- filter(annot.clinical, anatomic_organ_subdivision  == "Oropharynx")[[1]]
larynx <- filter(annot.clinical, anatomic_organ_subdivision  == "Larynx")[[1]]

# pathology years -------------------------------------------------------------------
pathafter2010.yes  <- filter(clinical.dat, initial_pathologic_dx_year >= 2010 )[[1]]
pathafter2010.no <- filter(clinical.dat, initial_pathologic_dx_year < 2010 )[[1]]

path.year.plot <- ggplot(clinical.dat, aes(x = initial_pathologic_dx_year , y = ..count..)) + 
  geom_bar(aes(fill = ..count..)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Initial pathologic dx year") + 
  ylab("Frequency") + 
  ggtitle("Distribution of Initial pathologic dx year")


# count functions ---------------------------------------------------------------------
percent.count <- function(set) {
  paste(round( length(set) / 229 * 100, 2), "%", sep = " ")
}

raw.count <- function(set) {
  length(set)
}

# feature counts -----------------------------------------------------------------------
features <- list(annot.clinical[[1]], progressor, nonprogressor, alcohol.yes, alcohol.no, alcohol.na, 
                 tobacco.yes, tobacco.no, tobacco.na, hpv.yes, hpv.no, hpv.na, oralcavity, oropharynx,
                 larynx, female, male)

names(features) <- c("all", "progressor", "nonprogressor", "alcohol.yes", "alcohol.no", "alcohol.na", 
                     "tobacco.yes", "tobacco.no", "tobacco.na", "hpv.yes", "hpv.no", "hpv.na", "oralcavity", "oropharynx",
                     "larynx", "female", "male")


cld <- cbind(names(features), melt(unlist(lapply(features, raw.count))), 
             melt(unlist(lapply(features, percent.count))))
names(cld) <- c("feature", "rawcount", "percent")

cld$feature <- factor(cld$feature, levels = cld$feature[order(cld$rawcount)])

cld.plt <- ggplot(cld, aes(x = rawcount, y = feature, size = rawcount)) + 
  geom_point(aes(colour = rawcount), stat = "identity") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + 
  xlab("Raw count") + 
  ylab("Clinical Feature") + 
  ggtitle("Clinical Features Data-Point Count")

#browseVignettes("gridExtra")
tt <- ttheme_default(core = list(fg_params = list(hjust = 0, x = 0.1, cex = 0.8)),
                     colhead = list(fg_params = list(hjust = 0, x = 0.1, cex = 0.8)),
                     rowhead = list(fg_params = list(hjust = 0, x = 0.1, cex = 0.8)))
clinical.tbl.1 <- tableGrob(slice(arrange(cld, desc(cld$feature)), 1:15), theme = tt)
clinical.tbl.2 <- tableGrob(slice(arrange(cld, desc(cld$feature)), 16:31), theme = tt)

pdf(file = "hnsc_annotated_counts.pdf")
grid.arrange(cld.plt)
grid.arrange(clinical.tbl.1, clinical.tbl.2, ncol = 2)
dev.off()
#----------------------------------------------------------------------------------------
