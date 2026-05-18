# LYNCH SYNDROME PROJECT, COLLABORATION AMSTERDAM GRONINGEN
# by Femke Prins, Q1 2024
# script for loading and preparing Lynch data

setwd("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data")

#### Loading libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(vegan)
library(ggplot2)
library(ggpubr)
library(tibble)
library(emmeans)
library(gtsummary)
library(table1)
library(data.table)
library(flextable)
library(ggrepel)
library(readr)

#### Loading data ----
Lynch_metaphlan <- readRDS("Lynch_metaphlan_complete.rds")
Lynch_reads <- read_csv("Lynch_reads_complete.csv")
Lynch_pathways <- readRDS("Lynch_pathways_complete.rds")
meta_Lynch <- read_excel("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Metadata_Lynch.xlsx")

#### Filter functions (from Ranko) ----

# function for filtering microbiome (metaphlan) results
# - takes dataframe (any)
# - replaces NAs with 0.0
# - filters OUT taxa present in less then presPerc samples
# - filters OUT taxa with mean releative abundance < minMRelAb
# - filters OUT taxa with median relative abundance < minMedRelAb
# - filters OUT taxonomic levels not in keepLevels ()
#   -> keepLevels should be vector of following: T = strain, S = species, G = Genera, F = families
#                                                O = Orders, C = classes, P = phyla, K = Kingdoms
#   -> example: keepLevels=c('S','G') keeps only species and genera
# - filters OUT domains not in keepDomains: 
#   -> input is vector of following: Eukaryota, Bacteria, Viruses, Archaea
#   -> example: keepDomains=c('B') keeps only bacteria, removes viruses, archaea and eukarya
# - if rescaleTaxa = True, rescales relative abundancies after filtering
# returns modified dataframe
# NOTES:
# - assumes metaphlan encoding (k__<kingdom>.o__<order> ... ); it will mess stuff
# up if non-metagenome rows are encoded like this!
# DEFAULTS: 
# - removes non-bacteria
# - keeps all except kingdoms

filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0,
                               rescaleTaxa=T,verbose=T,
                               keepDomains=c('Bacteria','Archaea'),
                               keepLevels=c('S','G','F','O','C','P'),
                               keepUnknown=F) {
  
  # -- drop unknown & rescale? --
  if (!keepUnknown) {
    inDF$UNKNOWN <- NULL
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  tColsNMG = grep('k__',colnames(inDF),invert = T) # colums with microbiome
  
  # replaces NAs in microbiome with 0s
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    nrnZ = as.numeric(sum(inDF[,c]!=0.0))
    if ( (nrnZ/as.numeric(nrow(inDF))) < presPerc) {
      # if (verbose) {
      #   print (paste('col',c,': ',colnames(inDF)[c],'; nr non Zero:',nrnZ,'=',nrnZ/as.numeric(nrow(inDF)),'>> Column removed!'))
      # }
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('^[dgtspcfko]__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = mean(inDF[,c])
    if ( mn < minMRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; mean rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = median(inDF[,c])
    if ( mn < minMedRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; median rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }  
  
  # keep domains
  # -----------------------------
  toKeep <- NULL
  if (length(keepDomains) == 1) {
    if (keepDomains == 'All' | keepDomains == "") {
      toKeep = grep('k__',colnames(inDF),invert = T)
    } else {
      toKeep = grep(paste('k__',d,sep=''),colnames(inDF))
    }
  } else {
    for (d in keepDomains) {
      #d print(d)
      toKeep = c(toKeep,grep(paste('k__',d,sep=''),colnames(inDF)))
    }
  }
  inDF <- inDF[,toKeep]
  
  # remove taxonomic levels
  # -----------------------------
  inDFnonTaxa <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=T)])
  colnames(inDFnonTaxa) <- colnames(inDF)[grep('k__',colnames(inDF),invert=T)]
  
  inDF2 <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=F)])
  # pick strains (T)
  taxaTCols <- grep('t__',colnames(inDF2))
  taxaT <- as.data.frame(inDF2[,taxaTCols])
  colnames(taxaT) <- colnames(inDF2)[grep('t__',colnames(inDF2))]
  if (length(taxaTCols) > 0) {inDF2 <- inDF2[,-taxaTCols]}
  # pick species (S)
  taxaSCols <- grep('s__',colnames(inDF2))
  if (length(taxaSCols) > 0) {
    taxaS <- as.data.frame(inDF2[,taxaSCols])
    colnames(taxaS) <- colnames(inDF2)[grep('s__',colnames(inDF2))]
  }
  if (length(taxaSCols) > 0) {inDF2 <- inDF2[,-taxaSCols]}
  # pick genera (G)
  taxaGCols <- grep('g__',colnames(inDF2))
  if (length(taxaGCols) > 0) {
    taxaG <- as.data.frame(inDF2[,taxaGCols])
    colnames(taxaG) <- colnames(inDF2)[grep('g__',colnames(inDF2))]
  }
  if (length(taxaGCols) > 0) {inDF2 <- inDF2[,-taxaGCols]}
  # pick families (F)
  taxaFCols <- grep('f__',colnames(inDF2))
  if (length(taxaFCols) > 0) {    
    taxaF <- as.data.frame(inDF2[,taxaFCols])
    colnames(taxaF) <- colnames(inDF2)[grep('f__',colnames(inDF2))]
  }
  if (length(taxaFCols) > 0) {inDF2 <- inDF2[,-taxaFCols]}
  # pick orders (O)
  taxaOCols <- grep('o__',colnames(inDF2))
  if (length(taxaOCols) > 0) {
    taxaO <- as.data.frame(inDF2[,taxaOCols])
    colnames(taxaO) <- colnames(inDF2)[grep('o__',colnames(inDF2))]
  }
  if (length(taxaOCols) > 0) {inDF2 <- inDF2[,-taxaOCols]}
  # pick classes (C)
  taxaCCols <- grep('c__',colnames(inDF2))
  if (length(taxaCCols) > 0) {
    taxaC <- as.data.frame(inDF2[,taxaCCols])
    colnames(taxaC) <- colnames(inDF2)[grep('c__',colnames(inDF2))]
  }
  if (length(taxaCCols) > 0) {
    if( length(colnames(inDF2)) - length(taxaCCols) == 1) {
      ccn <- colnames(inDF2)[grep('c__',colnames(inDF2))]
      tempN <- colnames(inDF2)[!(colnames(inDF2) %in% ccn)]
      inDF2 <- as.data.frame(inDF2[,-taxaCCols])
      colnames(inDF2) <- tempN
    } else {
      inDF2 <- inDF2[,-taxaCCols]
    }
  }
  # pick phyla (P)
  taxaPColsKeepN <- NULL
  taxaPCols <- grep('p__',colnames(inDF2))
  if (length(taxaPCols) > 0) {
    taxaPColsN <- colnames(inDF2)[grep('p__',colnames(inDF2))]
    taxaPColsKeep <- grep('p__',colnames(inDF2),invert = T)
    taxaPColsKeepN <- colnames(inDF2)[grep('p__',colnames(inDF2),invert = T)]
    taxaP <- as.data.frame(inDF2[,taxaPCols])
    colnames(taxaP) <- taxaPColsN
  }
  # pick kingdoms (K)
  if (length(taxaPCols) > 0) {inDF2 <- as.data.frame(inDF2[,taxaPColsKeep])}
  if (!is.null(taxaPColsKeepN)){
    colnames(inDF2) <- taxaPColsKeepN
  }
  taxaK <- inDF2
  # pick 
  oDF <- inDFnonTaxa
  if (verbose) {print ('Keeping following taxonomic levels:'); print(keepLevels)}
  if (ncol(taxaK) > 0) {  
    if ('K' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaK),'Kingdoms'))}
      if (rescaleTaxa) {taxaK <- taxaK/rowSums(taxaK)}
      oDF <- cbind(oDF,taxaK)}
  }
  if (ncol(taxaP) > 0) {
    if ('P' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaP),'Phyla'))} 
      if (rescaleTaxa) {taxaP <- taxaP/rowSums(taxaP)}
      oDF <- cbind(oDF,taxaP)}
  }
  if (ncol(taxaC) > 0) {
    if ('C' %in% keepLevels) {
      if (verbose) {print(paste(' -> kept',ncol(taxaC),'Classes'))} 
      if (rescaleTaxa) {taxaC <- taxaC/rowSums(taxaC)}
      oDF <- cbind(oDF,taxaC)}
  }
  if (ncol(taxaO) > 0) {
    if ('O' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaO),'Orders'))}
      if (rescaleTaxa) {taxaO <- taxaO/rowSums(taxaO)}
      oDF <- cbind(oDF,taxaO)}
  }
  if (ncol(taxaF) > 0) {
    if ('F' %in% keepLevels) {
      if (verbose){print(paste(' -> kept',ncol(taxaF),'Families'))}
      if (rescaleTaxa) {taxaF <- taxaF/rowSums(taxaF)}
      oDF <- cbind(oDF,taxaF)}
  }
  if (ncol(taxaG) > 0) { 
    if ('G' %in% keepLevels) {if (verbose){ print(paste(' -> kept',ncol(taxaG),'Genera'))}
      if (rescaleTaxa) {taxaG <- taxaG/rowSums(taxaG)}
      oDF <- cbind(oDF,taxaG)}
  }
  if (ncol(taxaS) > 0) {
    if ('S' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaS),'Species'))}
      if (rescaleTaxa) {taxaS <- taxaS/rowSums(taxaS)}
      oDF <- cbind(oDF,taxaS)}
  }
  if (ncol(taxaT) > 0) {
    if ('T' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaT),'Strains'))}
      if (rescaleTaxa) {taxaT <- taxaT/rowSums(taxaT)}
      oDF <- cbind(oDF,taxaT)}
  }
  if (verbose) {print ('data processing done, returning Dataframe')}
  oDF
}

#Ranko's function, keep first rescaling to make it rel abundance (keeping unmapped/unintegrated in) and then removing unmapped/unintegrated
filterHumannDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,minMedRelAb=0.0,minSum=90.0, rescale=T,verbose=T,type='MetaCyc') {
  
  nonPWYpwys <- c("ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis",
                  "CHLOROPHYLL-SYN: chlorophyllide a biosynthesis I (aerobic, light-dependent)",
                  "GLYCOLYSIS-E-D: superpathway of glycolysis and Entner-Doudoroff",
                  "GLYCOLYSIS-TCA-GLYOX-BYPASS: superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
                  "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
                  "GLYOXYLATE-BYPASS: glyoxylate cycle",
                  "HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
                  "MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS: protein N-glycosylation (eukaryotic, high mannose)",
                  "NAD-BIOSYNTHESIS-II: NAD salvage pathway II",                  
                  "REDCITCYC: TCA cycle VIII (helicobacter)",
                  "TCA-GLYOX-BYPASS: superpathway of glyoxylate bypass and TCA",
                  "TCA: TCA cycle I (prokaryotic)")
  
  colnames(inDF)[colnames(inDF) %in% nonPWYpwys] <- paste0('PWY_',colnames(inDF)[colnames(inDF) %in% nonPWYpwys])
  
  if (type=='MetaCyc') {
    nonPWYdf <- as.data.frame(inDF[,-grep('PWY',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^EC_',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    nonPWYdf <- as.data.frame(inDF[,-grep('RXN',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^PF[01]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^GO',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^K[012]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^K[012]',colnames(inDF))] ])
  }
  colnames(nonPWYdf) <- cnsNonPWYdf
  if (type=='MetaCyc') {
    yesPWYdf <- as.data.frame(inDF[,grep('PWY',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    yesPWYdf <- as.data.frame(inDF[,grep('^EC_',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    yesPWYdf <- as.data.frame(inDF[,grep('RXN',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    yesPWYdf <- as.data.frame(inDF[,grep('^PF[01]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    yesPWYdf <- as.data.frame(inDF[,grep('^GO',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    yesPWYdf <- as.data.frame(inDF[,grep('^K[012]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^K[012]',colnames(inDF))] ])
  }
  
  if (verbose) {print (paste(' > nonPWY',colnames(nonPWYdf))) }
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  # remove Unmapped and unintegrated
  yesPWYdf$PWY_UNMAPPED <- NULL
  yesPWYdf$PWY_UNINTEGRATED <- NULL
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}


#### Making dataframes to work with ----
#Preparing metadata
meta_Lynch_baseline <- subset(meta_Lynch, !is.na(`Pre-scopie sample: Novogene sample name`)) #200
meta_Lynch_baseline <- meta_Lynch_baseline %>% rename(Participant_ID = `Pre-scopie sample: Novogene sample name`)
meta_Lynch_baseline <- subset(meta_Lynch_baseline, !grepl("\\b1\\b", `!! Participation in microbiome study:\r\n1= both pre-scopie and post-scopie sample excluded (patients are marked red)\r\n2= pre-scopie sample included, post-scopie sample excluded \r\n3 = both pre-scopie and post-scopie sample included`))

meta_Lynch_baseline <- merge(meta_Lynch_baseline, Lynch_reads, by = 'Participant_ID') #198 samples
meta_Lynch_baseline <- meta_Lynch_baseline %>% filter(reads >= 1000000) #no samples removed
meta_Lynch_baseline <- as.data.frame(meta_Lynch_baseline)
rownames(meta_Lynch_baseline) <- meta_Lynch_baseline$Participant_ID

#Taxa data: only select the rows that are in the baseline metadata
Lynch_metaphlan_2 <- Lynch_metaphlan %>% filter(rownames(Lynch_metaphlan) %in% meta_Lynch_baseline$Participant_ID)
Lynch_metaphlan_2 <- Lynch_metaphlan_2[,colSums(Lynch_metaphlan_2)>0] #keep only columns that are not only zeros, from 6529 to 5726

#Pathways data: keeping the samples that are in the baseline metadata
rownames(Lynch_pathways) <- sub("_kneaddata.*", "", rownames(Lynch_pathways))
Lynch_pathways_2 <- Lynch_pathways[rownames(Lynch_pathways) %in% meta_Lynch_baseline$Participant_ID, ]
pathways_Lynch_all <- Lynch_pathways_2[,colSums(Lynch_pathways_2)>0] #keep only columns that are not only zeros, from 509 to 482

#### Adding test variables to metadata ----
meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AllNeoplasia_controls = case_when(
    `Relevant lesion? = CRC, adenoma, adv serrated polyp (1=yes, 0=no)` == 1 ~ "Neoplasia",
    `Relevant lesion? = CRC, adenoma, adv serrated polyp (1=yes, 0=no)` == 0 ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AllNeoplasia_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AdvAdenomasCRC_controls = case_when(
    `CRC + Adv adenoma? (1=yes, 0=no, -98 = missing; non-adv adenoma en adv serrated polyp)` == 1 ~ "AdvAdenomasCRC",
    `CRC + Adv adenoma? (1=yes, 0=no, -98 = missing; non-adv adenoma en adv serrated polyp)` == 0 ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AdvAdenomasCRC_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AdvAdenomas_nonAdvAdenomas = case_when(
    `Type of lesion` == 1 ~ "Adv_adenoma",
    `Type of lesion` == 3 ~ "NonAdv_adenoma",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AdvAdenomas_nonAdvAdenomas)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AdvAdenomasCRCAdvSerr_controls = case_when(
    `Type of lesion` %in% c(0, 1, 2) ~ "AdvAdenomasCRCAdvSerr",
    `Type of lesion` %in% c(5, 6, 7) ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AdvAdenomasCRCAdvSerr_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(CRC_controls = case_when(
    `Type of lesion` == 0 ~ "CRC",
    `Type of lesion` %in% c(5, 6, 7) ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$CRC_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(nonAdvAdenomas_controls = case_when(
    `Type of lesion` == 3 ~ "NonAdv_adenoma",
    `Type of lesion` %in% c(5, 6, 7) ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$nonAdvAdenomas_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>% rename(Gene_mutation = `Gene mutation (0=MLH1, 1=MSH2, 2=MSH6, 3=PMS2, 4=non-pathogenic, 5=EPCAM, \r\n-98=unknown`)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(LowRisk_HighRisk = case_when(
    Gene_mutation == 2 | Gene_mutation == 3 ~ "Low-risk",
    Gene_mutation == 0 | Gene_mutation == 1 | Gene_mutation == 5 ~ "High-risk",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$LowRisk_HighRisk)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AllGenes = case_when(
    Gene_mutation == 0 ~ "MLH1",
    Gene_mutation == 1 | Gene_mutation == 5 ~ "MSH2/EPCAM",
    Gene_mutation == 2 ~ "MSH6",
    Gene_mutation == 3 ~ "PMS2",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AllGenes)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(LR_Neoplasia_controls = case_when(
    LowRisk_HighRisk == "High-risk" ~ NA,
    LowRisk_HighRisk == "Low-risk" ~ AllNeoplasia_controls,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$LR_Neoplasia_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(HR_Neoplasia_controls = case_when(
    LowRisk_HighRisk == "Low-risk" ~ NA,
    LowRisk_HighRisk == "High-risk" ~ AllNeoplasia_controls,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$HR_Neoplasia_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NoNeo_High_Lowrisk = case_when(
    AllNeoplasia_controls == "Neoplasia" ~ NA,
    AllNeoplasia_controls == "Control" ~ LowRisk_HighRisk,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NoNeo_High_Lowrisk)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NoNeo_Mutation = case_when(
    AllNeoplasia_controls == "Neoplasia" ~ NA,
    AllNeoplasia_controls == "Control" ~ AllGenes,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NoNeo_Mutation)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NohistoryCRC_LR_Neoplasia_controls = case_when(
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 0 & LowRisk_HighRisk == "High-risk" ~ NA_character_,
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 0 & LowRisk_HighRisk == "Low-risk" ~ AllNeoplasia_controls,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NohistoryCRC_LR_Neoplasia_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NohistoryCRC_HR_Neoplasia_controls = case_when(
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 0 & LowRisk_HighRisk == "Low-risk" ~ NA_character_,
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 0 & LowRisk_HighRisk == "High-risk" ~ AllNeoplasia_controls,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NohistoryCRC_HR_Neoplasia_controls)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(Lynch_historyCRC = case_when(
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 0 ~ "No",
    `Pre-scopie sample:  history CRC (1=yes, 0=no)` == 1 ~ NA,
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$Lynch_historyCRC)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NeoplasiaType_control = case_when(
    `Type of lesion` == 0 ~ "CRC",
    `Type of lesion` == 1 ~ "Adv_adenoma",
    `Type of lesion` == 2 ~ "Adv_serrated",
    `Type of lesion` == 3 ~ "NonAdv_adenoma",
    `Type of lesion` %in% c(5, 6, 7) ~ "Control",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NeoplasiaType_control)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(NohistoryCRC_endoscopy1 = case_when(
    Lynch_historyCRC == "No" & `Feces pre-scopie 1: antal eerdere scopieen (0=0, 1=1, 2=2+)` <= 1 ~ "NoCRC_endo1",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$NohistoryCRC_endoscopy1)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(UMCGLynch = case_when(
    str_detect(`Study number`, "^006") ~ "UMCG",
    TRUE ~ "Non-UMCG"))
table(meta_Lynch_baseline$UMCGLynch)

#### Preparing some variables ----
meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(Bowel_Resection = case_when(
    `Pre-scopie sample:  \r\nhistory bowel resection` == 1 | `Pre-scopie sample:  \r\nhistory bowel resection` == 2 | `Pre-scopie sample:  \r\nhistory bowel resection` == 3 | `Pre-scopie sample:  \r\nhistory bowel resection` == "1 + 3 (counting as 1)" | `Pre-scopie sample:  \r\nhistory bowel resection` == "3 + 4 (couting as 3)"~ "Resection",
    `Pre-scopie sample:  \r\nhistory bowel resection` == 0 ~ "No_resection",
    TRUE ~ NA_character_))

meta_Lynch_baseline <- meta_Lynch_baseline %>% mutate(`gender (1=female)` = ifelse(`gender (1=female)` == 1, "F", "M"))
meta_Lynch_baseline <- meta_Lynch_baseline %>% mutate(`Pre-scopie sample:   smoking` = ifelse(`Pre-scopie sample:   smoking` == 0, "current_smoker", ifelse(`Pre-scopie sample:   smoking` == 1, "ex_smoker", ifelse(`Pre-scopie sample:   smoking` == 2, "never_smoked", NA))))
meta_Lynch_baseline$`Pre-scopie sample: BMI` <- as.numeric(meta_Lynch_baseline$`Pre-scopie sample: BMI`)

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  rename(
    Age = `Pre-scopie sample: age patient`,
    Sex = `gender (1=female)`,
    BMI = `Pre-scopie sample: BMI`,
    Smoking = `Pre-scopie sample:   smoking`,
    Bristol_score = `Pre-scopie sample: BSC (-98 missing)`)

#adjust bristol stool scale 0 = 1
meta_Lynch_baseline$Bristol_score <- as.numeric(meta_Lynch_baseline$Bristol_score)
meta_Lynch_baseline$Bristol_score <- meta_Lynch_baseline$Bristol_score + 1
meta_Lynch_baseline$Bristol_score[meta_Lynch_baseline$Bristol_score == -97] <- 4
table(meta_Lynch_baseline$Bristol_score)

library(table1)
table1(~ Sex + Age + BMI + Smoking + 
         reads + Bristol_score | AllNeoplasia_controls, data=meta_Lynch_baseline)

