# LYNCH SYNDROME PROJECT, COLLABORATION AMSTERDAM GRONINGEN
# by Femke Prins, Q2 2024
# script for filtering and preparing dataframes for analyses

setwd("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data")

#### Input data ----
# Lynch controls
Control_samples <- rownames(meta_Lynch_baseline[meta_Lynch_baseline$AllNeoplasia_controls == "Control", ])
Controls_Lynch_metaphlan <- Lynch_metaphlan_2[rownames(Lynch_metaphlan_2) %in% Control_samples,]
Controls_Lynch_metaphlan <- Controls_Lynch_metaphlan[,colSums(Controls_Lynch_metaphlan)>0]
Controls_Lynch_pathways <- pathways_Lynch_all[rownames(pathways_Lynch_all) %in% Control_samples,]
Controls_Lynch_pathways <- Controls_Lynch_pathways[,colSums(Controls_Lynch_pathways)>0]

# Lynch cases
Neoplasia_samples <- rownames(meta_Lynch_baseline[meta_Lynch_baseline$AllNeoplasia_controls == "Neoplasia", ])
Neoplasia_Lynch_metaphlan <- Lynch_metaphlan_2[rownames(Lynch_metaphlan_2) %in% Neoplasia_samples,]
Neoplasia_Lynch_metaphlan <- Neoplasia_Lynch_metaphlan[,colSums(Neoplasia_Lynch_metaphlan)>0]
Neoplasia_Lynch_pathways <- pathways_Lynch_all[rownames(pathways_Lynch_all) %in% Neoplasia_samples,]
Neoplasia_Lynch_pathways <- Neoplasia_Lynch_pathways[,colSums(Neoplasia_Lynch_pathways)>0]

# General population
GP_metaphlan <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_metaphlan_complete.rds")
GP_pathways <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_pathways_complete.rds")

# LifeLines controls
Controls_LL_metaphlan <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_metaphlan_complete.rds")
Controls_LL_pathways <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_pathways_complete.rds")

# LifeLines cases
Cases_LL_metaphlan <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_metaphlan_complete.rds")
Cases_LL_pathways <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_pathways_complete.rds")

#### Functions for filtering (from Ranko) ----
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

#### Filtering species per group ----
species_Controls_Lynch_filt <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                         keepDomains=c("Bacteria", "NA"),keepLevels=c('S'), 
                                         keepUnknown=F) #290

species_Neoplasia_Lynch_filt <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                       keepDomains=c("Bacteria", "NA"),keepLevels=c('S'), 
                                       keepUnknown=F) #260

species_GP_filt <- filterMetaGenomeDF(GP_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"),keepLevels=c('S'), 
                                                   keepUnknown=F) #259

species_Controls_LL_filt <- filterMetaGenomeDF(Controls_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"),keepLevels=c('S'), 
                                                   keepUnknown=F) #261

species_Cases_LL_filt <- filterMetaGenomeDF(Cases_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"),keepLevels=c('S'), 
                                                   keepUnknown=F) #258

#### Filtering strains per group ----
strains_Controls_Lynch_filt <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                  keepDomains=c("Bacteria", "NA"),keepLevels=c('T'), 
                                                  keepUnknown=F) #307

strains_Neoplasia_Lynch_filt <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"),keepLevels=c('T'), 
                                                   keepUnknown=F) #274

strains_GP_filt <- filterMetaGenomeDF(GP_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                      keepDomains=c("Bacteria", "NA"),keepLevels=c('T'), 
                                      keepUnknown=F) #267

strains_Controls_LL_filt <- filterMetaGenomeDF(Controls_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                               keepDomains=c("Bacteria", "NA"),keepLevels=c('T'), 
                                               keepUnknown=F) #268

strains_Cases_LL_filt <- filterMetaGenomeDF(Cases_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                            keepDomains=c("Bacteria", "NA"),keepLevels=c('T'), 
                                            keepUnknown=F) #265

#### Filtering genera per group ----
genera_Controls_Lynch_filt <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                  keepDomains=c("Bacteria", "NA"),keepLevels=c('G'), 
                                                  keepUnknown=F) #195

genera_Neoplasia_Lynch_filt <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"),keepLevels=c('G'), 
                                                   keepUnknown=F) #179

genera_GP_filt <- filterMetaGenomeDF(GP_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                      keepDomains=c("Bacteria", "NA"),keepLevels=c('G'), 
                                      keepUnknown=F) #175

genera_Controls_LL_filt <- filterMetaGenomeDF(Controls_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                               keepDomains=c("Bacteria", "NA"),keepLevels=c('G'), 
                                               keepUnknown=F) #176

genera_Cases_LL_filt <- filterMetaGenomeDF(Cases_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                            keepDomains=c("Bacteria", "NA"),keepLevels=c('G'), 
                                            keepUnknown=F) #177

#### Filtering phyla per group ----
phyla_Controls_Lynch_filt <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                 keepDomains=c("Bacteria", "NA"),keepLevels=c('P'), 
                                                 keepUnknown=F) #10

phyla_Neoplasia_Lynch_filt <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                  keepDomains=c("Bacteria", "NA"),keepLevels=c('P'), 
                                                  keepUnknown=F) #10

phyla_GP_filt <- filterMetaGenomeDF(GP_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                     keepDomains=c("Bacteria", "NA"),keepLevels=c('P'), 
                                     keepUnknown=F) #9

phyla_Controls_LL_filt <- filterMetaGenomeDF(Controls_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                              keepDomains=c("Bacteria", "NA"),keepLevels=c('P'), 
                                              keepUnknown=F) #9

phyla_Cases_LL_filt <- filterMetaGenomeDF(Cases_LL_metaphlan,presPerc = 0.20,minMRelAb = 0.01,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                           keepDomains=c("Bacteria", "NA"),keepLevels=c('P'), 
                                           keepUnknown=F) #9

#### Filtering pathways per group ----
pathways_Controls_Lynch_filt <- filterHumannDF(Controls_Lynch_pathways,presPerc = 0.20,minMRelAb = 0.0001,minMedRelAb=0.0, rescale=T, verbose=T) #146 pathways

pathways_Neoplasia_Lynch_filt <- filterHumannDF(Neoplasia_Lynch_pathways,presPerc = 0.20,minMRelAb = 0.0001,minMedRelAb=0.0, rescale=T, verbose=T) #147 pathways

pathways_GP_filt <- filterHumannDF(GP_pathways,presPerc = 0.20,minMRelAb = 0.0001,minMedRelAb=0.0, rescale=T, verbose=T) #128 pathways

pathways_Controls_LL_filt <- filterHumannDF(Controls_LL_pathways,presPerc = 0.20,minMRelAb = 0.0001,minMedRelAb=0.0, rescale=T, verbose=T) #128 pathways

pathways_Cases_LL_filt <- filterHumannDF(Cases_LL_pathways,presPerc = 0.20,minMRelAb = 0.0001,minMedRelAb=0.0, rescale=T, verbose=T) #133 pathways

#### Making lists of features to include in each comparison ----
#Species
Controls_Lynch_s <- c(colnames(species_Controls_Lynch_filt)) #290

Neoplasia_Lynch_s <- c(colnames(species_Neoplasia_Lynch_filt)) #260

All_Lynch_s <- c(colnames(species_Controls_Lynch_filt), colnames(species_Neoplasia_Lynch_filt))
All_Lynch_s <- unique(All_Lynch_s) #303

Lynch_GP_s <- c(colnames(species_Controls_Lynch_filt), colnames(species_Neoplasia_Lynch_filt), colnames(species_GP_filt))
Lynch_GP_s <- unique(Lynch_GP_s) #337

Lynch_LLcontrols_s <- c(colnames(species_Controls_Lynch_filt), colnames(species_Controls_LL_filt))
Lynch_LLcontrols_s <- unique(Lynch_LLcontrols_s) #329

Lynch_LLcases_s <- c(colnames(species_Neoplasia_Lynch_filt), colnames(species_Cases_LL_filt))
Lynch_LLcases_s <- unique(Lynch_LLcases_s) #326

#checking unique species in one group
setdiff(colnames(species_GP_filt), colnames(Lynch_metaphlan_2)) -> b #0 only in LifeLines
setdiff(All_Lynch_s, colnames(GP_metaphlan)) -> c #1 only in Lynch

#Strains
Controls_Lynch_t <- c(colnames(strains_Controls_Lynch_filt)) #307

Neoplasia_Lynch_t <- c(colnames(strains_Neoplasia_Lynch_filt)) #274

All_Lynch_t <- c(colnames(strains_Controls_Lynch_filt), colnames(strains_Neoplasia_Lynch_filt))
All_Lynch_t <- unique(All_Lynch_t) #320

Lynch_GP_t <- c(colnames(strains_Controls_Lynch_filt), colnames(strains_Neoplasia_Lynch_filt), colnames(strains_GP_filt))
Lynch_GP_t <- unique(Lynch_GP_t) #352

Lynch_LLcontrols_t <- c(colnames(strains_Controls_Lynch_filt), colnames(strains_Controls_LL_filt))
Lynch_LLcontrols_t <- unique(Lynch_LLcontrols_t) #345

Lynch_LLcases_t <- c(colnames(strains_Neoplasia_Lynch_filt), colnames(strains_Cases_LL_filt))
Lynch_LLcases_t <- unique(Lynch_LLcases_t) #342

#checking unique strains in one group
setdiff(colnames(strains_GP_filt), colnames(Lynch_metaphlan_2)) -> b #0 only in LifeLines
setdiff(All_Lynch_t, colnames(GP_metaphlan)) -> c #1 only in Lynch

#Genera
Controls_Lynch_g <- c(colnames(genera_Controls_Lynch_filt)) #195

Neoplasia_Lynch_g <- c(colnames(genera_Neoplasia_Lynch_filt)) #179

All_Lynch_g <- c(colnames(genera_Controls_Lynch_filt), colnames(genera_Neoplasia_Lynch_filt))
All_Lynch_g <- unique(All_Lynch_g) #204

Lynch_GP_g <- c(colnames(genera_Controls_Lynch_filt), colnames(genera_Neoplasia_Lynch_filt), colnames(genera_GP_filt))
Lynch_GP_g <- unique(Lynch_GP_g) #230

Lynch_LLcontrols_g <- c(colnames(genera_Controls_Lynch_filt), colnames(genera_Controls_LL_filt))
Lynch_LLcontrols_g <- unique(Lynch_LLcontrols_g) #225

Lynch_LLcases_g <- c(colnames(genera_Neoplasia_Lynch_filt), colnames(genera_Cases_LL_filt))
Lynch_LLcases_g <- unique(Lynch_LLcases_g) #222

#checking unique genera in one group
setdiff(colnames(genera_GP_filt), colnames(Lynch_metaphlan_2)) -> b #0 only in LifeLines
setdiff(All_Lynch_g, colnames(GP_metaphlan)) -> c #0 only in Lynch

#Pathways
Controls_Lynch_p <- c(colnames(pathways_Controls_Lynch_filt)) #146

Neoplasia_Lynch_p <- c(colnames(pathways_Neoplasia_Lynch_filt)) #147

All_Lynch_p <- c(colnames(pathways_Controls_Lynch_filt), colnames(pathways_Neoplasia_Lynch_filt))
All_Lynch_p <- unique(All_Lynch_p) #149

Lynch_GP_p <- c(colnames(pathways_Controls_Lynch_filt), colnames(pathways_Neoplasia_Lynch_filt), colnames(pathways_GP_filt))
Lynch_GP_p <- unique(Lynch_GP_p) #153

Lynch_LLcontrols_p <- c(colnames(pathways_Controls_Lynch_filt), colnames(pathways_Controls_LL_filt))
Lynch_LLcontrols_p <- unique(Lynch_LLcontrols_p) #149

Lynch_LLcases_p <- c(colnames(pathways_Neoplasia_Lynch_filt), colnames(pathways_Cases_LL_filt))
Lynch_LLcases_p <- unique(Lynch_LLcases_p) #152

#checking unique pathways in one cohort
setdiff(colnames(pathways_GP_filt), colnames(pathways_Lynch_all)) -> b #0 only in LifeLines
setdiff(All_Lynch_p, colnames(GP_pathways)) -> c #0 only in Lynch

#Phyla
Controls_Lynch_phyla <- c(colnames(phyla_Controls_Lynch_filt)) #10

Neoplasia_Lynch_phyla <- c(colnames(phyla_Neoplasia_Lynch_filt)) #10

All_Lynch_phyla <- c(colnames(phyla_Controls_Lynch_filt), colnames(phyla_Neoplasia_Lynch_filt))
All_Lynch_phyla <- unique(All_Lynch_phyla) #10

Lynch_GP_phyla <- c(colnames(phyla_Controls_Lynch_filt), colnames(phyla_Neoplasia_Lynch_filt), colnames(phyla_GP_filt))
Lynch_GP_phyla <- unique(Lynch_GP_phyla) #10

Lynch_LLcontrols_phyla <- c(colnames(phyla_Controls_Lynch_filt), colnames(phyla_Controls_LL_filt))
Lynch_LLcontrols_phyla <- unique(Lynch_LLcontrols_phyla) #10

Lynch_LLcases_phyla <- c(colnames(phyla_Neoplasia_Lynch_filt), colnames(phyla_Cases_LL_filt))
Lynch_LLcases_phyla <- unique(Lynch_LLcases_phyla) #10

#checking unique phyla in one group
setdiff(colnames(phyla_GP_filt), colnames(Lynch_metaphlan_2)) -> b #0 only in LifeLines
setdiff(All_Lynch_phyla, colnames(GP_metaphlan)) -> c #0 only in Lynch

#### a. Dataframes Lynch (neoplasia + cases) ----
#metadata
meta_Lynch_baseline
write.csv(meta_Lynch_baseline, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_metadata.csv")

#species
All_Lynch_species <- filterMetaGenomeDF(Lynch_metaphlan_2,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                    keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                    keepUnknown=F) #1773

All_Lynch_species_filt <- All_Lynch_species[, colnames(All_Lynch_species) %in% All_Lynch_s] #203 species
All_Lynch_species_filt <- All_Lynch_species_filt / ifelse(rowSums(All_Lynch_species_filt) == 0, 1, rowSums(All_Lynch_species_filt))
write.csv(All_Lynch_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_species.csv")
write.csv(All_Lynch_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_speciesfiltered.csv")

#strains
All_Lynch_strains <- filterMetaGenomeDF(Lynch_metaphlan_2,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                        keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                        keepUnknown=F) #1884

All_Lynch_strains_filt <- All_Lynch_strains[, colnames(All_Lynch_strains) %in% All_Lynch_t] #320 strains
All_Lynch_strains_filt <- All_Lynch_strains_filt / ifelse(rowSums(All_Lynch_strains_filt) == 0, 1, rowSums(All_Lynch_strains_filt))
write.csv(All_Lynch_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_strains.csv")
write.csv(All_Lynch_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_strainsfiltered.csv")

#genera
All_Lynch_genera <- filterMetaGenomeDF(Lynch_metaphlan_2,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                        keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                        keepUnknown=F) #1077

All_Lynch_genera_filt <- All_Lynch_genera[, colnames(All_Lynch_genera) %in% All_Lynch_g] #204
All_Lynch_genera_filt <- All_Lynch_genera_filt / ifelse(rowSums(All_Lynch_genera_filt) == 0, 1, rowSums(All_Lynch_genera_filt))
write.csv(All_Lynch_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_genera.csv")
write.csv(All_Lynch_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_generafiltered.csv")

#phyla
All_Lynch_phyla_all <- filterMetaGenomeDF(Lynch_metaphlan_2,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                       keepDomains=c("Bacteria", "NA"), keepLevels=c('P'), 
                                       keepUnknown=F) #14

All_Lynch_phyla_filt <- All_Lynch_phyla_all[, colnames(All_Lynch_phyla_all) %in% All_Lynch_phyla] #204
All_Lynch_phyla_filt <- All_Lynch_phyla_filt / ifelse(rowSums(All_Lynch_phyla_filt) == 0, 1, rowSums(All_Lynch_phyla_filt))
write.csv(All_Lynch_phyla_all, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_phyla.csv")
write.csv(All_Lynch_phyla_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_phylafiltered.csv")

#pathways
All_Lynch_pathways <- filterHumannDF(pathways_Lynch_all,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #480 pathways, rescaled

All_Lynch_pathways_filt <- All_Lynch_pathways[, colnames(All_Lynch_pathways) %in% All_Lynch_p] #149
All_Lynch_pathways_filt <- All_Lynch_pathways_filt / ifelse(rowSums(All_Lynch_pathways_filt) == 0, 1, rowSums(All_Lynch_pathways_filt))
write.csv(All_Lynch_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_pathways.csv")
write.csv(All_Lynch_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_pathwaysfiltered.csv")

#### b. Dataframes Lynch (neoplasia) ----
#metadata
meta_Lynch_baseline_neoplasia <- meta_Lynch_baseline %>% filter(AllNeoplasia_controls == "Neoplasia")
write.csv(meta_Lynch_baseline_neoplasia, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_metadata.csv")

#species
Neoplasia_Lynch_species <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                        keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                        keepUnknown=F) #1299

Neoplasia_Lynch_species_filt <- Neoplasia_Lynch_species[, colnames(Neoplasia_Lynch_species) %in% Neoplasia_Lynch_s] #260 species
Neoplasia_Lynch_species_filt <- Neoplasia_Lynch_species_filt / ifelse(rowSums(Neoplasia_Lynch_species_filt) == 0, 1, rowSums(Neoplasia_Lynch_species_filt))
write.csv(Neoplasia_Lynch_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_species.csv")
write.csv(Neoplasia_Lynch_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_speciesfiltered.csv")

#strains
Neoplasia_Lynch_strains <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                        keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                        keepUnknown=F) #1355

Neoplasia_Lynch_strains_filt <- Neoplasia_Lynch_strains[, colnames(Neoplasia_Lynch_strains) %in% Neoplasia_Lynch_t] #274 strains
Neoplasia_Lynch_strains_filt <- Neoplasia_Lynch_strains_filt / ifelse(rowSums(Neoplasia_Lynch_strains_filt) == 0, 1, rowSums(Neoplasia_Lynch_strains_filt))
write.csv(Neoplasia_Lynch_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_strains.csv")
write.csv(Neoplasia_Lynch_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_strainsfiltered.csv")

#genera
Neoplasia_Lynch_genera <- filterMetaGenomeDF(Neoplasia_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                       keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                       keepUnknown=F) #789

Neoplasia_Lynch_genera_filt <- Neoplasia_Lynch_genera[, colnames(Neoplasia_Lynch_genera) %in% Neoplasia_Lynch_g] #179
Neoplasia_Lynch_genera_filt <- Neoplasia_Lynch_genera_filt / ifelse(rowSums(Neoplasia_Lynch_genera_filt) == 0, 1, rowSums(Neoplasia_Lynch_genera_filt))
write.csv(Neoplasia_Lynch_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_genera.csv")
write.csv(Neoplasia_Lynch_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_generafiltered.csv")

#pathways
Lynch_Neoplasia_pathways <- filterHumannDF(Neoplasia_Lynch_pathways,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #451 pathways, rescaled

Neoplasia_Lynch_pathways_filt <- Lynch_Neoplasia_pathways[, colnames(Lynch_Neoplasia_pathways) %in% Neoplasia_Lynch_p] #147
Neoplasia_Lynch_pathways_filt <- Neoplasia_Lynch_pathways_filt / ifelse(rowSums(Neoplasia_Lynch_pathways_filt) == 0, 1, rowSums(Neoplasia_Lynch_pathways_filt))
write.csv(Lynch_Neoplasia_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_pathways.csv")
write.csv(Neoplasia_Lynch_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_pathwaysfiltered.csv")

#### c. Dataframes Lynch (controls) ----
#metadata
meta_Lynch_baseline_control <- meta_Lynch_baseline %>% filter(AllNeoplasia_controls == "Control")
write.csv(meta_Lynch_baseline_control, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_metadata.csv")

#species
Control_Lynch_species <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                            keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                            keepUnknown=F) #1677

Control_Lynch_species_filt <- Control_Lynch_species[, colnames(Control_Lynch_species) %in% Controls_Lynch_s] #290 species
Control_Lynch_species_filt <- Control_Lynch_species_filt / ifelse(rowSums(Control_Lynch_species_filt) == 0, 1, rowSums(Control_Lynch_species_filt))
write.csv(Control_Lynch_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_species.csv")
write.csv(Control_Lynch_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_speciesfiltered.csv")

#strains
Control_Lynch_strains <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                            keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                            keepUnknown=F) #1739

Control_Lynch_strains_filt <- Control_Lynch_strains[, colnames(Control_Lynch_strains) %in% Controls_Lynch_t] #307 strains
Control_Lynch_strains_filt <- Control_Lynch_strains_filt / ifelse(rowSums(Control_Lynch_strains_filt) == 0, 1, rowSums(Control_Lynch_strains_filt))
write.csv(Control_Lynch_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_strains.csv")
write.csv(Control_Lynch_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_strainsfiltered.csv")

#genera
Control_Lynch_genera <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                           keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                           keepUnknown=F) #1036

Control_Lynch_genera_filt <- Control_Lynch_genera[, colnames(Control_Lynch_genera) %in% Controls_Lynch_g] #195
Control_Lynch_genera_filt <- Control_Lynch_genera_filt / ifelse(rowSums(Control_Lynch_genera_filt) == 0, 1, rowSums(Control_Lynch_genera_filt))
write.csv(Control_Lynch_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_genera.csv")
write.csv(Control_Lynch_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_generafiltered.csv")

#phyla
Control_Lynch_Phyla <- filterMetaGenomeDF(Controls_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                           keepDomains=c("Bacteria", "NA"), keepLevels=c('P'), 
                                           keepUnknown=F) #14

Control_Lynch_phyla_filt <- Control_Lynch_Phyla[, colnames(Control_Lynch_Phyla) %in% Controls_Lynch_phyla] #10
Control_Lynch_phyla_filt <- Control_Lynch_phyla_filt / ifelse(rowSums(Control_Lynch_phyla_filt) == 0, 1, rowSums(Control_Lynch_phyla_filt))
write.csv(Control_Lynch_Phyla, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_phyla.csv")
write.csv(Control_Lynch_phyla_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_phylafiltered.csv")


#pathways
Lynch_Control_pathways <- filterHumannDF(Controls_Lynch_pathways,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #476 pathways, rescaled

Control_Lynch_pathways_filt <- Lynch_Control_pathways[, colnames(Lynch_Control_pathways) %in% Controls_Lynch_p] #146
Control_Lynch_pathways_filt <- Control_Lynch_pathways_filt / ifelse(rowSums(Control_Lynch_pathways_filt) == 0, 1, rowSums(Control_Lynch_pathways_filt))
write.csv(Lynch_Control_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_pathways.csv")
write.csv(Control_Lynch_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_pathwaysfiltered.csv")

#### d. Dataframes LifeLines (general population) + Lynch (neoplasia + controls) ----
#metadata
meta_GP <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_metadata_complete.rds")
meta_GP <- select(meta_GP, c(DAG3_sampleID, ANTHRO.Sex, ANTHRO.AGE, Smoking, ANTHRO.BMI, MED.MEDS.PPIs_ATC_A02BC, reads, Bristol_score))
meta_GP <- mutate(meta_GP, historyCRC = NA)
meta_GP <- mutate(meta_GP, NohistoryCRC_endoscopy1 = NA)
meta_GP <- mutate(meta_GP, UMCGLynch = NA)
meta_GP <- mutate(meta_GP, cohort = 'General_population')
colnames(meta_GP) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads", "Bristol_score", "historyCRC", "historyCRCendo","UMCGLynch" ,"cohort")
meta_GP$Age <- as.numeric(meta_GP$Age)
meta_GP$reads <- as.numeric(meta_GP$reads)
str(meta_GP)

meta_Lynch_GP <- select(meta_Lynch_baseline, c(Participant_ID, Sex, Age, Smoking, BMI, `Pre-scopie sample: PPI (1=yes, 0=no)`, reads, Bristol_score, Lynch_historyCRC, NohistoryCRC_endoscopy1, UMCGLynch))
meta_Lynch_GP <- mutate(meta_Lynch_GP, cohort = 'Lynch')
colnames(meta_Lynch_GP) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads", "Bristol_score","historyCRC","historyCRCendo","UMCGLynch" ,"cohort")
meta_Lynch_GP$PPI_use <- ifelse(meta_Lynch_GP$PPI_use == 1, "Y", "N")
str(meta_Lynch_GP)

rbind(meta_GP, meta_Lynch_GP) -> meta_GPandLynch
meta_GPandLynch[sapply(meta_GPandLynch, is.character)] <- lapply(meta_GPandLynch[sapply(meta_GPandLynch, is.character)], as.factor)
str(meta_GPandLynch)

meta_GPandLynch <- meta_GPandLynch %>%
  mutate(NohistoryCRC_cohort = case_when(
    historyCRC == "No" & cohort == "Lynch" ~ "Lynch",
    is.na(historyCRC) & cohort == "General_population" ~ "General_population",
    TRUE ~ as.character(NA)))
table(meta_GPandLynch$NohistoryCRC_cohort)

meta_GPandLynch <- meta_GPandLynch %>%
  mutate(NohistoryCRCendo_cohort = case_when(
    historyCRCendo == "NoCRC_endo1" & cohort == "Lynch" ~ "Lynch",
    is.na(historyCRCendo) & cohort == "General_population" ~ "General_population",
    TRUE ~ as.character(NA)))
table(meta_GPandLynch$NohistoryCRCendo_cohort)

meta_GPandLynch <- meta_GPandLynch %>%
  mutate(UMCGcases = case_when(
    UMCGLynch == "UMCG" & cohort == "Lynch" ~ "Lynch",
    is.na(UMCGLynch) & cohort == "General_population" ~ "General_population",
    TRUE ~ as.character(NA)))
table(meta_GPandLynch$UMCGcases)

write.csv(meta_GPandLynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_metadata.csv")

# Merging taxonomic data
GP_metaphlan #general population #7140 variables
Lynch_metaphlan_2 #all of Lynch #5726 variables

intersect(colnames(GP_metaphlan), colnames(Lynch_metaphlan_2)) -> a #5258 overlapping
setdiff(colnames(GP_metaphlan), colnames(Lynch_metaphlan_2)) -> b #1882 only in LifeLines
setdiff(colnames(Lynch_metaphlan_2), colnames(GP_metaphlan)) -> c #468 only in Lynch

taxa_Lynch <- Lynch_metaphlan_2
taxa_Lynch$ID <- rownames(taxa_Lynch)
taxa_GP <- GP_metaphlan
taxa_GP$ID <- rownames(taxa_GP)

GP_Lynch_metaphlan <- full_join(taxa_GP, taxa_Lynch)
rownames(GP_Lynch_metaphlan) <- GP_Lynch_metaphlan$ID
GP_Lynch_metaphlan$ID <- NULL
GP_Lynch_metaphlan[is.na(GP_Lynch_metaphlan)] <- 0 #make all NAs a zero

# Merging pathways
GP_pathways #general population #7140 variables
pathways_Lynch_all #all of Lynch #5726 variables

intersect(colnames(GP_pathways), colnames(pathways_Lynch_all)) -> a #468
setdiff(colnames(GP_pathways), colnames(pathways_Lynch_all)) -> b #106 only in LifeLines
setdiff(colnames(pathways_Lynch_all), colnames(GP_pathways)) -> c #14 only in Lynch

pathway_Lynch <- pathways_Lynch_all
pathway_Lynch$ID <- rownames(pathways_Lynch_all)
pathway_GP <- GP_pathways
pathway_GP$ID <- rownames(GP_pathways)

GP_Lynch_pathways <- full_join(pathway_GP, pathway_Lynch)
rownames(GP_Lynch_pathways) <- GP_Lynch_pathways$ID
GP_Lynch_pathways$ID <- NULL
GP_Lynch_pathways[is.na(GP_Lynch_pathways)] <- 0 #make all NAs a zero

#species
Lynch_GP_species <- filterMetaGenomeDF(GP_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                       keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                       keepUnknown=F) #2373

Lynch_GP_species_filt <- Lynch_GP_species[, colnames(Lynch_GP_species) %in% Lynch_GP_s] #337 species
Lynch_GP_species_filt <- Lynch_GP_species_filt / ifelse(rowSums(Lynch_GP_species_filt) == 0, 1, rowSums(Lynch_GP_species_filt))
write.csv(Lynch_GP_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_species.csv")
write.csv(Lynch_GP_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_speciesfiltered.csv")

#strains
Lynch_GP_strains <- filterMetaGenomeDF(GP_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                       keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                       keepUnknown=F) #2461

Lynch_GP_strains_filt <- Lynch_GP_strains[, colnames(Lynch_GP_strains) %in% Lynch_GP_t] #307 strains
Lynch_GP_strains_filt <- Lynch_GP_strains_filt / ifelse(rowSums(Lynch_GP_strains_filt) == 0, 1, rowSums(Lynch_GP_strains_filt))
write.csv(Lynch_GP_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_strains.csv")
write.csv(Lynch_GP_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_strainsfiltered.csv")

#genera
Lynch_GP_genera <- filterMetaGenomeDF(GP_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                      keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                      keepUnknown=F) #1430

Lynch_GP_genera_filt <- Lynch_GP_genera[, colnames(Lynch_GP_genera) %in% Lynch_GP_g] #195
Lynch_GP_genera_filt <- Lynch_GP_genera_filt / ifelse(rowSums(Lynch_GP_genera_filt) == 0, 1, rowSums(Lynch_GP_genera_filt))
write.csv(Lynch_GP_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_genera.csv")
write.csv(Lynch_GP_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_generafiltered.csv")

#phylum
Lynch_GP_phyla_all <- filterMetaGenomeDF(GP_Lynch_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                      keepDomains=c("Bacteria", "NA"), keepLevels=c('P'), 
                                      keepUnknown=F) #16

Lynch_GP_phyla_filt <- Lynch_GP_phyla_all[, colnames(Lynch_GP_phyla_all) %in% Lynch_GP_phyla] #10
Lynch_GP_phyla_filt <- Lynch_GP_phyla_filt / ifelse(rowSums(Lynch_GP_phyla_filt) == 0, 1, rowSums(Lynch_GP_phyla_filt))
write.csv(Lynch_GP_phyla_all, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_phyla.csv")
write.csv(Lynch_GP_phyla_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_phylafiltered.csv")


#pathways
Lynch_GP_pathways <- filterHumannDF(GP_Lynch_pathways,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #586 pathways, rescaled

Lynch_GP_pathways_filt <- Lynch_GP_pathways[, colnames(Lynch_GP_pathways) %in% Lynch_GP_p] #153
Lynch_GP_pathways_filt <- Lynch_GP_pathways_filt / ifelse(rowSums(Lynch_GP_pathways_filt) == 0, 1, rowSums(Lynch_GP_pathways_filt))
write.csv(Lynch_GP_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_pathways.csv")
write.csv(Lynch_GP_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_pathwaysfiltered.csv")

#### e. Dataframes LifeLines (cases) + Lynch (neoplasia) ----
#metadata
meta_cases <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_metadata_complete.rds")

meta_cases <- meta_cases %>%
  mutate(Smoking = case_when(
    EXP.SMOKING.Smoker.stopped == "Y" ~ "ex_smoker",
    EXP.SMOKING.Smoker.Now == "Y" ~ "current_smoker",
    TRUE ~ "never_smoked"))

meta_cases <- select(meta_cases, c(DAG3_sampleID, ANTHRO.Sex, ANTHRO.AGE, Smoking, ANTHRO.BMI, MED.MEDS.PPIs_ATC_A02BC, reads, new_neoplasia_type, Bristol_score))
meta_cases <- mutate(meta_cases, cohort = 'LifeLines')
colnames(meta_cases) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads", "Neoplasia_type", "Bristol_score","cohort")
meta_cases$reads <- as.numeric(meta_cases$reads)
meta_cases <- meta_cases %>% mutate(Neoplasia_type = recode(Neoplasia_type,
                                 "Adv.Adenoma" = "Adv_adenoma",
                                 "Adv.serrated" = "Adv_serrated",
                                 "Adenoma" = "NonAdv_adenoma",
                                 "CRC" = "CRC"))
str(meta_cases)

meta_neoplasia <- meta_Lynch_baseline %>% filter(AllNeoplasia_controls == "Neoplasia")
meta_neoplasia <- select(meta_neoplasia, c(Participant_ID, Sex, Age, Smoking, BMI, `Pre-scopie sample: PPI (1=yes, 0=no)`, reads, NeoplasiaType_control, Bristol_score))
meta_neoplasia <- mutate(meta_neoplasia, cohort = 'Lynch')
colnames(meta_neoplasia) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads", "Neoplasia_type", "Bristol_score","cohort")
meta_neoplasia$PPI_use <- ifelse(meta_neoplasia$PPI_use == 1, "Y", "N")
str(meta_neoplasia)

rbind(meta_neoplasia, meta_cases) -> meta_casesLynchGP
meta_casesLynchGP[sapply(meta_casesLynchGP, is.character)] <- lapply(meta_casesLynchGP[sapply(meta_casesLynchGP, is.character)], as.factor)
str(meta_casesLynchGP)

#adding variables
meta_casesLynchGP <- meta_casesLynchGP %>%
  mutate(LynchLL_AdvAdenomas = case_when(
    Neoplasia_type == "Adv_adenoma" ~ cohort,
    TRUE ~ NA_character_))
table(meta_casesLynchGP$LynchLL_AdvAdenomas)

meta_casesLynchGP <- meta_casesLynchGP %>%
  mutate(LynchLL_AdvAdenomasCRC = case_when(
    Neoplasia_type == "Adv_adenoma" ~ cohort,
    Neoplasia_type == "CRC" ~ cohort,
    TRUE ~ NA_character_))
table(meta_casesLynchGP$LynchLL_AdvAdenomasCRC)

meta_casesLynchGP <- meta_casesLynchGP %>%
  mutate(LynchLL_AdvNeoplasia = case_when(
    Neoplasia_type == "Adv_adenoma" ~ cohort,
    Neoplasia_type == "CRC" ~ cohort,
    Neoplasia_type == "Adv_serrated" ~ cohort,
    TRUE ~ NA_character_))
table(meta_casesLynchGP$LynchLL_AdvNeoplasia)

meta_casesLynchGP <- meta_casesLynchGP %>%
  mutate(LynchLL_NonAdvAdenomas = case_when(
    Neoplasia_type == "NonAdv_adenoma" ~ cohort,
    TRUE ~ NA_character_))
table(meta_casesLynchGP$LynchLL_NonAdvAdenomas)

meta_casesLynchGP <- meta_casesLynchGP %>%
  mutate(LynchLL_Neoplasia = cohort)
table(meta_casesLynchGP$LynchLL_Neoplasia)

write.csv(meta_casesLynchGP, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_metadata.csv")

# Merging taxonomic data
Cases_LL_metaphlan #general population #3864 variables
Neoplasia_Lynch_metaphlan #all of Lynch #4230 variables

intersect(colnames(Cases_LL_metaphlan), colnames(Neoplasia_Lynch_metaphlan)) -> a #3153 overlapping
setdiff(colnames(Cases_LL_metaphlan), colnames(Neoplasia_Lynch_metaphlan)) -> b #711 only in LifeLines
setdiff(colnames(Neoplasia_Lynch_metaphlan), colnames(Cases_LL_metaphlan)) -> c #1077 only in Lynch

taxa_Lynch_neoplasia <- Neoplasia_Lynch_metaphlan
taxa_Lynch_neoplasia$ID <- rownames(taxa_Lynch_neoplasia)
taxa_LifeLines_cases <- Cases_LL_metaphlan
taxa_LifeLines_cases$ID <- rownames(taxa_LifeLines_cases)

LynchLifeLines_cases_metaphlan <- full_join(taxa_LifeLines_cases, taxa_Lynch_neoplasia)
rownames(LynchLifeLines_cases_metaphlan) <- LynchLifeLines_cases_metaphlan$ID
LynchLifeLines_cases_metaphlan$ID <- NULL
LynchLifeLines_cases_metaphlan[is.na(LynchLifeLines_cases_metaphlan)] <- 0 #make all NAs a zero

# Merging pathways
Cases_LL_pathways #general population #465 variables
Neoplasia_Lynch_pathways #all of Lynch #453 variables

intersect(colnames(Cases_LL_pathways), colnames(Neoplasia_Lynch_pathways)) -> a #433
setdiff(colnames(Cases_LL_pathways), colnames(Neoplasia_Lynch_pathways)) -> b #22 only in LifeLines
setdiff(colnames(Neoplasia_Lynch_pathways), colnames(Cases_LL_pathways)) -> c #10 only in Lynch

pathway_Lynch_neoplasia <- Neoplasia_Lynch_pathways
pathway_Lynch_neoplasia$ID <- rownames(pathway_Lynch_neoplasia)
pathway_LifeLines_cases <- Cases_LL_pathways
pathway_LifeLines_cases$ID <- rownames(pathway_LifeLines_cases)

LynchLifeLines_cases_pathways <- full_join(pathway_Lynch_neoplasia, pathway_LifeLines_cases)
rownames(LynchLifeLines_cases_pathways) <- LynchLifeLines_cases_pathways$ID
LynchLifeLines_cases_pathways$ID <- NULL
LynchLifeLines_cases_pathways[is.na(LynchLifeLines_cases_pathways)] <- 0 #make all NAs a zero

#species
LynchLifeLines_cases_species <- filterMetaGenomeDF(LynchLifeLines_cases_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                                   keepUnknown=F) #1510

LynchLifeLines_cases_species_filt <- LynchLifeLines_cases_species[, colnames(LynchLifeLines_cases_species) %in% Lynch_LLcases_s] #326 species
LynchLifeLines_cases_species_filt <- LynchLifeLines_cases_species_filt / ifelse(rowSums(LynchLifeLines_cases_species_filt) == 0, 1, rowSums(LynchLifeLines_cases_species_filt))
write.csv(LynchLifeLines_cases_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_species.csv")
write.csv(LynchLifeLines_cases_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_speciesfiltered.csv")

#strains
LynchLifeLines_cases_strains <- filterMetaGenomeDF(LynchLifeLines_cases_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                   keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                                   keepUnknown=F) #1570

LynchLifeLines_cases_strains_filt <- LynchLifeLines_cases_strains[, colnames(LynchLifeLines_cases_strains) %in% Lynch_LLcases_t] #342 strains
LynchLifeLines_cases_strains_filt <- LynchLifeLines_cases_strains_filt / ifelse(rowSums(LynchLifeLines_cases_strains_filt) == 0, 1, rowSums(LynchLifeLines_cases_strains_filt))
write.csv(LynchLifeLines_cases_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_strains.csv")
write.csv(LynchLifeLines_cases_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_strainsfiltered.csv")

#genera
LynchLifeLines_cases_genera <- filterMetaGenomeDF(LynchLifeLines_cases_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                  keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                                  keepUnknown=F) #934

LynchLifeLines_cases_genera_filt <- LynchLifeLines_cases_genera[, colnames(LynchLifeLines_cases_genera) %in% Lynch_LLcases_g] #222
LynchLifeLines_cases_genera_filt <- LynchLifeLines_cases_genera_filt / ifelse(rowSums(LynchLifeLines_cases_genera_filt) == 0, 1, rowSums(LynchLifeLines_cases_genera_filt))
write.csv(LynchLifeLines_cases_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_genera.csv")
write.csv(LynchLifeLines_cases_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_generafiltered.csv")

#pathways
Cases_LynchLifelines_pathways <- filterHumannDF(LynchLifeLines_cases_pathways,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #473 pathways, rescaled

LynchLifeLines_cases_pathways_filt <- Cases_LynchLifelines_pathways[, colnames(Cases_LynchLifelines_pathways) %in% Lynch_LLcases_p] #152
LynchLifeLines_cases_pathways_filt <- LynchLifeLines_cases_pathways_filt / ifelse(rowSums(LynchLifeLines_cases_pathways_filt) == 0, 1, rowSums(LynchLifeLines_cases_pathways_filt))
write.csv(Cases_LynchLifelines_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_pathways.csv")
write.csv(LynchLifeLines_cases_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_pathwaysfiltered.csv")

#phylum
LynchLifeLines_cases_phyla <- filterMetaGenomeDF(LynchLifeLines_cases_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                  keepDomains=c("Bacteria", "NA"), keepLevels=c('P'), 
                                                  keepUnknown=F) #15

LynchLifeLines_cases_phyla_filt <- LynchLifeLines_cases_phyla[, colnames(LynchLifeLines_cases_phyla) %in% Lynch_LLcases_phyla] #10
LynchLifeLines_cases_phyla_filt <- LynchLifeLines_cases_phyla_filt / ifelse(rowSums(LynchLifeLines_cases_phyla_filt) == 0, 1, rowSums(LynchLifeLines_cases_phyla_filt))
write.csv(LynchLifeLines_cases_phyla, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_phyla.csv")
write.csv(LynchLifeLines_cases_phyla_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_phylafiltered.csv")


#### f. Dataframes LifeLines (controls) + Lynch (controls) ----
#metadata
meta_controls_LL <- readRDS("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_metadata_complete.rds")

meta_controls_LL <- select(meta_controls_LL, c(DAG3_sampleID, ANTHRO.Sex, ANTHRO.AGE, Smoking, ANTHRO.BMI, MED.MEDS.PPIs_ATC_A02BC, reads, Bristol_score))
meta_controls_LL$historyCRC <- NA
meta_controls_LL$historyCRCendo <- NA
meta_controls_LL <- mutate(meta_controls_LL, cohort = 'LifeLines')
colnames(meta_controls_LL) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads","Bristol_score", "historyCRC", "historyCRCendo","cohort")
meta_controls_LL$reads <- as.numeric(meta_controls_LL$reads)
str(meta_controls_LL)

meta_control_Lynch <- meta_Lynch_baseline %>% filter(AllNeoplasia_controls == "Control")
meta_control_Lynch <- select(meta_control_Lynch, c(Participant_ID, Sex, Age, Smoking, BMI, `Pre-scopie sample: PPI (1=yes, 0=no)`, reads, Bristol_score, Lynch_historyCRC, NohistoryCRC_endoscopy1))
meta_control_Lynch <- mutate(meta_control_Lynch, cohort = 'Lynch')
colnames(meta_control_Lynch) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "PPI_use", "reads", "Bristol_score", "historyCRC", "historyCRCendo" ,"cohort")
meta_control_Lynch$PPI_use <- ifelse(meta_control_Lynch$PPI_use == 1, "Y", "N")
str(meta_control_Lynch)

rbind(meta_controls_LL, meta_control_Lynch) -> meta_controlsLynchGP
meta_controlsLynchGP[sapply(meta_controlsLynchGP, is.character)] <- lapply(meta_controlsLynchGP[sapply(meta_controlsLynchGP, is.character)], as.factor)
str(meta_controlsLynchGP)

meta_controlsLynchGP <- meta_controlsLynchGP %>%
  mutate(LynchLL_Controls = cohort)
table(meta_controlsLynchGP$LynchLL_Controls)

meta_controlsLynchGP <- meta_controlsLynchGP %>%
  mutate(LynchLL_NohistoryCRC_cohort = case_when(
    historyCRC == "No" & cohort == "Lynch" ~ "Lynch",
    is.na(historyCRC) & cohort == "LifeLines" ~ "LifeLines",
    TRUE ~ as.character(NA)))
table(meta_controlsLynchGP$LynchLL_NohistoryCRC_cohort)

meta_controlsLynchGP <- meta_controlsLynchGP %>%
  mutate(LynchLL_NohistoryCRCendo_cohort = case_when(
    historyCRCendo == "NoCRC_endo1" & cohort == "Lynch" ~ "Lynch",
    is.na(historyCRCendo) & cohort == "LifeLines" ~ "LifeLines",
    TRUE ~ as.character(NA)))
table(meta_controlsLynchGP$LynchLL_NohistoryCRCendo_cohort)

#write.csv(meta_controlsLynchGP, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_metadata.csv")

# Merging taxonomic data
Controls_LL_metaphlan #general population #6554 variables
Controls_Lynch_metaphlan #all of Lynch 5456 variables

intersect(colnames(Controls_LL_metaphlan), colnames(Controls_Lynch_metaphlan)) -> a #4942 overlapping
setdiff(colnames(Controls_LL_metaphlan), colnames(Controls_Lynch_metaphlan)) -> b #1612 only in LifeLines
setdiff(colnames(Controls_Lynch_metaphlan), colnames(Controls_LL_metaphlan)) -> c #514 only in Lynch

taxa_Lynch_controls <- Controls_Lynch_metaphlan
taxa_Lynch_controls$ID <- rownames(taxa_Lynch_controls)
taxa_LifeLines_controls <- Controls_LL_metaphlan
taxa_LifeLines_controls$ID <- rownames(taxa_LifeLines_controls)

LynchLifeLines_controls_metaphlan <- full_join(taxa_LifeLines_controls, taxa_Lynch_controls)
rownames(LynchLifeLines_controls_metaphlan) <- LynchLifeLines_controls_metaphlan$ID
LynchLifeLines_controls_metaphlan$ID <- NULL
LynchLifeLines_controls_metaphlan[is.na(LynchLifeLines_controls_metaphlan)] <- 0 #make all NAs a zero

# Merging pathways
Controls_LL_pathways #general population #498 variables
Controls_Lynch_pathways #all of Lynch #478 variables

intersect(colnames(Controls_LL_pathways), colnames(Controls_Lynch_pathways)) -> a #459
setdiff(colnames(Controls_LL_pathways), colnames(Controls_Lynch_pathways)) -> b #39 only in LifeLines
setdiff(colnames(Controls_Lynch_pathways), colnames(Controls_LL_pathways)) -> c #19 only in Lynch

pathway_Lynch_controls <- Controls_Lynch_pathways
pathway_Lynch_controls$ID <- rownames(pathway_Lynch_controls)
pathway_LifeLines_controls <- Controls_LL_pathways
pathway_LifeLines_controls$ID <- rownames(pathway_LifeLines_controls)

LynchLifeLines_controls_pathways <- full_join(pathway_Lynch_controls, pathway_LifeLines_controls)
rownames(LynchLifeLines_controls_pathways) <- LynchLifeLines_controls_pathways$ID
LynchLifeLines_controls_pathways$ID <- NULL
LynchLifeLines_controls_pathways[is.na(LynchLifeLines_controls_pathways)] <- 0 #make all NAs a zero

#species
LynchLifeLines_controls_species <- filterMetaGenomeDF(LynchLifeLines_controls_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                      keepDomains=c("Bacteria", "NA"), keepLevels=c('S'), 
                                                      keepUnknown=F) #2200

LynchLifeLines_controls_species_filt <- LynchLifeLines_controls_species[, colnames(LynchLifeLines_controls_species) %in% Lynch_LLcontrols_s] #329 species
LynchLifeLines_controls_species_filt <- LynchLifeLines_controls_species_filt / ifelse(rowSums(LynchLifeLines_controls_species_filt) == 0, 1, rowSums(LynchLifeLines_controls_species_filt))
#write.csv(LynchLifeLines_controls_species, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_species.csv")
#write.csv(LynchLifeLines_controls_species_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_speciesfiltered.csv")

#strains
LynchLifeLines_controls_strains <- filterMetaGenomeDF(LynchLifeLines_controls_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                      keepDomains=c("Bacteria", "NA"), keepLevels=c('T'), 
                                                      keepUnknown=F) #2279

LynchLifeLines_controls_strains_filt <- LynchLifeLines_controls_strains[, colnames(LynchLifeLines_controls_strains) %in% Lynch_LLcontrols_t] #345 strains
LynchLifeLines_controls_strains_filt <- LynchLifeLines_controls_strains_filt / ifelse(rowSums(LynchLifeLines_controls_strains_filt) == 0, 1, rowSums(LynchLifeLines_controls_strains_filt))
#write.csv(LynchLifeLines_controls_strains, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_strains.csv")
#write.csv(LynchLifeLines_controls_strains_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_strainsfiltered.csv")

#genera
LynchLifeLines_controls_genera <- filterMetaGenomeDF(LynchLifeLines_controls_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                     keepDomains=c("Bacteria", "NA"), keepLevels=c('G'), 
                                                     keepUnknown=F) #1328

LynchLifeLines_controls_genera_filt <- LynchLifeLines_controls_genera[, colnames(LynchLifeLines_controls_genera) %in% Lynch_LLcontrols_g] #225
LynchLifeLines_controls_genera_filt <- LynchLifeLines_controls_genera_filt / ifelse(rowSums(LynchLifeLines_controls_genera_filt) == 0, 1, rowSums(LynchLifeLines_controls_genera_filt))
#write.csv(LynchLifeLines_controls_genera, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_genera.csv")
#write.csv(LynchLifeLines_controls_genera_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_generafiltered.csv")

#pathways
controls_LynchLifelines_pathways <- filterHumannDF(LynchLifeLines_controls_metaphlan,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb=0.0, rescale=T, verbose=T) #515 pathways, rescaled

LynchLifeLines_controls_pathways_filt <- controls_LynchLifelines_pathways[, colnames(controls_LynchLifelines_pathways) %in% Lynch_LLcontrols_p] #149
LynchLifeLines_controls_pathways_filt <- LynchLifeLines_controls_pathways_filt / ifelse(rowSums(LynchLifeLines_controls_pathways_filt) == 0, 1, rowSums(LynchLifeLines_controls_pathways_filt))
#write.csv(controls_LynchLifelines_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_pathways.csv")
#write.csv(LynchLifeLines_controls_pathways_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_pathwaysfiltered.csv")

#phyla
LynchLifeLines_controls_phyla <- filterMetaGenomeDF(LynchLifeLines_controls_metaphlan,presPerc = 0.00,minMRelAb = 0.00,minMedRelAb=0.0,rescaleTaxa=T,verbose=T,
                                                     keepDomains=c("Bacteria", "NA"), keepLevels=c('P'), 
                                                     keepUnknown=F) #16

LynchLifeLines_controls_phyla_filt <- LynchLifeLines_controls_phyla[, colnames(LynchLifeLines_controls_phyla) %in% Lynch_LLcontrols_phyla] #10
LynchLifeLines_controls_phyla_filt <- LynchLifeLines_controls_phyla_filt / ifelse(rowSums(LynchLifeLines_controls_phyla_filt) == 0, 1, rowSums(LynchLifeLines_controls_phyla_filt))
#write.csv(LynchLifeLines_controls_phyla, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_phyla.csv")
#write.csv(LynchLifeLines_controls_phyla_filt, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_phylafiltered.csv")
