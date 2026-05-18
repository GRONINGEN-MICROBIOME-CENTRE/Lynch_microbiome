# LYNCH SYNDROME PROJECT, COLLABORATION AMSTERDAM GRONINGEN
# by Femke Prins, Q2 2024
# script for performing beta-diversity analyses

setwd("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data")

#Load values ----
totest_Lynch <- c("AllNeoplasia_controls",
                  "AdvAdenomasCRC_controls",
                  "AdvAdenomas_nonAdvAdenomas",
                  "CRC_controls",
                  "AdvAdenomasCRCAdvSerr_controls",
                  "nonAdvAdenomas_controls",
                  "LR_Neoplasia_controls",
                  "NohistoryCRC_LR_Neoplasia_controls",
                  "HR_Neoplasia_controls",
                  "NohistoryCRC_HR_Neoplasia_controls",
                  "NoNeo_High_Lowrisk",
                  "NoNeo_Mutation")

totest_Lynch_a <- c("AllNeoplasia_controls",
                    "AdvAdenomasCRC_controls",
                    "CRC_controls",
                    "AdvAdenomasCRCAdvSerr_controls",
                    "nonAdvAdenomas_controls",
                    "LR_Neoplasia_controls",
                    "NohistoryCRC_LR_Neoplasia_controls",
                    "HR_Neoplasia_controls",
                    "NohistoryCRC_HR_Neoplasia_controls")

totest_Lynch_b <- c("AdvAdenomas_nonAdvAdenomas", "AdvAdenomasCRCAdvSerr_nonAdvAdenomas")


totest_Lynch_c <- c("NoNeo_High_Lowrisk",
                    "NoNeo_Mutation")

totest_LynchGP_d <- c("cohort", "NohistoryCRC_cohort", "NohistoryCRCendo_cohort")

totest_LynchGP_e <- c("LynchLL_AdvAdenomas", "LynchLL_AdvAdenomasCRC", 'LynchLL_NonAdvAdenomas', "LynchLL_Neoplasia")

totest_LynchGP_f <- c("LynchLL_Controls", "LynchLL_NohistoryCRC_cohort", "LynchLL_NohistoryCRCendo_cohort")

totest_Lynch_neoplasia <- c("AllNeoplasia_controls","AdvAdenomasCRC_controls","AdvAdenomas_nonAdvAdenomas","CRC_controls", "AdvAdenomasCRCAdvSerr_controls", "nonAdvAdenomas_controls")
totest_Lynch_genes <- c("NoNeo_High_Lowrisk","NoNeo_Mutation")
totest_Lynch_risk <- c("AllNeoplasia_controls","LR_Neoplasia_controls","NohistoryCRC_LR_Neoplasia_controls","HR_Neoplasia_controls","NohistoryCRC_HR_Neoplasia_controls")

covariates_Lynch <- c("Sex", 'reads', "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
covariates_LifeLines <- c("Sex", 'reads', "Age", "BMI", "Smoking", "Bristol_score")

outline_colors <- c("Neoplasia"="#AA1963", "Control"="#41C5B0", 
                    "AdvAdenomasCRC"="#db107b", "NonAdv_adenoma"="#A496CF", 
                    "CRC" = "#B12f73", "Adv_adenoma" = "#c90076",
                    "AdvAdenomasCRCAdvSerr" = "#cb75a1", "Adv_serrated" = "lightpink",
                    "Low-risk"= "#21560a", "High-risk"= "#0b5394", 
                    "MLH1"="#00BFC4", "MSH2/EPCAM"="#719FA6", 
                    "MSH6" = "#7CAE00", "PMS2" = "#B6D7A8",
                    "General_population" = "orange", "Lynch" = "darkblue",
                    "LifeLines" = "darkorange")
my_colours <- outline_colors

adonis_colors <- c("FDR"="darkolivegreen3", "No"="indianred3", "Nominal" = "darkolivegreen1", "Yes" = "darkolivegreen3")

#Load functions ----
perform_beta_analysis_Lynch <- function(variable) {
  # Calculate centroids
  ref <- reformulate(variable, paste("cbind(", paste0("V", 1:5, collapse = ","), ")"))
  centroids <- aggregate(ref, PCoA_meta, mean)
  
  # Create PCoA plot
  PCoA1_2 <- ggplot(PCoA_meta, aes_string(x = "V1", y = "V2", color = variable)) + 
    xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
    ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
    geom_point(size = 3, alpha = 0.9) + 
    theme_light() + 
    stat_ellipse(geom = "polygon", alpha = 0.1) +
    theme(legend.title = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", color = "black"),
          axis.title.x = element_text(size = 12, face = "bold", color = "black"),
          legend.position = "bottom",
          strip.text = element_text(size = 12)) +
    scale_fill_manual(values = my_colours) +
    scale_color_manual(values = outline_colors) +
    geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
    geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = variable), alpha = 0.8, colour = "black") +
    ggtitle(paste("Beta Diversity:", variable)) 
  
  # Return results
  return(list(beta_diversity = PCoA1_2))
}

perform_beta_analysis_Lynch <- function(variable) {
  # Remove rows with NA in the specified variable
  PCoA_meta_filtered <- PCoA_meta[!is.na(PCoA_meta[[variable]]), ]
  
  # Calculate centroids on the filtered data
  ref <- reformulate(variable, paste("cbind(", paste0("V", 1:5, collapse = ","), ")"))
  centroids <- aggregate(ref, PCoA_meta_filtered, mean)
  
  # Create PCoA plot on the filtered data
  PCoA1_2 <- ggplot(PCoA_meta_filtered, aes_string(x = "V1", y = "V2", color = variable)) + 
    xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
    ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
    geom_point(size = 3, alpha = 0.9) + 
    theme_light() + 
    stat_ellipse(geom = "polygon", alpha = 0.1) +
    theme(legend.title = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", color = "black"),
          axis.title.x = element_text(size = 12, face = "bold", color = "black"),
          legend.position = "bottom",
          strip.text = element_text(size = 12)) +
    scale_fill_manual(values = my_colours) +
    scale_color_manual(values = outline_colors) +
    geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
    geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = variable), alpha = 0.8, colour = "black") +
    ggtitle(paste("Beta Diversity:", variable)) 
  
  # Return results
  return(list(beta_diversity = PCoA1_2))
}

perform_multiadonis_analysis_Lynch <- function(variable, meta_data, clr_data, covariates = NULL) {
  # Select relevant covariates
  test_covar <- covariates[!covariates %in% variable]
  
  # Merge phenos and taxa
  allDF <- merge(meta_data %>% select(all_of(c(variable, test_covar, 'Participant_ID'))), clr_data, by = 'row.names')
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  
  # Remove rows with NA in the tested variable
  allDF <- allDF[complete.cases(allDF[, variable]), ]
  
  # Select taxa and phenos
  ad_taxa <- allDF[, grep("^k__", colnames(allDF))]
  ad_test <- allDF[, intersect(colnames(allDF), c(variable, test_covar))]
  
  # Calculate distance matrix
  inBC <- vegdist(ad_taxa, method = "euclidean", parallel = 4) #=Aitchison distance because CLR data
  
  # Perform adonis permanova
  formula <- as.formula(paste("inBC ~ ", paste(c(variable, test_covar), collapse = ' + ')))
  ad <- adonis2(formula, data = ad_test, permutations = permNR, by = "margin")
  adonis_mv <- model.frame(ad)
  
  # Add additional columns
  formula_text <- paste(deparse(formula), collapse = "")
  nrRows <- nrow(ad_test)
  adonis_mv$Variable <- variable
  adonis_mv$Formula <- as.character(formula_text)
  adonis_mv$Significance <- ifelse(adonis_mv$`Pr(>F)` < 0.05, "Yes", "No")
  adonis_mv$nrRows <- nrRows
  
  # Make plot
  adonis_plot <- ggplot(adonis_mv, aes(reorder(row.names(adonis_mv), R2), R2, fill = Significance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    ggtitle(paste("Multivariable: species,", variable)) +
    ylab("Explained variance (R^2)") + xlab("Factor") + 
    scale_fill_manual(values = adonis_colors) +
    theme(text = element_text(size = 10), legend.position = "bottom")
  
  # Store results and plot in lists
  results <- list(results = adonis_mv, plot = adonis_plot)
  return(results)
}
perform_multiadonis_analysis_pathways_Lynch <- function(variable, meta_data, clr_data, covariates = NULL) {
  # Select relevant covariates
  test_covar <- covariates[!covariates %in% variable]
  
  # Merge phenos and taxa
  allDF <- merge(meta_data %>% select(all_of(c(variable, test_covar, 'Participant_ID'))), clr_data, by = 'row.names')
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  
  # Remove rows with NA in the tested variable
  allDF <- allDF[complete.cases(allDF[, variable]), ]
  
  # Select taxa and phenos
  ad_pathway <- allDF[, grep("PWY", colnames(allDF))]
  ad_test <- allDF[, intersect(colnames(allDF), c(variable, test_covar))]
  
  # Calculate distance matrix
  inBC <- vegdist(ad_pathway, method = "euclidean", parallel = 4) #=Aitchison distance because CLR data
  
  # Perform adonis permanova
  formula <- as.formula(paste("inBC ~ ", paste(c(variable, test_covar), collapse = ' + ')))
  ad <- adonis2(formula, data = ad_test, permutations = permNR, by = "margin")
  adonis_mv <- model.frame(ad)
  
  # Add additional columns
  formula_text <- paste(deparse(formula), collapse = "")
  nrRows <- nrow(ad_test)
  adonis_mv$Variable <- variable
  adonis_mv$Formula <- as.character(formula_text)
  adonis_mv$Significance <- ifelse(adonis_mv$`Pr(>F)` < 0.05, "Yes", "No")
  adonis_mv$nrRows <- nrRows
  
  # Make plot
  adonis_plot <- ggplot(adonis_mv, aes(reorder(row.names(adonis_mv), R2), R2, fill = Significance)) + 
    geom_bar(stat = "identity") + coord_flip() + theme_bw() +
    ggtitle(paste("Multivariable: species,", variable)) +
    ylab("Explained variance (R^2)") + xlab("Factor") + 
    scale_fill_manual(values = adonis_colors) +
    theme(text = element_text(size = 10), legend.position = "bottom")
  
  # Store results and plot in lists
  results <- list(results = adonis_mv, plot = adonis_plot)
  return(results)
}

#### a.Beta-diversity species ----
Lynch_all_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_speciesfiltered.csv")
Lynch_all_species_filt <- as.data.frame(Lynch_all_species_filt)
rownames(Lynch_all_species_filt) <- Lynch_all_species_filt[, 1]
Lynch_all_species_filt <- Lynch_all_species_filt[, -1]
Lynch_all_species_filt_clr <- decostand(Lynch_all_species_filt, method = 'clr', pseudocount = min(Lynch_all_species_filt[Lynch_all_species_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_all_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# a. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_a <- list()

for (variable in totest_Lynch_a) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_a[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_a, ncol = 3, nrow = 3)

#Polished plot neoplasia types versus controls
library(ggExtra)

ref <- reformulate("NeoplasiaType_control", paste("cbind(", paste0("V", 1:5, collapse = ","), ")"))
centroids <- aggregate(ref, PCoA_meta, mean)
PCoA1_2_extra <- ggplot(PCoA_meta, aes_string(x = "V1", y = "V2", color = "NeoplasiaType_control")) + 
  xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
  ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = "AllNeoplasia_controls"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: NeoplasiaType_control"))
print(PCoA1_2_extra)

plot1 <- ggMarginal(PCoA1_2_extra, type="boxplot", groupColour = TRUE, groupFill = TRUE)
print(plot1)

# a. PERMANOVAs ----
#univariable
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_a, covariates_adonis, "NeoplasiaType_control")

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_all_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:16]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_a <- adonisResults

adonis_plot_uni_Lynch_a <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_a)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_a <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_a[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_a)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_a <- results_df_excluded

adonis_plot_multi_Lynch_a <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_a)

# a. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_a) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_a <- merge(anova_results_df, permutest_results_df)

#### Beta-diversity pathways ----
Lynch_all_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_pathwaysfiltered.csv")
Lynch_all_pathways_filt <- as.data.frame(Lynch_all_pathways_filt)
rownames(Lynch_all_pathways_filt) <- Lynch_all_pathways_filt[, 1]
Lynch_all_pathways_filt <- Lynch_all_pathways_filt[, -1]
Lynch_all_pathways_filt_clr <- decostand(Lynch_all_pathways_filt, method = 'clr', pseudocount = min(Lynch_all_pathways_filt[Lynch_all_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_all_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# a. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_a <- list()

for (variable in totest_Lynch_a) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_a[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_a, ncol = 3, nrow = 3)

#Polished plot neoplasia types versus controls
library(ggExtra)

ref <- reformulate("NeoplasiaType_control", paste("cbind(", paste0("V", 1:5, collapse = ","), ")"))
centroids <- aggregate(ref, PCoA_meta, mean)
PCoA1_2_extra <- ggplot(PCoA_meta, aes_string(x = "V1", y = "V2", color = "NeoplasiaType_control")) + 
  xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
  ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = "AllNeoplasia_controls"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: NeoplasiaType_control"))
print(PCoA1_2_extra)

plot1 <- ggMarginal(PCoA1_2_extra, type="boxplot", groupColour = TRUE, groupFill = TRUE)
print(plot1)

# a. PERMANOVAs ----
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_a, covariates_adonis, "NeoplasiaType_control")

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_all_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:16]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_a <- adonisResults

adonis_plot_uni_Lynch_pwy_a <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_a)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_a <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_a[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_a)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_a <- results_df_excluded

adonis_plot_multi_Lynch_pwy_a <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_a)

# a. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_a) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_a <- merge(anova_results_df, permutest_results_df)






#### b.Beta-diversity species ----
Lynch_neoplasia_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_speciesfiltered.csv")
Lynch_neoplasia_species_filt <- as.data.frame(Lynch_neoplasia_species_filt)
rownames(Lynch_neoplasia_species_filt) <- Lynch_neoplasia_species_filt[, 1]
Lynch_neoplasia_species_filt <- Lynch_neoplasia_species_filt[, -1]
Lynch_neoplasia_species_filt_clr <- decostand(Lynch_neoplasia_species_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_species_filt[Lynch_neoplasia_species_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_neoplasia_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# b. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_b <- list()

for (variable in totest_Lynch_b) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_b[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_b, ncol = 1, nrow = 1)

# b. PERMANOVAs ----
#univariable
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_b, covariates_adonis)

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_neoplasia_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:7]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_b <- adonisResults

adonis_plot_uni_Lynch_b <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_b)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_b <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_b[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_b)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", "Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_b <- results_df_excluded

adonis_plot_multi_Lynch_b <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_b)

# b. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_b) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_b <- merge(anova_results_df, permutest_results_df)

#### b. Beta-diversity pathways ----
Lynch_neoplasia_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_pathwaysfiltered.csv")
Lynch_neoplasia_pathways_filt <- as.data.frame(Lynch_neoplasia_pathways_filt)
rownames(Lynch_neoplasia_pathways_filt) <- Lynch_neoplasia_pathways_filt[, 1]
Lynch_neoplasia_pathways_filt <- Lynch_neoplasia_pathways_filt[, -1]
Lynch_neoplasia_pathways_filt_clr <- decostand(Lynch_neoplasia_pathways_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_pathways_filt[Lynch_neoplasia_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_neoplasia_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# b. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_b <- list()

for (variable in totest_Lynch_b) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_b[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_b, ncol = 1, nrow = 2)

# b. PERMANOVAs ----
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_b, covariates_adonis)

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_neoplasia_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:7]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_b <- adonisResults

adonis_plot_uni_Lynch_pwy_b <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_b)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_b <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_b[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_b)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_b <- results_df_excluded

adonis_plot_multi_Lynch_pwy_b <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_b)

# b. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_b) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_b <- merge(anova_results_df, permutest_results_df)







#### c.Beta-diversity species ----
Lynch_control_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_speciesfiltered.csv")
Lynch_control_species_filt <- as.data.frame(Lynch_control_species_filt)
rownames(Lynch_control_species_filt) <- Lynch_control_species_filt[, 1]
Lynch_control_species_filt <- Lynch_control_species_filt[, -1]
Lynch_control_species_filt_clr <- decostand(Lynch_control_species_filt, method = 'clr', pseudocount = min(Lynch_control_species_filt[Lynch_control_species_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_control_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# c. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_c <- list()

for (variable in totest_Lynch_c) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_c[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_c, ncol = 2, nrow = 1)

# c. PERMANOVAs ----
#univariable
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_c, covariates_adonis)

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_control_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_c <- adonisResults

adonis_plot_uni_Lynch_c <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_c)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_c <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_c[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_c)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", "Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_c <- results_df_excluded

adonis_plot_multi_Lynch_c <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_c)

# c. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_c) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_c <- merge(anova_results_df, permutest_results_df)

#### c. Beta-diversity pathways ----
Lynch_control_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_pathwaysfiltered.csv")
Lynch_control_pathways_filt <- as.data.frame(Lynch_control_pathways_filt)
rownames(Lynch_control_pathways_filt) <- Lynch_control_pathways_filt[, 1]
Lynch_control_pathways_filt <- Lynch_control_pathways_filt[, -1]
Lynch_control_pathways_filt_clr <- decostand(Lynch_control_pathways_filt, method = 'clr', pseudocount = min(Lynch_control_pathways_filt[Lynch_control_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(Lynch_control_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# c. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, meta_Lynch_baseline, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_c <- list()

for (variable in totest_Lynch_c) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_c[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_c, ncol = 2, nrow = 1)

# c. PERMANOVAs ----
covariates_adonis <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")
adonisVarsTouse <- c(totest_Lynch_c, covariates_adonis)

adonisResults <- NULL
permNR <- 1000
set.seed(1456)

inMB <- Lynch_control_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- meta_Lynch_baseline
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_c <- adonisResults

adonis_plot_uni_Lynch_pwy_c <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_c)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_c <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_c[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_c)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bowel_Resection", "Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_c <- results_df_excluded

adonis_plot_multi_Lynch_pwy_c <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_c)

# c. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_Lynch_c) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_c <- merge(anova_results_df, permutest_results_df)







#### Saving plots and dataframes together ----
#species
PCOAs_species_Lynch <- c(plot_list_beta_Lynch_a, plot_list_beta_Lynch_b, plot_list_beta_Lynch_c)
ggarrange(plotlist = PCOAs_species_Lynch, ncol = 3, nrow = 4)

univariable_adonis_species_b <- univariable_adonis_species_b %>% filter(!Var %in% covariates_adonis)
univariable_adonis_species_c <- univariable_adonis_species_c %>% filter(!Var %in% covariates_adonis)
univariable_species_Lynch <- bind_rows(univariable_adonis_species_a, univariable_adonis_species_b, univariable_adonis_species_c)
write.csv(univariable_species_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Univariable_species_Lynch.csv")

adonis_plot_uni_Lynch_species <- ggplot(univariable_species_Lynch, aes(reorder(univariable_species_Lynch$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_species)

multivariable_species_Lynch <- bind_rows(multivariable_adonis_species_a, multivariable_adonis_species_b, multivariable_adonis_species_c)
write.csv(multivariable_species_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Multivariable_species_Lynch.csv")

adonis_plot_uni_Lynch_species_multi <- ggplot(multivariable_species_Lynch, aes(reorder(multivariable_species_Lynch$Var, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_species_multi)

ggarrange(adonis_plot_uni_Lynch_species, adonis_plot_uni_Lynch_species_multi, ncol = 2, nrow = 1)

betadisper_species_Lynch <- bind_rows(betadisper_results_df_a, betadisper_results_df_b, betadisper_results_df_c)
write.csv(betadisper_species_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Betadisper_species_Lynch.csv")

#pathways
PCOAs_pathways_Lynch <- c(plot_list_beta_Lynch_pwy_a, plot_list_beta_Lynch_pwy_b, plot_list_beta_Lynch_pwy_c)
ggarrange(plotlist = PCOAs_pathways_Lynch, ncol = 3, nrow = 4)

univariable_adonis_pathways_b <- univariable_adonis_pathways_b %>% filter(!Var %in% covariates_adonis)
univariable_adonis_pathways_c <- univariable_adonis_pathways_c %>% filter(!Var %in% covariates_adonis)
univariable_pathways_Lynch <- bind_rows(univariable_adonis_pathways_a, univariable_adonis_pathways_b, univariable_adonis_pathways_c)
write.csv(univariable_pathways_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Univariable_pathways_Lynch.csv")

adonis_plot_uni_Lynch_pathways <- ggplot(univariable_pathways_Lynch, aes(reorder(univariable_pathways_Lynch$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pathways)

multivariable_pathways_Lynch <- bind_rows(multivariable_adonis_pathways_a, multivariable_adonis_pathways_b, multivariable_adonis_pathways_c)
write.csv(multivariable_pathways_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Multivariable_pathways_Lynch.csv")

adonis_plot_uni_Lynch_pathways_multi <- ggplot(multivariable_pathways_Lynch, aes(reorder(multivariable_pathways_Lynch$Var, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pathways_multi)

ggarrange(adonis_plot_uni_Lynch_pathways, adonis_plot_uni_Lynch_pathways_multi, ncol = 2, nrow = 1)

betadisper_pathways_Lynch <- bind_rows(betadisper_results_df_pwy_a, betadisper_results_df_pwy_b, betadisper_results_df_pwy_c)
write.csv(betadisper_pathways_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Betadisper_pathways_Lynch.csv")



#### d.Beta-diversity species ----
LynchGP_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_speciesfiltered.csv")
LynchGP_species_filt <- as.data.frame(LynchGP_species_filt)
rownames(LynchGP_species_filt) <- LynchGP_species_filt[, 1]
LynchGP_species_filt <- LynchGP_species_filt[, -1]
LynchGP_species_filt_clr <- decostand(LynchGP_species_filt, method = 'clr', pseudocount = min(LynchGP_species_filt[LynchGP_species_filt > 0])/2)

LynchGP_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_metadata.csv")
LynchGP_meta <- as.data.frame(LynchGP_meta)
rownames(LynchGP_meta) <- LynchGP_meta$Participant_ID
LynchGP_meta <- LynchGP_meta[, -1]

# Calculate aitchison distance
vegdist(LynchGP_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# d. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchGP_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_d <- list()

for (variable in totest_LynchGP_d) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_d[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_d, ncol = 3, nrow = 1)

#Polished plot neoplasia types versus controls
library(ggExtra)

ref <- reformulate("cohort", paste("cbind(", paste0("V", 1:5, collapse = ","), ")"))
centroids <- aggregate(ref, PCoA_meta, mean)
PCoA1_2_extra <- ggplot(PCoA_meta, aes_string(x = "V1", y = "V2", color = "cohort")) + 
  xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
  ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = "cohort"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: cohort"))
print(PCoA1_2_extra)

plot1 <- ggMarginal(PCoA1_2_extra, type="boxplot", groupColour = TRUE, groupFill = TRUE)
print(plot1)

# d. PERMANOVAs ----
#univariable
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_d, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchGP_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchGP_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_d <- adonisResults

adonis_plot_uni_Lynch_d <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_d)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_d <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_d[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_d)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_d <- results_df_excluded

adonis_plot_multi_Lynch_d <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_d)

# d. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_d) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_d <- merge(anova_results_df, permutest_results_df)

#### d. Beta-diversity pathways ----
LynchGP_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_pathwaysfiltered.csv")
LynchGP_pathways_filt <- as.data.frame(LynchGP_pathways_filt)
rownames(LynchGP_pathways_filt) <- LynchGP_pathways_filt[, 1]
LynchGP_pathways_filt <- LynchGP_pathways_filt[, -1]
LynchGP_pathways_filt_clr <- decostand(LynchGP_pathways_filt, method = 'clr', pseudocount = min(LynchGP_pathways_filt[LynchGP_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(LynchGP_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# d. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchGP_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_d <- list()

for (variable in totest_LynchGP_d) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_d[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_d, ncol = 3, nrow = 1)

# d. PERMANOVAs ----
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_d, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchGP_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchGP_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_d <- adonisResults

adonis_plot_uni_Lynch_pwy_d <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_d)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_d <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_d[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_d)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_d <- results_df_excluded

adonis_plot_multi_Lynch_pwy_d<- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_d)

# d. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_d) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_d <- merge(anova_results_df, permutest_results_df)

#### e.Beta-diversity species ----
LynchLL_neoplasia_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_speciesfiltered.csv")
LynchLL_neoplasia_species_filt <- as.data.frame(LynchLL_neoplasia_species_filt)
rownames(LynchLL_neoplasia_species_filt) <- LynchLL_neoplasia_species_filt[, 1]
LynchLL_neoplasia_species_filt <- LynchLL_neoplasia_species_filt[, -1]
LynchLL_neoplasia_species_filt_clr <- decostand(LynchLL_neoplasia_species_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_species_filt[LynchLL_neoplasia_species_filt > 0])/2)

LynchLL_neoplasia_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_metadata.csv")
LynchLL_neoplasia_meta <- as.data.frame(LynchLL_neoplasia_meta)
rownames(LynchLL_neoplasia_meta) <- LynchLL_neoplasia_meta$Participant_ID
LynchLL_neoplasia_meta <- LynchLL_neoplasia_meta[, -1]

# Calculate aitchison distance
vegdist(LynchLL_neoplasia_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# e. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchLL_neoplasia_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_e <- list()

for (variable in totest_LynchGP_e) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_e[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_e, ncol = 2, nrow = 2)

# e. PERMANOVAs ----
#univariable
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_e, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchLL_neoplasia_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchLL_neoplasia_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:9]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_e <- adonisResults

adonis_plot_uni_Lynch_e <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_e)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_e <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_e[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_e)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_e <- results_df_excluded

adonis_plot_multi_Lynch_e <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_e)

# e. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_e) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_e <- merge(anova_results_df, permutest_results_df)

#### e. Beta-diversity pathways ----
LynchLL_neoplasia_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_pathwaysfiltered.csv")
LynchLL_neoplasia_pathways_filt <- as.data.frame(LynchLL_neoplasia_pathways_filt)
rownames(LynchLL_neoplasia_pathways_filt) <- LynchLL_neoplasia_pathways_filt[, 1]
LynchLL_neoplasia_pathways_filt <- LynchLL_neoplasia_pathways_filt[, -1]
LynchLL_neoplasia_pathways_filt_clr <- decostand(LynchLL_neoplasia_pathways_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_pathways_filt[LynchLL_neoplasia_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(LynchLL_neoplasia_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# e. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchLL_neoplasia_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_e <- list()

for (variable in totest_LynchGP_e) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_e[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_e, ncol = 2, nrow = 2)

# e. PERMANOVAs ----
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_e, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchLL_neoplasia_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchLL_neoplasia_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:9]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_e <- adonisResults

adonis_plot_uni_Lynch_pwy_e <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_e)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_e <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_e[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_e)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_e <- results_df_excluded

adonis_plot_multi_Lynch_pwy_e<- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_e)

# e. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_e) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_e <- merge(anova_results_df, permutest_results_df)

#### f.Beta-diversity species ----
LynchLL_controls_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_speciesfiltered.csv")
LynchLL_controls_species_filt <- as.data.frame(LynchLL_controls_species_filt)
rownames(LynchLL_controls_species_filt) <- LynchLL_controls_species_filt[, 1]
LynchLL_controls_species_filt <- LynchLL_controls_species_filt[, -1]
LynchLL_controls_species_filt_clr <- decostand(LynchLL_controls_species_filt, method = 'clr', pseudocount = min(LynchLL_controls_species_filt[LynchLL_controls_species_filt > 0])/2)

LynchLL_controls_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_metadata.csv")
LynchLL_controls_meta <- as.data.frame(LynchLL_controls_meta)
rownames(LynchLL_controls_meta) <- LynchLL_controls_meta$Participant_ID
LynchLL_controls_meta <- LynchLL_controls_meta[, -1]

# Calculate aitchison distance
vegdist(LynchLL_controls_species_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# f. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchLL_controls_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_f <- list()

for (variable in totest_LynchGP_f) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_e[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_e, ncol = 3, nrow = 1)

# f. PERMANOVAs ----
#univariable
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_f, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchLL_controls_species_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchLL_controls_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_species_f <- adonisResults

adonis_plot_uni_Lynch_f <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_f)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_f <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_f[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_f)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]

multivariable_adonis_species_f <- results_df_excluded

adonis_plot_multi_Lynch_f <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_f)

# f. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_f) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df_f <- merge(anova_results_df, permutest_results_df)

#### f. Beta-diversity pathways ----
LynchLL_controls_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_pathwaysfiltered.csv")
LynchLL_controls_pathways_filt <- as.data.frame(LynchLL_controls_pathways_filt)
rownames(LynchLL_controls_pathways_filt) <- LynchLL_controls_pathways_filt[, 1]
LynchLL_controls_pathways_filt <- LynchLL_controls_pathways_filt[, -1]
LynchLL_controls_pathways_filt_clr <- decostand(LynchLL_controls_pathways_filt, method = 'clr', pseudocount = min(LynchLL_controls_pathways_filt[LynchLL_controls_pathways_filt > 0])/2)

# Calculate aitchison distance
vegdist(LynchLL_controls_pathways_filt_clr, method = "euclidean") -> Beta_diversity #=Aitchison distance because CLR transformed

# f. Plot PCoAs ----
cmdscale(Beta_diversity, k=5, eig = TRUE) -> my_pcoa
PC = as.matrix(my_pcoa$points)
var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)
PCoA_meta <- NULL
PCoA_meta <- merge(PC, LynchLL_controls_meta, by.x = 'row.names', by.y = 'Participant_ID')

# Make PCOA plots
results <- list()
plot_list_beta_Lynch_pwy_f <- list()

for (variable in totest_LynchGP_f) {
  # Perform alpha analysis
  results[[variable]] <- perform_beta_analysis_Lynch(variable)
  # Store the alpha diversity plot in a list
  plot_list_beta_Lynch_pwy_f[[variable]] <- results[[variable]]$beta_diversity
}

ggarrange(plotlist = plot_list_beta_Lynch_pwy_f, ncol = 3, nrow = 1)

# f. PERMANOVAs ----
covariates_adonis_GP <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")
adonisVarsTouse <- c(totest_LynchGP_f, covariates_adonis_GP)

adonisResults <- NULL
permNR <- 100
set.seed(1456)

inMB <- LynchLL_controls_pathways_filt_clr
inMB <- rownames_to_column(inMB, "Participant_ID")
inPhenos <- LynchLL_controls_meta
inPhenos[sapply(inPhenos, is.character)] <- lapply(inPhenos[sapply(inPhenos, is.character)], as.factor)

for (i in adonisVarsTouse[1:8]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  #print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"Participant_ID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="Participant_ID",y=inMB,by.y="Participant_ID")
  rownames(allDF) <- allDF$Participant_ID
  allDF$Participant_ID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating Aitchison distance')
  inBC <- vegdist(allDF,method = "euclidean",parallel=4)
  #print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis2(inBC ~ av,permutations=permNR)
    model.frame(ad) -> aov.table
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov.table[1,1],
                         SumsOfSqs=aov.table[1,2],
                         FModel=aov.table[1,4],
                         R2=aov.table[1,3],
                         pval=aov.table[1,5],
                         FDR.BH=NA,
                         Significant_nominal=NA,
                         Significant_FDR=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    print (paste0('--- ',i,' DONE! ---'))
  }
}

rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant_nominal <- ifelse(adonisResults$pval < 0.05, "Yes", "No")
adonisResults$Significant_FDR <- ifelse(adonisResults$FDR.BH < 0.05, "Yes", "No")
adonisResults$significance <- ifelse(adonisResults$Significant_FDR == "Yes", "FDR", ifelse(adonisResults$Significant_nominal == "Yes", "Nominal", "No"))
rownames(adonisResults) <- NULL
univariable_adonis_pathways_f <- adonisResults

adonis_plot_uni_Lynch_pwy_f <- ggplot(adonisResults, aes(reorder(adonisResults$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_pwy_f)

# PERMANOVA multivariable
results_list <- list()
plot_list_adonis_all_pwy_f <- list()
rownames(inMB) <- inMB$Participant_ID

# Loop through variables #takes a long time
for (i in adonisVarsTouse) {
  print(paste(' >>> ANALYSING VARIABLE <', i, '>    <<<'))
  # Perform Adonis analysis
  analysis_results <- perform_multiadonis_analysis_pathways_Lynch(variable = i, meta_data = inPhenos, clr_data = inMB, covariates = covariates_adonis_GP)
  # Store results and plot in lists
  results_list[[i]] <- analysis_results$results
  plot_list_adonis_all_pwy_f[[i]] <- analysis_results$plot
  print(paste0('--- ', i, ' DONE! ---'))
}

# Combine all results into one dataframe and print plots
results_df_all <- do.call(rbind, results_list)
ggarrange(plotlist = plot_list_adonis_all_pwy_f)

excluded_suffixes <- c(".Sex", ".Age", ".BMI", ".Smoking", ".Bristol_score")
rows_to_exclude <- grepl(paste(excluded_suffixes, collapse = "|"), rownames(results_df_all))
results_df_excluded <- results_df_all[!rows_to_exclude, ]
multivariable_adonis_pathways_f <- results_df_excluded

adonis_plot_multi_Lynch_pwy_f <- ggplot(results_df_excluded, aes(reorder(results_df_excluded$Variable, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_multi_Lynch_pwy_f)

# f. Betadisper for homogeneity of group dispersion ----
anova_results_df <- data.frame()
permutest_results_df <- data.frame()
betadisper_results <- list()
anova_results <- list()
permutest_results <- list()
plots <- list()

# Loop through each variable in totest_Lynch
for (variable in totest_LynchGP_f) {
  
  # Perform betadisper
  mod <- betadisper(Beta_diversity, PCoA_meta[[variable]])
  betadisper_results[[variable]] <- mod
  
  # Perform ANOVA
  anova_res <- anova(mod)
  anova_results[[variable]] <- anova_res
  
  # Perform permutest
  pmod <- permutest(mod, permutations = 99, pairwise = TRUE)
  permutest_results[[variable]] <- pmod
  
  # Extract ANOVA results and add to dataframe if valid
  if (!is.null(anova_res) && nrow(anova_res) > 0) {
    anova_results_df <- rbind(anova_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(anova_res$Df[1]), anova_res$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(anova_res$`Sum Sq`[1]), anova_res$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(anova_res$`Mean Sq`[1]), anova_res$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(anova_res$F[1]), anova_res$F[1], NA),
      PValue = ifelse(!is.null(anova_res$`Pr(>F)`[1]), anova_res$`Pr(>F)`[1], NA)
    ))
  }
  
  # Extract permutest results and add to dataframe if valid
  if (!is.null(pmod$tab) && nrow(pmod$tab) > 0) {
    permutest_results_df <- rbind(permutest_results_df, data.frame(
      Variable = variable,
      Df = ifelse(!is.null(pmod$tab$Df[1]), pmod$tab$Df[1], NA),
      SumsOfSqs = ifelse(!is.null(pmod$tab$`Sum Sq`[1]), pmod$tab$`Sum Sq`[1], NA),
      MeanSqs = ifelse(!is.null(pmod$tab$`Mean Sq`[1]), pmod$tab$`Mean Sq`[1], NA),
      FValue = ifelse(!is.null(pmod$tab$`F`[1]), pmod$tab$`F`[1], NA),
      PValue = ifelse(!is.null(pmod$tab$`Pr(>F)`[1]), pmod$tab$`Pr(>F)`[1], NA),
      Permutations = ifelse(!is.null(pmod$tab$`N.Perm`[1]), pmod$tab$`N.Perm`[1], NA)
    ))
  }
  
  # Plot the results and store the plot
  plot(mod, main = paste("Beta Dispersion for", variable))
  plots[[variable]] <- recordPlot()
}

# Print or inspect results for each variable
betadisper_results
colnames(anova_results_df)[colnames(anova_results_df) == "PValue"] <- "anova_pvalue"
colnames(permutest_results_df)[colnames(permutest_results_df) == "PValue"] <- "permutest_pvalue"
betadisper_results_df <- merge(anova_results_df, permutest_results_df)
betadisper_results_df_pwy_f <- merge(anova_results_df, permutest_results_df)
#### Saving plots and dataframes together ----
#species
PCOAs_species_Lynch_LL <- c(plot_list_beta_Lynch_d, plot_list_beta_Lynch_e, plot_list_beta_Lynch_f)
ggarrange(plotlist = PCOAs_species_Lynch_LL, ncol = 4, nrow = 3)

PCOAs_species_Lynch_LL_controlcases <- c(plot_list_beta_Lynch_e, plot_list_beta_Lynch_f)
ggarrange(plotlist = PCOAs_species_Lynch_LL_controlcases, ncol = 4, nrow = 2)

univariable_adonis_species_e$Var <- ifelse(
  univariable_adonis_species_e$Var %in% covariates_adonis_GP,
  paste0(univariable_adonis_species_e$Var, "_cases"),
  univariable_adonis_species_e$Var)

univariable_adonis_species_f$Var <- ifelse(
  univariable_adonis_species_f$Var %in% covariates_adonis_GP,
  paste0(univariable_adonis_species_f$Var, "_controls"),
  univariable_adonis_species_f$Var)

univariable_species_Lynch_LL <- bind_rows(univariable_adonis_species_d, univariable_adonis_species_e, univariable_adonis_species_f)
write.csv(univariable_species_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Univariable_species_LynchLifeLines.csv")

adonis_plot_uni_Lynch_LL_species <- ggplot(univariable_species_Lynch_LL, aes(reorder(univariable_species_Lynch_LL$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_LL_species)

multivariable_species_Lynch_LL <- bind_rows(multivariable_adonis_species_d, multivariable_adonis_species_e, multivariable_adonis_species_f)
write.csv(multivariable_species_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Multivariable_species_LynchLifeLines_allcovariates.csv")

adonis_plot_uni_Lynch_LL_species_multi <- ggplot(multivariable_species_Lynch_LL, aes(reorder(multivariable_species_Lynch_LL$Var, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: species Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_LL_species_multi)

ggarrange(adonis_plot_uni_Lynch_LL_species, adonis_plot_uni_Lynch_LL_species_multi, ncol = 2, nrow = 1)

betadisper_species_Lynch_LL <- bind_rows(betadisper_results_df_d, betadisper_results_df_e, betadisper_results_df_f)
write.csv(betadisper_species_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Betadisper_species_LynchLifeLines.csv")

#pathways
PCOAs_pathways_Lynch_LL <- c(plot_list_beta_Lynch_pwy_d, plot_list_beta_Lynch_pwy_e, plot_list_beta_Lynch_pwy_f)
ggarrange(plotlist = PCOAs_pathways_Lynch_LL, ncol = 4, nrow = 3)

PCOAs_pathways_Lynch_LL_casescontrols <- c(plot_list_beta_Lynch_pwy_e, plot_list_beta_Lynch_pwy_f)
ggarrange(plotlist = PCOAs_pathways_Lynch_LL_casescontrols, ncol = 4, nrow = 2)

univariable_adonis_pathways_e$Var <- ifelse(
  univariable_adonis_pathways_e$Var %in% covariates_adonis_GP,
  paste0(univariable_adonis_pathways_e$Var, "_cases"),
  univariable_adonis_pathways_e$Var)

univariable_adonis_pathways_f$Var <- ifelse(
  univariable_adonis_pathways_f$Var %in% covariates_adonis_GP,
  paste0(univariable_adonis_pathways_f$Var, "_controls"),
  univariable_adonis_pathways_f$Var)

univariable_pathways_Lynch_LL <- bind_rows(univariable_adonis_pathways_d, univariable_adonis_pathways_e, univariable_adonis_pathways_f)
write.csv(univariable_pathways_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Univariable_pathways_LynchLifeLines.csv")

adonis_plot_uni_Lynch_LL_pathways <- ggplot(univariable_pathways_Lynch_LL, aes(reorder(univariable_pathways_Lynch_LL$Var, R2), R2, fill=significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Univariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_LL_pathways)

multivariable_pathways_Lynch_LL <- bind_rows(multivariable_adonis_pathways_d, multivariable_adonis_pathways_e, multivariable_adonis_pathways_f)
write.csv(multivariable_pathways_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Multivariable_pathways_LynchLifeLines_allcovariates.csv")

adonis_plot_uni_Lynch_LL_pathways_multi <- ggplot(multivariable_pathways_Lynch_LL, aes(reorder(multivariable_pathways_Lynch_LL$Var, R2), R2, fill=Significance)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + ggtitle("Multivariable: pathways Aitchison distance") +
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11)) +
  scale_fill_manual(values = adonis_colors) +  theme(text = element_text(size = 11), legend.position = "bottom") 
plot(adonis_plot_uni_Lynch_LL_pathways_multi)

ggarrange(adonis_plot_uni_Lynch_LL_pathways, adonis_plot_uni_Lynch_LL_pathways_multi, ncol = 2, nrow = 1)

betadisper_pathways_Lynch_LL <- bind_rows(betadisper_results_df_pwy_d, betadisper_results_df_pwy_e, betadisper_results_df_pwy_f)
write.csv(betadisper_pathways_Lynch_LL, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Beta_diversity/Betadisper_pathways_LynchLifeLines.csv")


#### Extra tests ----
PCoA1_2_extra <- ggplot(PCoA_meta, aes_string(x = "V1", y = "V2", color = "NeoplasiaType_control")) + 
  xlab(paste0("PCo1: ", var_expl[1], "% variance")) +
  ylab(paste0("PCo2: ", var_expl[2], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V1", y = "V2"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V1", y = "V2", col = "AllNeoplasia_controls"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: NeoplasiaType_control"))
print(PCoA1_2_extra)

PCoA2_3_extra <- ggplot(PCoA_meta, aes_string(x = "V2", y = "V3", color = "NeoplasiaType_control")) + 
  xlab(paste0("PCo2: ", var_expl[2], "% variance")) +
  ylab(paste0("PCo3: ", var_expl[3], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V2", y = "V3"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V2", y = "V3", col = "AllNeoplasia_controls"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: NeoplasiaType_control"))
print(PCoA2_3_extra)

PCoA3_4_extra <- ggplot(PCoA_meta, aes_string(x = "V3", y = "V4", color = "NeoplasiaType_control")) + 
  xlab(paste0("PCo3: ", var_expl[3], "% variance")) +
  ylab(paste0("PCo4: ", var_expl[4], "% variance")) +
  geom_point(size = 3, alpha = 0.9) + 
  theme_light() + 
  stat_ellipse(geom = "polygon", alpha = 0.05, linetype = "dotted") +
  #stat_ellipse(geom = "polygon", alpha = 0.05) +
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold", color = "black"),
        legend.position = "bottom",
        strip.text = element_text(size = 12)) +
  scale_fill_manual(values = my_colours) +
  scale_color_manual(values = outline_colors) +
  geom_point(data = centroids, shape = 16, stroke = 3, size = 4, aes_string(x = "V3", y = "V4"), alpha = 1) +
  geom_point(data = centroids, shape = 21, stroke = 2, size = 4, aes_string(x = "V3", y = "V4", col = "AllNeoplasia_controls"), alpha = 0.8, colour = "black") +
  ggtitle(paste("Beta Diversity: NeoplasiaType_control"))
print(PCoA3_4_extra)

ggarrange(PCoA1_2_extra, PCoA2_3_extra, PCoA3_4_extra, ncol = 3, nrow = 1)

#Ordination
dist_pcoa <- pcoa(dist(pathways_Lynch_filt_clr))
biplot.pcoa(dist_pcoa, pathways_Lynch_filt_clr, dir.axis2 = -1)

#test NDMS
Beta_diversity %>%
  metaMDS(trace = F) %>%
  ordiplot(type = "none") %>%
  text("sites")

#Test PCA
PCA <- rda(pathways_Lynch_filt_clr, scale = FALSE)
barplot(as.vector(PCA$CA$eig)/sum(PCA$CA$eig)) 
sum((as.vector(PCA$CA$eig)/sum(PCA$CA$eig))[1:2])
plot(PCA)
test_object <- biplot(PCA, choices = c(1,2), type = c("text", "points"), xlim = c(-5,10))

