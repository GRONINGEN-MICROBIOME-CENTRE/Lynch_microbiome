# LYNCH SYNDROME PROJECT, COLLABORATION AMSTERDAM GRONINGEN
# by Femke Prins, Q2 2024
# script for performing alpha-diversity analyses

setwd("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data")

# Load values ----
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
                    "AdvAdenomasCRCAdvSerr" = "#cb75a1",
                    "Low-risk"= "#21560a", "High-risk"= "#0b5394", 
                    "MLH1"="#00BFC4", "MSH2/EPCAM"="#719FA6", 
                    "MSH6" = "#7CAE00", "PMS2" = "#B6D7A8",
                    "General_population" = "orange", "Lynch" = "darkblue",
                    "LifeLines" = "darkorange")
my_colours <- outline_colors

# Load functions ----
calculate_alpha <- function(inDF,IDcol="RowNames",
                            metrics=c("shannon","simpson","invsimpson","richness"),
                            DIVlvls=c("taxS")) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over metrics, calculate each
  # NOTE: richness is not implemented in vegan, requires special treatment
  for (l in DIVlvls) {
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l) } }
  for (m in metrics) {
    print(paste0('  > calculating ',m,'[',l,']'))
    if (m=="richness") {
      inDFpa <- inDF
      inDFpa[inDFpa > 0] <- 1
      dv <- rowSums(inDFpa)
    } else {
      dv <- diversity(inDF,index = m)
    }
    DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
    colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
  }
  return(DIVMatrix)
} #made by Ranko

perform_alpha_analysis <- function(variable, Data_div_total, covariates = NULL) {
  # Filter out NA rows in the variable being tested
  div_complete <- Data_div_total[!is.na(Data_div_total[[variable]]), ]
  
  # Make a plot for Alpha diversity
  levels <- as.character(unique(div_complete[[variable]]))
  levels_without_na <- levels[!is.na(levels)]
  my_comparisons_all <- combn(levels_without_na, 2, simplify = FALSE)
  
  alpha_diversity <- ggplot(div_complete, aes_string(x = variable, y = "DIV.S.shannon", fill = variable, color = variable)) +
    geom_jitter(alpha = 1, width = 0.1) +
    geom_violin(trim = FALSE, position = position_dodge(0.9), alpha = 0.5, linewidth = 0.8) +
    geom_boxplot(alpha = 0.3, width = 0.3, size = 0.8) +
    theme_light() +
    theme(
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 0),
      legend.position = "none",
      strip.text = element_text(size = 12)) +
    ylab("Shannon diversity index") +
    scale_fill_manual(values = my_colours) +
    scale_color_manual(values = outline_colors) +
    stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", paired = FALSE, hide.ns = FALSE) +
    ggtitle(paste("Alpha Diversity:", variable)) 
  
  # Perform linear model analyses
  if (!is.null(covariates)) {
    # Remove covariates with only one level
    valid_covariates <- covariates[sapply(covariates, function(cov) length(unique(div_complete[[cov]])) > 1)]
    
    if (length(valid_covariates) > 0) {
      # Create the formula with valid covariates
      formula_str <- paste("DIV.S.shannon ~", variable, "+", paste(valid_covariates, collapse = ' + '))
    } else {
      # If no valid covariates, only include the variable
      formula_str <- paste("DIV.S.shannon ~", variable)
    }
    
    formula <- as.formula(formula_str)
    lm_model <- lm(formula, data = Data_div_total)
    lm_summary <- summary(lm_model)
    lm_results <- data.frame(
      Variable = variable,
      Factor = rownames(lm_summary$coefficients),
      Coefficients = lm_summary$coefficients[, "Estimate"],
      P_values = lm_summary$coefficients[, "Pr(>|t|)"],
      Formula = formula_str  # Add the formula used to the dataframe
    )
  } else {
    lm_results <- NULL
  }
  
  # Return results
  return(list(alpha_diversity = alpha_diversity, lm_results = lm_results))
}

richness_pwys <- function(variable, Data_div_total, covariates = NULL) {
  # Filter out NA rows in the variable being tested
  div_complete <- Data_div_total[!is.na(Data_div_total[[variable]]), ]
  
  # Make a plot for Alpha diversity
  levels <- as.character(unique(div_complete[[variable]]))
  levels_without_na <- levels[!is.na(levels)]
  my_comparisons_all <- combn(levels_without_na, 2, simplify = FALSE)
  
  alpha_diversity <- ggplot(div_complete, aes_string(x = variable, y = "DIV.S.richness", fill = variable, color = variable)) +
    geom_jitter(alpha = 1, width = 0.1) +
    geom_violin(trim = FALSE, position = position_dodge(0.9), alpha = 0.5, linewidth = 0.8) +
    geom_boxplot(alpha = 0.3, width = 0.3, size = 0.8) +
    theme_light() +
    theme(
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black", hjust = 0),
      legend.position = "none",
      strip.text = element_text(size = 12)) +
    ylab("Richness") +
    scale_fill_manual(values = my_colours) +
    scale_color_manual(values = outline_colors) +
    stat_compare_means(comparisons = my_comparisons_all, method = "wilcox.test", paired = FALSE, hide.ns = FALSE) +
    ggtitle(paste("Richness:", variable)) 
  
  # Perform linear model analyses
  if (!is.null(covariates)) {
    # Remove covariates with only one level
    valid_covariates <- covariates[sapply(covariates, function(cov) length(unique(div_complete[[cov]])) > 1)]
    
    if (length(valid_covariates) > 0) {
      # Create the formula with valid covariates
      formula_str <- paste("DIV.S.richness ~", variable, "+", paste(valid_covariates, collapse = ' + '))
    } else {
      # If no valid covariates, only include the variable
      formula_str <- paste("DIV.S.richness ~", variable)
    }
    
    formula <- as.formula(formula_str)
    lm_model <- lm(formula, data = Data_div_total)
    lm_summary <- summary(lm_model)
    lm_results <- data.frame(
      Variable = variable,
      Factor = rownames(lm_summary$coefficients),
      Coefficients = lm_summary$coefficients[, "Estimate"],
      P_values = lm_summary$coefficients[, "Pr(>|t|)"],
      Formula = formula_str  # Add the formula used to the dataframe
    )
  } else {
    lm_results <- NULL
  }
  
  # Return results
  return(list(alpha_diversity = alpha_diversity, lm_results = lm_results))
}

#### a. Alpha-diversity: species Shannon index----
Lynch_all_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_species.csv")
Lynch_all_species <- as.data.frame(Lynch_all_species)
rownames(Lynch_all_species) <- Lynch_all_species[, 1]
Lynch_all_species <- Lynch_all_species[, -1]

Alpha_metrices <- calculate_alpha(Lynch_all_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_a <- list()
all_lm_results_Lynch_a <- data.frame()

for (variable in totest_Lynch_a) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_a[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_a <- rbind(all_lm_results_Lynch_a, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_a$P_values <- format(all_lm_results_Lynch_a$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_a) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_a, ncol = 3, nrow = 3)

#### b. Alpha-diversity: species Shannon index----
Lynch_neoplasia_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_species.csv")
Lynch_neoplasia_species <- as.data.frame(Lynch_neoplasia_species)
rownames(Lynch_neoplasia_species) <- Lynch_neoplasia_species[, 1]
Lynch_neoplasia_species <- Lynch_neoplasia_species[, -1]

meta_Lynch_baseline <- meta_Lynch_baseline %>%
  mutate(AdvAdenomasCRCAdvSerr_nonAdvAdenomas = case_when(
    `Type of lesion` %in% c(0, 1, 2) ~ "AdvAdenomasCRCAdvSerr",
    `Type of lesion` == 3 ~ "NonAdv_adenoma",
    TRUE ~ NA_character_))
table(meta_Lynch_baseline$AdvAdenomasCRCAdvSerr_nonAdvAdenomas)

Alpha_metrices <- calculate_alpha(Lynch_neoplasia_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_b <- list()
all_lm_results_Lynch_b <- data.frame()

for (variable in totest_Lynch_b) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_b[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_b <- rbind(all_lm_results_Lynch_b, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_b$P_values <- format(all_lm_results_Lynch_b$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_b) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_b, ncol = 1, nrow = 2)

#### c. Alpha-diversity: species Shannon index----
Lynch_control_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_species.csv")
Lynch_control_species <- as.data.frame(Lynch_control_species)
rownames(Lynch_control_species) <- Lynch_control_species[, 1]
Lynch_control_species <- Lynch_control_species[, -1]

Alpha_metrices <- calculate_alpha(Lynch_control_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_c <- list()
all_lm_results_Lynch_c <- data.frame()

for (variable in totest_Lynch_c) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_c[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_c <- rbind(all_lm_results_Lynch_c, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_c$P_values <- format(all_lm_results_Lynch_c$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_c) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_c, ncol = 2, nrow = 1)

modelx <- lm(DIV.S.shannon ~ NoNeo_Mutation + Sex + reads + Age + BMI + Smoking + Bowel_Resection + Bristol_score, data = Data_div_total)
summary(modelx)

emm <- emmeans(modelx, ~ NoNeo_Mutation)

# Show pairwise comparisons between levels of NoNeo_Mutation
pairwise_differences <- pairs(emm)
pairwise_differences

#### saving linear models and plotting Lynch alpha diversity ----
all_lm_results_Lynch <- bind_rows(all_lm_results_Lynch_a, all_lm_results_Lynch_b, all_lm_results_Lynch_c)
write.csv(all_lm_results_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Alpha_diversity/Shannon_linearmodel.csv")

all_plot_list_alpha_Lynch <- c(plot_list_alpha_Lynch_a, plot_list_alpha_Lynch_b, plot_list_alpha_Lynch_c)
ggarrange(plotlist = all_plot_list_alpha_Lynch, ncol = 3, nrow = 4)

#### a. Alpha-diversity: pathways richness ----
Lynch_all_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_pathways.csv")
Lynch_all_pathways <- as.data.frame(Lynch_all_pathways)
rownames(Lynch_all_pathways) <- Lynch_all_pathways[, 1]
Lynch_all_pathways <- Lynch_all_pathways[, -1]

Alpha_metrices <- calculate_alpha(Lynch_all_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_a <- list()
all_lm_richness_Lynch_a <- data.frame()

for (variable in totest_Lynch_a) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_a[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_a <- rbind(all_lm_richness_Lynch_a, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_a$P_values <- format(all_lm_richness_Lynch_a$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_a) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_a, ncol = 3, nrow = 3)

#### b. Alpha-diversity: pathways richness ----
Lynch_neoplasia_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_pathways.csv")
Lynch_neoplasia_pathways <- as.data.frame(Lynch_neoplasia_pathways)
rownames(Lynch_neoplasia_pathways) <- Lynch_neoplasia_pathways[, 1]
Lynch_neoplasia_pathways <- Lynch_neoplasia_pathways[, -1]

Alpha_metrices <- calculate_alpha(Lynch_neoplasia_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_b <- list()
all_lm_richness_Lynch_b <- data.frame()

for (variable in totest_Lynch_b) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_b[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_b <- rbind(all_lm_richness_Lynch_b, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_b$P_values <- format(all_lm_richness_Lynch_b$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_b) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_b, ncol = 1, nrow = 2)

#### c. Alpha-diversity: pathways richness ----
Lynch_control_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_pathways.csv")
Lynch_control_pathways <- as.data.frame(Lynch_control_pathways)
rownames(Lynch_control_pathways) <- Lynch_control_pathways[, 1]
Lynch_control_pathways <- Lynch_control_pathways[, -1]

Alpha_metrices <- calculate_alpha(Lynch_control_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, meta_Lynch_baseline, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_c <- list()
all_lm_richness_Lynch_c <- data.frame()

for (variable in totest_Lynch_c) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_Lynch)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_c[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_c <- rbind(all_lm_richness_Lynch_c, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_c$P_values <- format(all_lm_richness_Lynch_c$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_c) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_c, ncol = 2, nrow = 1)


modelx <- lm(DIV.S.richness ~ NoNeo_Mutation + Sex + reads + Age + BMI + Smoking + Bowel_Resection + Bristol_score, data = Data_div_total)
summary(modelx)
emm <- emmeans(modelx, ~ NoNeo_Mutation)

# Show pairwise comparisons between levels of NoNeo_Mutation
pairwise_differences <- pairs(emm)
pairwise_differences

#### saving linear models and plotting Lynch richness ----
all_lm_richness_Lynch <- bind_rows(all_lm_richness_Lynch_a, all_lm_richness_Lynch_b, all_lm_richness_Lynch_c)
write.csv(all_lm_richness_Lynch, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Alpha_diversity/Richness_linearmodel.csv")

all_plot_list_richness_Lynch <- c(plot_list_richness_Lynch_a, plot_list_richness_Lynch_b, plot_list_richness_Lynch_c)
ggarrange(plotlist = all_plot_list_richness_Lynch, ncol = 3, nrow = 4)

#### d. Alpha-diversity: species Shannon index----
LynchGP_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_species.csv")
LynchGP_species <- as.data.frame(LynchGP_species)
rownames(LynchGP_species) <- LynchGP_species[, 1]
LynchGP_species <- LynchGP_species[, -1]

LynchGP_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_metadata.csv")
LynchGP_meta <- as.data.frame(LynchGP_meta)
rownames(LynchGP_meta) <- LynchGP_meta$Participant_ID
LynchGP_meta <- LynchGP_meta[, -1]

#calculate shannon
Alpha_metrices <- calculate_alpha(LynchGP_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchGP_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_d <- list()
all_lm_results_Lynch_d <- data.frame()

for (variable in totest_LynchGP_d) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_d[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_d <- rbind(all_lm_results_Lynch_d, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_d$P_values <- format(all_lm_results_Lynch_d$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_d) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_d, ncol = 3, nrow = 1)

#### d. Alpha-diversity: pathways richness ----
LynchGP_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_pathways.csv")
LynchGP_pathways <- as.data.frame(LynchGP_pathways)
rownames(LynchGP_pathways) <- LynchGP_pathways[, 1]
LynchGP_pathways <- LynchGP_pathways[, -1]

Alpha_metrices <- calculate_alpha(LynchGP_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchGP_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_d <- list()
all_lm_richness_Lynch_d <- data.frame()

for (variable in totest_LynchGP_d) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_d[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_d <- rbind(all_lm_richness_Lynch_d, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_d$P_values <- format(all_lm_richness_Lynch_d$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_d) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_d, ncol = 3, nrow = 1)

#### e. Alpha-diversity: species Shannon index----
LynchLL_neoplasia_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_species.csv")
LynchLL_neoplasia_species <- as.data.frame(LynchLL_neoplasia_species)
rownames(LynchLL_neoplasia_species) <- LynchLL_neoplasia_species[, 1]
LynchLL_neoplasia_species <- LynchLL_neoplasia_species[, -1]

LynchLL_neoplasia_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_metadata.csv")
LynchLL_neoplasia_meta <- as.data.frame(LynchLL_neoplasia_meta)
rownames(LynchLL_neoplasia_meta) <- LynchLL_neoplasia_meta$Participant_ID
LynchLL_neoplasia_meta <- LynchLL_neoplasia_meta[, -1]

#calculate shannon
Alpha_metrices <- calculate_alpha(LynchLL_neoplasia_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchLL_neoplasia_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_e <- list()
all_lm_results_Lynch_e <- data.frame()

for (variable in totest_LynchGP_e) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_e[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_e <- rbind(all_lm_results_Lynch_e, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_e$P_values <- format(all_lm_results_Lynch_e$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_e) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_e, ncol = 2, nrow = 2)

#### e. Alpha-diversity: pathways richness ----
LynchLL_neoplasia_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_pathways.csv")
LynchLL_neoplasia_pathways <- as.data.frame(LynchLL_neoplasia_pathways)
rownames(LynchLL_neoplasia_pathways) <- LynchLL_neoplasia_pathways[, 1]
LynchLL_neoplasia_pathways <- LynchLL_neoplasia_pathways[, -1]

Alpha_metrices <- calculate_alpha(LynchLL_neoplasia_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchLL_neoplasia_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_e <- list()
all_lm_richness_Lynch_e <- data.frame()

for (variable in totest_LynchGP_e) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_e[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_e <- rbind(all_lm_richness_Lynch_e, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_e$P_values <- format(all_lm_richness_Lynch_e$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_e) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_e, ncol = 2, nrow = 2)

#### f. Alpha-diversity: species Shannon index----
LynchLL_controls_species <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_species.csv")
LynchLL_controls_species <- as.data.frame(LynchLL_controls_species)
rownames(LynchLL_controls_species) <- LynchLL_controls_species[, 1]
LynchLL_controls_species <- LynchLL_controls_species[, -1]

LynchLL_controls_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_metadata.csv")
LynchLL_controls_meta <- as.data.frame(LynchLL_controls_meta)
rownames(LynchLL_controls_meta) <- LynchLL_controls_meta$Participant_ID
LynchLL_controls_meta <- LynchLL_controls_meta[, -1]

#calculate shannon
Alpha_metrices <- calculate_alpha(LynchLL_controls_species, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchLL_controls_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)

results <- list()
plot_list_alpha_Lynch_f <- list()
all_lm_results_Lynch_f <- data.frame()

for (variable in totest_LynchGP_f) {
  # Perform alpha analysis
  results[[variable]] <- perform_alpha_analysis(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_alpha_Lynch_f[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_results_Lynch_f <- rbind(all_lm_results_Lynch_f, results[[variable]]$lm_results)
  }
}

all_lm_results_Lynch_f$P_values <- format(all_lm_results_Lynch_f$P_values, scientific = FALSE)
rownames(all_lm_results_Lynch_f) <- NULL
ggarrange(plotlist = plot_list_alpha_Lynch_f, ncol = 3, nrow = 1)

#### f. Alpha-diversity: pathways richness ----
LynchLL_controls_pathways <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_pathways.csv")
LynchLL_controls_pathways <- as.data.frame(LynchLL_controls_pathways)
rownames(LynchLL_controls_pathways) <- LynchLL_controls_pathways[, 1]
LynchLL_controls_pathways <- LynchLL_controls_pathways[, -1]

Alpha_metrices <- calculate_alpha(LynchLL_controls_pathways, IDcol = "RowNames", metrics = c("shannon", "simpson", "invsimpson", "richness"), DIVlvls = c("taxS"))
Data_div_total <- merge(Alpha_metrices, LynchLL_controls_meta, by.x = 'RN', by.y = 'Participant_ID')
rownames(Data_div_total) <- Data_div_total$RN
Data_div_total$RN <- NULL
Data_div_total[sapply(Data_div_total, is.character)] <- lapply(Data_div_total[sapply(Data_div_total, is.character)], as.factor)
Data_div_total <- subset(Data_div_total, DIV.S.richness >= 3) #remove outlier

results <- list()
plot_list_richness_Lynch_f <- list()
all_lm_richness_Lynch_f <- data.frame()

for (variable in totest_LynchGP_f) {
  # Perform alpha analysis
  results[[variable]] <- richness_pwys(variable, Data_div_total, covariates_LifeLines)
  # Store the alpha diversity plot in a list
  plot_list_richness_Lynch_f[[variable]] <- results[[variable]]$alpha_diversity
  # Get the linear model results
  if (!is.null(results[[variable]]$lm_results)) {
    all_lm_richness_Lynch_f <- rbind(all_lm_richness_Lynch_f, results[[variable]]$lm_results)
  }
}

all_lm_richness_Lynch_f$P_values <- format(all_lm_richness_Lynch_f$P_values, scientific = FALSE)
rownames(all_lm_richness_Lynch_f) <- NULL
ggarrange(plotlist = plot_list_richness_Lynch_f, ncol = 3, nrow = 1)

#### saving linear models and plotting Lynch + LifeLines alpha diversity ----
all_lm_results_Lynch_LifeLines <- bind_rows(all_lm_results_Lynch_d, all_lm_results_Lynch_e, all_lm_results_Lynch_f)
write.csv(all_lm_results_Lynch_LifeLines, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Alpha_diversity/LynchversusLifeLines_Shannon_linearmodel.csv")

all_plot_list_alpha_Lynch_LifeLines <- c(plot_list_alpha_Lynch_d, plot_list_alpha_Lynch_e, plot_list_alpha_Lynch_f)
ggarrange(plotlist = all_plot_list_alpha_Lynch_LifeLines, ncol = 4, nrow = 3)

alpha_Lynch_LifeLines_casescontrols <- c(plot_list_alpha_Lynch_e, plot_list_alpha_Lynch_f)
ggarrange(plotlist = alpha_Lynch_LifeLines_casescontrols, ncol = 4, nrow = 2)

all_lm_richness_Lynch_LifeLines <- bind_rows(all_lm_richness_Lynch_d, all_lm_richness_Lynch_e, all_lm_richness_Lynch_f)
write.csv(all_lm_richness_Lynch_LifeLines, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Alpha_diversity/LynchversusLifeLines_Richness_linearmodel.csv")

all_plot_list_richness_Lynch_LifeLines <- c(plot_list_richness_Lynch_d, plot_list_richness_Lynch_e, plot_list_richness_Lynch_f)
ggarrange(plotlist = all_plot_list_richness_Lynch_LifeLines, ncol = 4, nrow = 3)

richness_Lynch_LifeLines_casescontrols <- c(plot_list_richness_Lynch_e, plot_list_richness_Lynch_f)
ggarrange(plotlist = richness_Lynch_LifeLines_casescontrols, ncol = 4, nrow = 2)
