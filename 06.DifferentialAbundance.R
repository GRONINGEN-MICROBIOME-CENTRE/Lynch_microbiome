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

totest_Lynch_c <- c("NoNeo_High_Lowrisk")

totest_LynchGP_d <- c("cohort", "NohistoryCRC_cohort", "NohistoryCRCendo_cohort")

totest_LynchGP_e <- c("LynchLL_AdvAdenomas", "LynchLL_AdvAdenomasCRC", 'LynchLL_NonAdvAdenomas', "LynchLL_Neoplasia")

totest_LynchGP_f <- c("LynchLL_Controls", "LynchLL_NohistoryCRC_cohort", "LynchLL_NohistoryCRCendo_cohort")

totest_Lynch_neoplasia <- c("AllNeoplasia_controls","AdvAdenomasCRC_controls","AdvAdenomas_nonAdvAdenomas","CRC_controls", "AdvAdenomasCRCAdvSerr_controls", "nonAdvAdenomas_controls")
totest_Lynch_genes <- c("NoNeo_High_Lowrisk","NoNeo_Mutation")
totest_Lynch_risk <- c("AllNeoplasia_controls","LR_Neoplasia_controls","NohistoryCRC_LR_Neoplasia_controls","HR_Neoplasia_controls","NohistoryCRC_HR_Neoplasia_controls")

covariates_Lynch <- c("Sex", 'reads', "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

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
DAA_model_taxa_Lynch <- function(metadata, ID, CLR_transformed_data, phenotype, covariates) {
  df <- metadata
  row.names(df) <- df[,ID]
  df <- merge(df, CLR_transformed_data, by = 'row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Species <- colnames(CLR_transformed_data)
  Overall_result_phenos <- tibble() 
  
  for (Bug in Species) {
    if (! Bug %in% colnames(df)) { next }
    Bug2 <- paste(c("`", Bug, "`"), collapse = "")
    for (pheno in phenotype) {
      pheno2 <- paste(c("`", pheno, "`"), collapse = "")
      df$Participant_ID[!is.na(df[colnames(df) == pheno])] -> To_keep
      df_pheno_dda <- filter(df, Participant_ID %in% To_keep)
      length_df_phenodaa <- nrow(df_pheno_dda)
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(df_pheno_dda[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      
      formula_parts <- c(Bug2, "~", pheno2, "+", paste(valid_covariates, collapse = " + "))
      Model <- as.formula(paste(formula_parts, collapse = " "))
      lm(Model, df_pheno_dda) -> resultmodel
      as.data.frame(summary(resultmodel)$coefficients)[2,1:4] -> Summ_simple
      formula_text <- paste(deparse(Model), collapse = "")
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Bug = Bug, Pheno = pheno, Samples = length_df_phenodaa, Formula = as.character(formula_text) ) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p <- as.data.frame(Overall_result_phenos)
  p$FDR <- p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}

perform_DAA_taxa_Lynch <- function(metadata, ID, CLR_transformed_data, totest, covariates) {
  results_list <- list()
  for (variable in totest) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    result <- DAA_model_taxa_Lynch(metadata, ID, CLR_transformed_data, variable, covariates)
    results_list[[variable]] <- result
  }
  return(results_list)
}

#FDR nu ingesteld op 0.05
make_EMM_plots <- function(variables_to_test, results_list, CLR_transformed_data, metadata, covariates) {
  plot_list <- list()  # Initialize an empty list to store plots
  
  for (variable in variables_to_test) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    # Get the significant species
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(`Pr(>|t|)` < 0.05)
    sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.05)
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.01)
    
    # Check if there are significant species
    if (nrow(sig_diff_ab_bacteria) == 0) {
      cat("No significant findings for", variable, "\n")
    } else {
      # Get the estimated marginal means
      sig_species <- unique(sig_diff_ab_bacteria$Bug)
      sig_species <- gsub(".*s__", "s__", sig_species)
      df_species <- CLR_transformed_data
      colnames(df_species) <- gsub(".*s__", "s__", colnames(df_species))
      df_species_sig <- df_species[, sig_species, drop = FALSE]
      df <- merge(metadata, df_species_sig, by.x = "Participant_ID", by.y = "row.names")
      filtered_df <- df[!is.na(df[[variable]]), ]
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(filtered_df[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      test_variables <- c(variable, valid_covariates)
      
      create_lm <- function(species) {
        lm_formula <- as.formula(paste(species, "~", paste(test_variables, collapse = " + ")))
        #lm_formula <- as.formula(paste(species, "~", variable, "+ Age + Sex + BMI + Smoking + Bowel_Resection"))
        return(lm(lm_formula, data = df))
      }
      
      # Create a list of linear models
      lm_list <- lapply(sig_species, create_lm)
      
      # Calculate estimated marginal means for each linear model in lm_list
      formula <- as.formula(paste("~", variable))
      emm_list <- lapply(lm_list, emmeans, formula)
      
      # Combine the results into a single data frame
      result_df <- bind_rows(
        lapply(seq_along(sig_species), function(i) {
          emm_df <- as.data.frame(emm_list[[i]])
          emm_df$Species <- sig_species[i]
          return(emm_df)
        }))
      
      result_df <- result_df %>%
        group_by(Species) %>%
        mutate(Higher_in = levels(factor(!!sym(variable)))[which.max(emmean)]) %>%
        ungroup
      result_df$Higher_in <- as.factor(result_df$Higher_in)
      result_df <- result_df %>% arrange(Higher_in)
      
      # Reorder Species based on the rearranged result_df
      result_df$Species <- factor(result_df$Species, levels = unique(result_df$Species))
      
      # Make plot
      g <- ggplot(result_df, aes(x = emmean, y = Species, color = !!sym(variable))) +
        geom_point(size = 3) +
        geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
        theme_classic() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.text = element_text(size = 8, face = "bold", colour = "black"),
          legend.title = element_blank()) +
        scale_color_manual(values = my_colours) +
        labs(x = "EMM (CLR transformed data)", y = "Species") +
        ggtitle(paste("Differential abundance:", variable))
      
      # Print the plot
      print(g)
      
      # Store the plot in the list
      plot_list[[variable]] <- g
    }
  }
  
  return(plot_list)
}
make_EMM_plots_strains <- function(variables_to_test, results_list, CLR_transformed_data, metadata, covariates) {
  plot_list <- list()  # Initialize an empty list to store plots
  
  for (variable in variables_to_test) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    # Get the significant species
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(`Pr(>|t|)` < 0.05)
    sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.05)
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.01)
    
    # Check if there are significant species
    if (nrow(sig_diff_ab_bacteria) == 0) {
      cat("No significant findings for", variable, "\n")
    } else {
      # Get the estimated marginal means
      sig_species <- unique(sig_diff_ab_bacteria$Bug)
      sig_species <- gsub(".*t__", "t__", sig_species)
      df_species <- CLR_transformed_data
      colnames(df_species) <- gsub(".*t__", "t__", colnames(df_species))
      df_species_sig <- df_species[, sig_species, drop = FALSE]
      df <- merge(metadata, df_species_sig, by.x = "Participant_ID", by.y = "row.names")
      filtered_df <- df[!is.na(df[[variable]]), ]
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(filtered_df[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      test_variables <- c(variable, valid_covariates)
      
      create_lm <- function(species) {
        lm_formula <- as.formula(paste(species, "~", paste(test_variables, collapse = " + ")))
        #lm_formula <- as.formula(paste(species, "~", variable, "+ Age + Sex + BMI + Smoking + Bowel_Resection"))
        return(lm(lm_formula, data = df))
      }
      
      # Create a list of linear models
      lm_list <- lapply(sig_species, create_lm)
      
      # Calculate estimated marginal means for each linear model in lm_list
      formula <- as.formula(paste("~", variable))
      emm_list <- lapply(lm_list, emmeans, formula)
      
      # Combine the results into a single data frame
      result_df <- bind_rows(
        lapply(seq_along(sig_species), function(i) {
          emm_df <- as.data.frame(emm_list[[i]])
          emm_df$Species <- sig_species[i]
          return(emm_df)
        }))
      
      result_df <- result_df %>%
        group_by(Species) %>%
        mutate(Higher_in = levels(factor(!!sym(variable)))[which.max(emmean)]) %>%
        ungroup
      result_df$Higher_in <- as.factor(result_df$Higher_in)
      result_df <- result_df %>% arrange(Higher_in)
      
      # Reorder Species based on the rearranged result_df
      result_df$Species <- factor(result_df$Species, levels = unique(result_df$Species))
      
      # Make plot
      g <- ggplot(result_df, aes(x = emmean, y = Species, color = !!sym(variable))) +
        geom_point(size = 3) +
        geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
        theme_classic() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.text = element_text(size = 8, face = "bold", colour = "black"),
          legend.title = element_blank()) +
        scale_color_manual(values = my_colours) +
        labs(x = "EMM (CLR transformed data)", y = "Strains") +
        ggtitle(paste("Differential abundance:", variable))
      
      # Print the plot
      print(g)
      
      # Store the plot in the list
      plot_list[[variable]] <- g
    }
  }
  
  return(plot_list)
}
make_EMM_plots_genera <- function(variables_to_test, results_list, CLR_transformed_data, metadata, covariates) {
  plot_list <- list()  # Initialize an empty list to store plots
  
  for (variable in variables_to_test) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    # Get the significant species
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(`Pr(>|t|)` < 0.05)
    sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.05)
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.01)
    
    # Check if there are significant species
    if (nrow(sig_diff_ab_bacteria) == 0) {
      cat("No significant findings for", variable, "\n")
    } else {
      # Get the estimated marginal means
      sig_species <- unique(sig_diff_ab_bacteria$Bug)
      sig_species <- gsub(".*g__", "g__", sig_species)
      df_species <- CLR_transformed_data
      colnames(df_species) <- gsub(".*g__", "g__", colnames(df_species))
      df_species_sig <- df_species[, sig_species, drop = FALSE]
      df <- merge(metadata, df_species_sig, by.x = "Participant_ID", by.y = "row.names")
      filtered_df <- df[!is.na(df[[variable]]), ]
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(filtered_df[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      test_variables <- c(variable, valid_covariates)
      
      create_lm <- function(species) {
        lm_formula <- as.formula(paste(species, "~", paste(test_variables, collapse = " + ")))
        #lm_formula <- as.formula(paste(species, "~", variable, "+ Age + Sex + BMI + Smoking + Bowel_Resection"))
        return(lm(lm_formula, data = df))
      }
      
      # Create a list of linear models
      lm_list <- lapply(sig_species, create_lm)
      
      # Calculate estimated marginal means for each linear model in lm_list
      formula <- as.formula(paste("~", variable))
      emm_list <- lapply(lm_list, emmeans, formula)
      
      # Combine the results into a single data frame
      result_df <- bind_rows(
        lapply(seq_along(sig_species), function(i) {
          emm_df <- as.data.frame(emm_list[[i]])
          emm_df$Species <- sig_species[i]
          return(emm_df)
        }))
      
      result_df <- result_df %>%
        group_by(Species) %>%
        mutate(Higher_in = levels(factor(!!sym(variable)))[which.max(emmean)]) %>%
        ungroup
      result_df$Higher_in <- as.factor(result_df$Higher_in)
      result_df <- result_df %>% arrange(Higher_in)
      
      # Reorder Species based on the rearranged result_df
      result_df$Species <- factor(result_df$Species, levels = unique(result_df$Species))
      
      # Make plot
      g <- ggplot(result_df, aes(x = emmean, y = Species, color = !!sym(variable))) +
        geom_point(size = 3) +
        geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
        theme_classic() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.text = element_text(size = 8, face = "bold", colour = "black"),
          legend.title = element_blank()) +
        scale_color_manual(values = my_colours) +
        labs(x = "EMM (CLR transformed data)", y = "Strains") +
        ggtitle(paste("Differential abundance:", variable))
      
      # Print the plot
      print(g)
      
      # Store the plot in the list
      plot_list[[variable]] <- g
    }
  }
  
  return(plot_list)
}
make_EMM_plots_phyla <- function(variables_to_test, results_list, CLR_transformed_data, metadata, covariates) {
  plot_list <- list()  # Initialize an empty list to store plots
  
  for (variable in variables_to_test) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    # Get the significant species
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(`Pr(>|t|)` < 0.05)
    sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.05)
    #sig_diff_ab_bacteria <- results_list[[variable]] %>% filter(FDR < 0.01)
    
    # Check if there are significant species
    if (nrow(sig_diff_ab_bacteria) == 0) {
      cat("No significant findings for", variable, "\n")
    } else {
      # Get the estimated marginal means
      sig_species <- unique(sig_diff_ab_bacteria$Bug)
      sig_species <- gsub(".*p__", "p__", sig_species)
      df_species <- CLR_transformed_data
      colnames(df_species) <- gsub(".*p__", "p__", colnames(df_species))
      df_species_sig <- df_species[, sig_species, drop = FALSE]
      df <- merge(metadata, df_species_sig, by.x = "Participant_ID", by.y = "row.names")
      filtered_df <- df[!is.na(df[[variable]]), ]
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(filtered_df[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      test_variables <- c(variable, valid_covariates)
      
      create_lm <- function(species) {
        lm_formula <- as.formula(paste(species, "~", paste(test_variables, collapse = " + ")))
        #lm_formula <- as.formula(paste(species, "~", variable, "+ Age + Sex + BMI + Smoking + Bowel_Resection"))
        return(lm(lm_formula, data = df))
      }
      
      # Create a list of linear models
      lm_list <- lapply(sig_species, create_lm)
      
      # Calculate estimated marginal means for each linear model in lm_list
      formula <- as.formula(paste("~", variable))
      emm_list <- lapply(lm_list, emmeans, formula)
      
      # Combine the results into a single data frame
      result_df <- bind_rows(
        lapply(seq_along(sig_species), function(i) {
          emm_df <- as.data.frame(emm_list[[i]])
          emm_df$Species <- sig_species[i]
          return(emm_df)
        }))
      
      result_df <- result_df %>%
        group_by(Species) %>%
        mutate(Higher_in = levels(factor(!!sym(variable)))[which.max(emmean)]) %>%
        ungroup
      result_df$Higher_in <- as.factor(result_df$Higher_in)
      result_df <- result_df %>% arrange(Higher_in)
      
      # Reorder Species based on the rearranged result_df
      result_df$Species <- factor(result_df$Species, levels = unique(result_df$Species))
      
      # Make plot
      g <- ggplot(result_df, aes(x = emmean, y = Species, color = !!sym(variable))) +
        geom_point(size = 3) +
        geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
        theme_classic() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.text = element_text(size = 8, face = "bold", colour = "black"),
          legend.title = element_blank()) +
        scale_color_manual(values = my_colours) +
        labs(x = "EMM (CLR transformed data)", y = "Strains") +
        ggtitle(paste("Differential abundance:", variable))
      
      # Print the plot
      print(g)
      
      # Store the plot in the list
      plot_list[[variable]] <- g
    }
  }
  
  return(plot_list)
}


DAA_model_pwy <- function(metadata, ID, CLR_transformed_data, phenotype, covariates) {
  df <- metadata
  row.names(df) <- df[,ID]
  df <- merge(df, CLR_transformed_data, by = 'row.names')
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  Species <- colnames(CLR_transformed_data)
  Overall_result_phenos <- tibble() 
  
  for (Bug in Species) {
    if (! Bug %in% colnames(df)) { next }
    Bug2 <- paste(c("`", Bug, "`"), collapse = "")
    for (pheno in phenotype) {
      pheno2 <- paste(c("`", pheno, "`"), collapse = "")
      df$Participant_ID[!is.na(df[colnames(df) == pheno])] -> To_keep
      df_pheno_dda <- filter(df, Participant_ID %in% To_keep)
      length_df_phenodaa <- nrow(df_pheno_dda)
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(df_pheno_dda[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      
      formula_parts <- c(Bug2, "~", pheno2, "+", paste(valid_covariates, collapse = " + "))
      Model <- as.formula(paste(formula_parts, collapse = " "))
      lm(Model, df_pheno_dda) -> resultmodel
      as.data.frame(summary(resultmodel)$coefficients)[2,1:4] -> Summ_simple
      formula_text <- paste(deparse(Model), collapse = "")
      Summ_simple %>% rownames_to_column("Feature") %>% as_tibble() %>% mutate(Pathway = Bug, Pheno = pheno, Samples = length_df_phenodaa, Formula = as.character(formula_text)) -> temp_output
      rbind(Overall_result_phenos, temp_output) -> Overall_result_phenos
    }
  }
  
  p <- as.data.frame(Overall_result_phenos)
  p$FDR <- p.adjust(p$`Pr(>|t|)`, method = "BH")
  
  return(p)
}
perform_DAA_pwy <- function(metadata, ID, CLR_transformed_data, totest, covariates) {
  results_list <- list()
  for (variable in totest) {
    print(paste(' >>> ANALYSING VARIABLE <', variable, '>    <<<'))
    result <- DAA_model_pwy(metadata, ID, CLR_transformed_data, variable, covariates)
    results_list[[variable]] <- result
  }
  return(results_list)
}
make_EMM_plots_pwy <- function(variables_to_test, results_list, CLR_transformed_data, metadata, covariates) {
  plot_list <- list()  # Initialize an empty list to store plots
  
  for (variable in variables_to_test) {
    # Get the significant species
    #sig_diff_ab_pathways <- results_list[[variable]] %>% filter(`Pr(>|t|)` < 0.05)
    sig_diff_ab_pathways <- results_list[[variable]] %>% filter(FDR < 0.05)
    #sig_diff_ab_pathways <- results_list[[variable]] %>% filter(FDR < 0.01)
    
    # Check if there are significant pathways
    if (nrow(sig_diff_ab_pathways) == 0) {
      cat("No significant findings for", variable, "\n")
    } else {
      # Get the estimated marginal means
      sig_pathways <- unique(sig_diff_ab_pathways$Pathway)
      df_pathways <- CLR_transformed_data
      df_pathways_sig <- df_pathways[, sig_pathways, drop = FALSE]
      df <- merge(metadata, df_pathways_sig, by.x = "Participant_ID", by.y = "row.names")
      filtered_df <- df[!is.na(df[[variable]]), ]
      
      valid_covariates <- sapply(covariates, function(cov) length(unique(filtered_df[[cov]])) > 1)
      valid_covariates <- covariates[valid_covariates]
      test_variables <- c(variable, valid_covariates)
      
      create_lm <- function(pathways) {
        clean_pwyname <- paste("`", pathways, "`", sep = "")
        lm_formula <- as.formula(paste(clean_pwyname, "~", paste(test_variables, collapse = " + ")))
        #lm_formula <- as.formula(paste(clean_pwyname, "~", variable, "+ Age + Sex"))
        return(lm(lm_formula, data = df))
      }
      
      # Create a list of linear models
      lm_list <- lapply(sig_pathways, create_lm)
      
      # Calculate estimated marginal means
      formula <- as.formula(paste("~", variable))
      emm_list <- lapply(lm_list, emmeans, formula)
      
      # Combine the results into a single data frame
      result_df <- bind_rows(
        lapply(seq_along(sig_pathways), function(i) {
          emm_df <- as.data.frame(emm_list[[i]])
          emm_df$Pathway <- sig_pathways[i]
          return(emm_df)
        }))
      
      result_df <- result_df %>%
        group_by(Pathway) %>%
        mutate(Higher_in = levels(factor(!!sym(variable)))[which.max(emmean)]) %>%
        ungroup
      result_df$Higher_in <- as.factor(result_df$Higher_in)
      result_df <- result_df %>% arrange(Higher_in)
      
      # Reorder Pathway based on the rearranged result_df
      result_df$Pathway <- factor(result_df$Pathway, levels = unique(result_df$Pathway))
      
      # Make plot
      g <- ggplot(result_df, aes(x = emmean, y = Pathway, color = !!sym(variable))) +
        geom_point(size = 3) +
        geom_errorbar(aes(xmin = emmean - SE, xmax = emmean + SE), width = 0.2) +
        theme_classic() +
        theme(
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          legend.direction = "horizontal",
          legend.box = "vertical",
          legend.position = "top",
          legend.text = element_text(size = 8, face = "bold", colour = "black"),
          legend.title = element_blank()) +
        scale_color_manual(values = my_colours) +
        labs(x = "EMM (CLR transformed data)", y = "Pathway") +
        ggtitle(paste("Differential abundance:", variable))
      
      # Print the plot
      print(g)
      
      # Store the plot in the list
      plot_list[[variable]] <- g
    }
  }
  
  return(plot_list)
}

#### a. Differential abundance species ----
Lynch_all_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_speciesfiltered.csv")
Lynch_all_species_filt <- as.data.frame(Lynch_all_species_filt)
rownames(Lynch_all_species_filt) <- Lynch_all_species_filt[, 1]
Lynch_all_species_filt <- Lynch_all_species_filt[, -1]
Lynch_all_species_filt_clr <- decostand(Lynch_all_species_filt, method = 'clr', pseudocount = min(Lynch_all_species_filt[Lynch_all_species_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_a <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                        ID = "Participant_ID", 
                                        CLR_transformed_data = Lynch_all_species_filt_clr, 
                                        totest = totest_Lynch_a,
                                        covariates = covariates_daa)

view(results_Lynch_a$HR_Neoplasia_controls)
combined_results_species_a <- bind_rows(results_Lynch_a)

plot_list_Lynch_species_a <- make_EMM_plots(totest_Lynch_a, results_Lynch_a, Lynch_all_species_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_species_a, ncol = 4, nrow = 3)

#### a. Differential abundance strains ----
Lynch_all_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_strainsfiltered.csv")
Lynch_all_strains_filt <- as.data.frame(Lynch_all_strains_filt)
rownames(Lynch_all_strains_filt) <- Lynch_all_strains_filt[, 1]
Lynch_all_strains_filt <- Lynch_all_strains_filt[, -1]
Lynch_all_strains_filt_clr <- decostand(Lynch_all_strains_filt, method = 'clr', pseudocount = min(Lynch_all_strains_filt[Lynch_all_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_strains_a <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                        ID = "Participant_ID", 
                                        CLR_transformed_data = Lynch_all_strains_filt_clr, 
                                        totest = totest_Lynch_a,
                                        covariates = covariates_daa)

view(results_Lynch_strains_a$AllNeoplasia_controls)
combined_results_strains_a <- bind_rows(results_Lynch_strains_a)

plot_list_Lynch_strains_a <- make_EMM_plots_strains(totest_Lynch_a, results_Lynch_strains_a, Lynch_all_strains_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_strains_a, ncol = 4, nrow = 3)

#### a. Differential abundance genera ----
Lynch_all_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_generafiltered.csv")
Lynch_all_genera_filt <- as.data.frame(Lynch_all_genera_filt)
rownames(Lynch_all_genera_filt) <- Lynch_all_genera_filt[, 1]
Lynch_all_genera_filt <- Lynch_all_genera_filt[, -1]
Lynch_all_genera_filt_clr <- decostand(Lynch_all_genera_filt, method = 'clr', pseudocount = min(Lynch_all_genera_filt[Lynch_all_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_genera_a <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                ID = "Participant_ID", 
                                                CLR_transformed_data = Lynch_all_genera_filt_clr, 
                                                totest = totest_Lynch_a,
                                                covariates = covariates_daa)

view(results_Lynch_genera_a$AllNeoplasia_controls)
combined_results_genera_a <- bind_rows(results_Lynch_genera_a)

plot_list_Lynch_genera_a <- make_EMM_plots_genera(totest_Lynch_a, results_Lynch_genera_a, Lynch_all_genera_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_genera_a, ncol = 4, nrow = 3)

#### a. Differential abundance pathways ----
Lynch_all_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_pathwaysfiltered.csv")
Lynch_all_pathways_filt <- as.data.frame(Lynch_all_pathways_filt)
rownames(Lynch_all_pathways_filt) <- Lynch_all_pathways_filt[, 1]
Lynch_all_pathways_filt <- Lynch_all_pathways_filt[, -1]
Lynch_all_pathways_filt_clr <- decostand(Lynch_all_pathways_filt, method = 'clr', pseudocount = min(Lynch_all_pathways_filt[Lynch_all_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_pathways_a <- perform_DAA_pwy(metadata = meta_Lynch_baseline, 
                                               ID = "Participant_ID", 
                                               CLR_transformed_data = Lynch_all_pathways_filt_clr, 
                                               totest = totest_Lynch_a,
                                               covariates = covariates_daa)

view(results_Lynch_pathways_a$AllNeoplasia_controls)
combined_results_pathways_a <- bind_rows(results_Lynch_pathways_a)

plot_list_Lynch_pathways_a <- make_EMM_plots_pwy(totest_Lynch_a, results_Lynch_pathways_a, Lynch_all_pathways_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_pathways_a, ncol = 4, nrow = 3)


#### a. Differential abundance phyla ----
Lynch_all_phyla_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/a_Lynch_all_phylafiltered.csv")
Lynch_all_phyla_filt <- as.data.frame(Lynch_all_phyla_filt)
rownames(Lynch_all_phyla_filt) <- Lynch_all_phyla_filt[, 1]
Lynch_all_phyla_filt <- Lynch_all_phyla_filt[, -1]
Lynch_all_phyla_filt_clr <- decostand(Lynch_all_phyla_filt, method = 'clr', pseudocount = min(Lynch_all_phyla_filt[Lynch_all_phyla_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_phyla_a <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = Lynch_all_phyla_filt_clr, 
                                                 totest = totest_Lynch_a,
                                                 covariates = covariates_daa)

view(results_Lynch_phyla_a$AllNeoplasia_controls)
combined_results_phyla_a <- bind_rows(results_Lynch_phyla_a)

plot_list_Lynch_phyla_a <- make_EMM_plots_phyla(totest_Lynch_a, results_Lynch_phyla_a, Lynch_all_phyla_filt_clr, meta_Lynch_baseline, covariates_daa)

#### b. Differential abundance species ----
Lynch_neoplasia_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_speciesfiltered.csv")
Lynch_neoplasia_species_filt <- as.data.frame(Lynch_neoplasia_species_filt)
rownames(Lynch_neoplasia_species_filt) <- Lynch_neoplasia_species_filt[, 1]
Lynch_neoplasia_species_filt <- Lynch_neoplasia_species_filt[, -1]
Lynch_neoplasia_species_filt_clr <- decostand(Lynch_neoplasia_species_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_species_filt[Lynch_neoplasia_species_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_b <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                          ID = "Participant_ID", 
                                          CLR_transformed_data = Lynch_neoplasia_species_filt_clr, 
                                          totest = totest_Lynch_b,
                                          covariates = covariates_daa)

view(results_Lynch_b$AllNeoplasia_controls)
combined_results_species_b <- bind_rows(results_Lynch_b)

plot_list_Lynch_species_b <- make_EMM_plots(totest_Lynch_b, results_Lynch_b, Lynch_neoplasia_species_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_species_b, ncol = 1, nrow = 1)

#### b. Differential abundance strains ----
Lynch_neoplasia_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_strainsfiltered.csv")
Lynch_neoplasia_strains_filt <- as.data.frame(Lynch_neoplasia_strains_filt)
rownames(Lynch_neoplasia_strains_filt) <- Lynch_neoplasia_strains_filt[, 1]
Lynch_neoplasia_strains_filt <- Lynch_neoplasia_strains_filt[, -1]
Lynch_neoplasia_strains_filt_clr <- decostand(Lynch_neoplasia_strains_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_strains_filt[Lynch_neoplasia_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_strains_b <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                  ID = "Participant_ID", 
                                                  CLR_transformed_data = Lynch_neoplasia_strains_filt_clr, 
                                                  totest = totest_Lynch_b,
                                                  covariates = covariates_daa)

combined_results_strains_b <- bind_rows(results_Lynch_strains_b)

plot_list_Lynch_strains_b <- make_EMM_plots_strains(totest_Lynch_b, results_Lynch_strains_b, Lynch_neoplasia_strains_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_strains_b, ncol = 1, nrow = 1)

#### b. Differential abundance genera ----
Lynch_neoplasia_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_generafiltered.csv")
Lynch_neoplasia_genera_filt <- as.data.frame(Lynch_neoplasia_genera_filt)
rownames(Lynch_neoplasia_genera_filt) <- Lynch_neoplasia_genera_filt[, 1]
Lynch_neoplasia_genera_filt <- Lynch_neoplasia_genera_filt[, -1]
Lynch_neoplasia_genera_filt_clr <- decostand(Lynch_neoplasia_genera_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_genera_filt[Lynch_neoplasia_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_genera_b <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = Lynch_neoplasia_genera_filt_clr, 
                                                 totest = totest_Lynch_b,
                                                 covariates = covariates_daa)

combined_results_genera_b <- bind_rows(results_Lynch_genera_b)

plot_list_Lynch_genera_b <- make_EMM_plots_genera(totest_Lynch_b, results_Lynch_genera_b, Lynch_neoplasia_genera_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_genera_b, ncol = 1, nrow = 1)

#### b. Differential abundance pathways ----
Lynch_neoplasia_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/b_Lynch_neoplasia_pathwaysfiltered.csv")
Lynch_neoplasia_pathways_filt <- as.data.frame(Lynch_neoplasia_pathways_filt)
rownames(Lynch_neoplasia_pathways_filt) <- Lynch_neoplasia_pathways_filt[, 1]
Lynch_neoplasia_pathways_filt <- Lynch_neoplasia_pathways_filt[, -1]
Lynch_neoplasia_pathways_filt_clr <- decostand(Lynch_neoplasia_pathways_filt, method = 'clr', pseudocount = min(Lynch_neoplasia_pathways_filt[Lynch_neoplasia_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_pathways_b <- perform_DAA_pwy(metadata = meta_Lynch_baseline, 
                                            ID = "Participant_ID", 
                                            CLR_transformed_data = Lynch_neoplasia_pathways_filt_clr, 
                                            totest = totest_Lynch_b,
                                            covariates = covariates_daa)

combined_results_pathways_b <- bind_rows(results_Lynch_pathways_b)

plot_list_Lynch_pathways_b <- make_EMM_plots_pwy(totest_Lynch_b, results_Lynch_pathways_b, Lynch_neoplasia_pathways_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_pathways_b, ncol = 1, nrow = 1)


#### c. Differential abundance species ----
Lynch_control_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_speciesfiltered.csv")
Lynch_control_species_filt <- as.data.frame(Lynch_control_species_filt)
rownames(Lynch_control_species_filt) <- Lynch_control_species_filt[, 1]
Lynch_control_species_filt <- Lynch_control_species_filt[, -1]
Lynch_control_species_filt_clr <- decostand(Lynch_control_species_filt, method = 'clr', pseudocount = min(Lynch_control_species_filt[Lynch_control_species_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_c <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                          ID = "Participant_ID", 
                                          CLR_transformed_data = Lynch_control_species_filt_clr, 
                                          totest = totest_Lynch_c,
                                          covariates = covariates_daa)

view(results_Lynch_c$NoNeo_Mutation)
combined_results_species_c <- bind_rows(results_Lynch_c)

plot_list_Lynch_species_c <- make_EMM_plots(totest_Lynch_c, results_Lynch_c, Lynch_control_species_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_species_c, ncol = 2, nrow = 1)

# Also for NoNeo_Mutation (since this variable has four levels)
df_lm1 <- meta_Lynch_baseline[c("Participant_ID", "NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")]
df_lm1[sapply(df_lm1, is.character)] <- lapply(df_lm1[sapply(df_lm1, is.character)], as.factor)
df_lm1$Participant_ID <- as.character(df_lm1$Participant_ID)
levels(df_lm1$NoNeo_Mutation)

Lynch_control_species_filt_clr_a <- Lynch_control_species_filt_clr
Lynch_control_species_filt_clr_a$Participant_ID <- rownames(Lynch_control_species_filt_clr_a)
df_lm1 <- df_lm1 %>% filter(Participant_ID %in% Lynch_control_species_filt_clr_a$Participant_ID)
df_lm1 <- left_join(df_lm1, Lynch_control_species_filt_clr_a, by="Participant_ID")
Lynch_control_species_filt_clr_a$Participant_ID <- NULL

myContr1 <-list(MLH1_MSH2_EPCAM=c(-1,1,0,0),
                MLH1_MSH6=c(-1,0,1,0),
                MLH1_PMS2=c(-1,0,0,1),
                MSH2_EPCAM_MSH6=c(0,-1,1,0),
                MSH2_EPCAM_PMS2=c(0,-1,0,1),
                MSH6_PMS2=c(0,0,-1,1))

lm_contrasts <- vector("list", ncol(Lynch_control_species_filt_clr_a))
names(lm_contrasts) <- colnames(Lynch_control_species_filt_clr_a)
lm_models <- vector("list", ncol(Lynch_control_species_filt_clr_a))
names(lm_models) <- colnames(Lynch_control_species_filt_clr_a)

for(i in colnames(df_lm1)[grepl(x=colnames(df_lm1), pattern="k__")]) {
  #tryCatch({
  print(i)
  m_species <- lm(df_lm1[,i] ~ NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score, 
                  data=df_lm1[,c(i,"NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")])
  lm_models[[i]] <- m_species
  mm_species <- emmeans::emmeans(m_species,"NoNeo_Mutation")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_species, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="taxon_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",6)
names(lm_contrasts_sigFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")
lm_contrasts_allFDR <- vector("list",6)
names(lm_contrasts_allFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")

for(c in c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
  #mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
  # sign=ifelse(estimate>0, "positive", "negative"),
  #contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    mutate(#taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
      #sign=ifelse(estimate>0, "positive", "negative"),
      contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR %>% bind_rows() #all results

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR_contrasts %>%
  rename(
    Estimate = estimate,
    `Std. Error` = SE,
    `t value` = t.ratio,
    `Pr(>|t|)` = p.value,
    Feature = contrast,
    Bug = taxon_id,
    FDR = p_adj
  ) %>%
  mutate(
    Samples = NA,
    Pheno = "NoNeo_Mutation",
    Formula = "Bug + NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score")

combined_results_species_c <- full_join(combined_results_species_c, lm_contrasts_allFDR_contrasts)
combined_results_species_c$df <- NULL

#### c. Differential abundance strains ----
Lynch_control_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_strainsfiltered.csv")
Lynch_control_strains_filt <- as.data.frame(Lynch_control_strains_filt)
rownames(Lynch_control_strains_filt) <- Lynch_control_strains_filt[, 1]
Lynch_control_strains_filt <- Lynch_control_strains_filt[, -1]
Lynch_control_strains_filt_clr <- decostand(Lynch_control_strains_filt, method = 'clr', pseudocount = min(Lynch_control_strains_filt[Lynch_control_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_strains_c <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                  ID = "Participant_ID", 
                                                  CLR_transformed_data = Lynch_control_strains_filt_clr, 
                                                  totest = totest_Lynch_c,
                                                  covariates = covariates_daa)

combined_results_strains_c <- bind_rows(results_Lynch_strains_c)

plot_list_Lynch_strains_c <- make_EMM_plots_strains(totest_Lynch_c, results_Lynch_strains_c, Lynch_control_strains_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_strains_c, ncol = 2, nrow = 1)

# Also for NoNeo_Mutation (since this variable has four levels)
df_lm1 <- meta_Lynch_baseline[c("Participant_ID", "NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")]
df_lm1[sapply(df_lm1, is.character)] <- lapply(df_lm1[sapply(df_lm1, is.character)], as.factor)
df_lm1$Participant_ID <- as.character(df_lm1$Participant_ID)
levels(df_lm1$NoNeo_Mutation)

Lynch_control_strains_filt_clr_a <- Lynch_control_strains_filt_clr
Lynch_control_strains_filt_clr_a$Participant_ID <- rownames(Lynch_control_strains_filt_clr_a)
df_lm1 <- df_lm1 %>% filter(Participant_ID %in% Lynch_control_strains_filt_clr_a$Participant_ID)
df_lm1 <- left_join(df_lm1, Lynch_control_strains_filt_clr_a, by="Participant_ID")
Lynch_control_strains_filt_clr_a$Participant_ID <- NULL

myContr1 <-list(MLH1_MSH2_EPCAM=c(-1,1,0,0),
                MLH1_MSH6=c(-1,0,1,0),
                MLH1_PMS2=c(-1,0,0,1),
                MSH2_EPCAM_MSH6=c(0,-1,1,0),
                MSH2_EPCAM_PMS2=c(0,-1,0,1),
                MSH6_PMS2=c(0,0,-1,1))

lm_contrasts <- vector("list", ncol(Lynch_control_strains_filt_clr_a))
names(lm_contrasts) <- colnames(Lynch_control_strains_filt_clr_a)
lm_models <- vector("list", ncol(Lynch_control_strains_filt_clr_a))
names(lm_models) <- colnames(Lynch_control_strains_filt_clr_a)

for(i in colnames(df_lm1)[grepl(x=colnames(df_lm1), pattern="k__")]) {
  #tryCatch({
  print(i)
  m_species <- lm(df_lm1[,i] ~ NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score, 
                  data=df_lm1[,c(i,"NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")])
  lm_models[[i]] <- m_species
  mm_species <- emmeans::emmeans(m_species,"NoNeo_Mutation")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_species, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="taxon_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",6)
names(lm_contrasts_sigFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")
lm_contrasts_allFDR <- vector("list",6)
names(lm_contrasts_allFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")

for(c in c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
  #mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
  # sign=ifelse(estimate>0, "positive", "negative"),
  #contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    mutate(#taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
      #sign=ifelse(estimate>0, "positive", "negative"),
      contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR %>% bind_rows() #all results

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR_contrasts %>%
  rename(
    Estimate = estimate,
    `Std. Error` = SE,
    `t value` = t.ratio,
    `Pr(>|t|)` = p.value,
    Feature = contrast,
    Bug = taxon_id,
    FDR = p_adj
  ) %>%
  mutate(
    Samples = NA,
    Pheno = "NoNeo_Mutation",
    Formula = "Bug + NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score")

combined_results_strains_c <- full_join(combined_results_strains_c, lm_contrasts_allFDR_contrasts)
combined_results_strains_c$df <- NULL

#### c. Differential abundance genera ----
Lynch_control_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_generafiltered.csv")
Lynch_control_genera_filt <- as.data.frame(Lynch_control_genera_filt)
rownames(Lynch_control_genera_filt) <- Lynch_control_genera_filt[, 1]
Lynch_control_genera_filt <- Lynch_control_genera_filt[, -1]
Lynch_control_genera_filt_clr <- decostand(Lynch_control_genera_filt, method = 'clr', pseudocount = min(Lynch_control_genera_filt[Lynch_control_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_genera_c <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = Lynch_control_genera_filt_clr, 
                                                 totest = totest_Lynch_c,
                                                 covariates = covariates_daa)

combined_results_genera_c <- bind_rows(results_Lynch_genera_c)

plot_list_Lynch_genera_c <- make_EMM_plots_genera(totest_Lynch_c, results_Lynch_genera_c, Lynch_control_genera_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_genera_c, ncol = 2, nrow = 1)

# Also for NoNeo_Mutation (since this variable has four levels)
df_lm1 <- meta_Lynch_baseline[c("Participant_ID", "NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")]
df_lm1[sapply(df_lm1, is.character)] <- lapply(df_lm1[sapply(df_lm1, is.character)], as.factor)
df_lm1$Participant_ID <- as.character(df_lm1$Participant_ID)
levels(df_lm1$NoNeo_Mutation)

Lynch_control_genera_filt_clr_a <- Lynch_control_genera_filt_clr
Lynch_control_genera_filt_clr_a$Participant_ID <- rownames(Lynch_control_genera_filt_clr_a)
df_lm1 <- df_lm1 %>% filter(Participant_ID %in% Lynch_control_genera_filt_clr_a$Participant_ID)
df_lm1 <- left_join(df_lm1, Lynch_control_genera_filt_clr_a, by="Participant_ID")
Lynch_control_genera_filt_clr_a$Participant_ID <- NULL

myContr1 <-list(MLH1_MSH2_EPCAM=c(-1,1,0,0),
                MLH1_MSH6=c(-1,0,1,0),
                MLH1_PMS2=c(-1,0,0,1),
                MSH2_EPCAM_MSH6=c(0,-1,1,0),
                MSH2_EPCAM_PMS2=c(0,-1,0,1),
                MSH6_PMS2=c(0,0,-1,1))

lm_contrasts <- vector("list", ncol(Lynch_control_genera_filt_clr_a))
names(lm_contrasts) <- colnames(Lynch_control_genera_filt_clr_a)
lm_models <- vector("list", ncol(Lynch_control_genera_filt_clr_a))
names(lm_models) <- colnames(Lynch_control_genera_filt_clr_a)

for(i in colnames(df_lm1)[grepl(x=colnames(df_lm1), pattern="k__")]) {
  #tryCatch({
  print(i)
  m_species <- lm(df_lm1[,i] ~ NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score, 
                  data=df_lm1[,c(i,"NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")])
  lm_models[[i]] <- m_species
  mm_species <- emmeans::emmeans(m_species,"NoNeo_Mutation")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_species, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="taxon_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",6)
names(lm_contrasts_sigFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")
lm_contrasts_allFDR <- vector("list",6)
names(lm_contrasts_allFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")

for(c in c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
  #mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
  # sign=ifelse(estimate>0, "positive", "negative"),
  #contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    mutate(#taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
      #sign=ifelse(estimate>0, "positive", "negative"),
      contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR %>% bind_rows() #all results

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR_contrasts %>%
  rename(
    Estimate = estimate,
    `Std. Error` = SE,
    `t value` = t.ratio,
    `Pr(>|t|)` = p.value,
    Feature = contrast,
    Bug = taxon_id,
    FDR = p_adj
  ) %>%
  mutate(
    Samples = NA,
    Pheno = "NoNeo_Mutation",
    Formula = "Bug + NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score")

combined_results_genera_c <- full_join(combined_results_genera_c, lm_contrasts_allFDR_contrasts)
combined_results_genera_c$df <- NULL

#### c. Differential abundance pathways ----
Lynch_control_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_pathwaysfiltered.csv")
Lynch_control_pathways_filt <- as.data.frame(Lynch_control_pathways_filt)
rownames(Lynch_control_pathways_filt) <- Lynch_control_pathways_filt[, 1]
Lynch_control_pathways_filt <- Lynch_control_pathways_filt[, -1]
Lynch_control_pathways_filt_clr <- decostand(Lynch_control_pathways_filt, method = 'clr', pseudocount = min(Lynch_control_pathways_filt[Lynch_control_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_pathways_c <- perform_DAA_pwy(metadata = meta_Lynch_baseline, 
                                            ID = "Participant_ID", 
                                            CLR_transformed_data = Lynch_control_pathways_filt_clr, 
                                            totest = totest_Lynch_c,
                                            covariates = covariates_daa)

combined_results_pathways_c <- bind_rows(results_Lynch_pathways_c)

plot_list_Lynch_pathways_c <- make_EMM_plots_pwy(totest_Lynch_c, results_Lynch_pathways_c, Lynch_control_pathways_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_pathways_c, ncol = 2, nrow = 1)

# Also for NoNeo_Mutation (since this variable has four levels)
df_lm1 <- meta_Lynch_baseline[c("Participant_ID", "NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")]
df_lm1[sapply(df_lm1, is.character)] <- lapply(df_lm1[sapply(df_lm1, is.character)], as.factor)
df_lm1$Participant_ID <- as.character(df_lm1$Participant_ID)
levels(df_lm1$NoNeo_Mutation)

Lynch_control_pathways_filt_clr_a <- Lynch_control_pathways_filt_clr
Lynch_control_pathways_filt_clr_a$Participant_ID <- rownames(Lynch_control_pathways_filt_clr_a)
df_lm1 <- df_lm1 %>% filter(Participant_ID %in% Lynch_control_pathways_filt_clr_a$Participant_ID)
df_lm1 <- left_join(df_lm1, Lynch_control_pathways_filt_clr_a, by="Participant_ID")
Lynch_control_pathways_filt_clr_a$Participant_ID <- NULL

myContr1 <-list(MLH1_MSH2_EPCAM=c(-1,1,0,0),
                MLH1_MSH6=c(-1,0,1,0),
                MLH1_PMS2=c(-1,0,0,1),
                MSH2_EPCAM_MSH6=c(0,-1,1,0),
                MSH2_EPCAM_PMS2=c(0,-1,0,1),
                MSH6_PMS2=c(0,0,-1,1))

lm_contrasts <- vector("list", ncol(Lynch_control_pathways_filt_clr_a))
names(lm_contrasts) <- colnames(Lynch_control_pathways_filt_clr_a)
lm_models <- vector("list", ncol(Lynch_control_pathways_filt_clr_a))
names(lm_models) <- colnames(Lynch_control_pathways_filt_clr_a)

for (i in colnames(df_lm1)[grepl("PWY", colnames(df_lm1))]) {
  #tryCatch({
  print(i)
  m_species <- lm(df_lm1[,i] ~ NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score, 
                  data=df_lm1[,c(i,"NoNeo_Mutation", "Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")])
  lm_models[[i]] <- m_species
  mm_species <- emmeans::emmeans(m_species,"NoNeo_Mutation")
  lm_contrasts[[i]] <- as.data.frame(emmeans::contrast(mm_species, myContr1)) } #transform the emmGrid object into a data frame, which will be recognized as a vector in tibble

#Get all results in one dataframe
lm_raw_results <- lm_contrasts %>% bind_rows(.id="taxon_id") 

#Calculate FDR per contrast
lm_contrasts_sigFDR <- vector("list",6)
names(lm_contrasts_sigFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")
lm_contrasts_allFDR <- vector("list",6)
names(lm_contrasts_allFDR) <- c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")

for(c in c("MLH1_MSH2_EPCAM", "MLH1_MSH6","MLH1_PMS2","MSH2_EPCAM_MSH6","MSH2_EPCAM_PMS2","MSH6_PMS2")) {
  df <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    filter(p_adj < 0.05) #%>%
  #mutate(taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
  # sign=ifelse(estimate>0, "positive", "negative"),
  #contrast=factor(contrast))
  lm_contrasts_sigFDR[[c]] <- df
  
  df1 <- lm_raw_results %>%
    filter(contrast %in% c) %>%
    mutate(p_adj=p.adjust(p.value, "BH")) %>%
    mutate(#taxon_id=paste0(str_split_fixed(taxon_id,"\\.",7)[,7]),
      #sign=ifelse(estimate>0, "positive", "negative"),
      contrast=factor(contrast))
  lm_contrasts_allFDR[[c]] <- df1
}

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR %>% bind_rows() #all results

lm_contrasts_allFDR_contrasts <- lm_contrasts_allFDR_contrasts %>%
  rename(
    Estimate = estimate,
    `Std. Error` = SE,
    `t value` = t.ratio,
    `Pr(>|t|)` = p.value,
    Feature = contrast,
    Pathway = taxon_id,
    FDR = p_adj
  ) %>%
  mutate(
    Samples = NA,
    Pheno = "NoNeo_Mutation",
    Formula = "Bug + NoNeo_Mutation + Sex + Age + BMI + Smoking + Bowel_Resection + Bristol_score")

combined_results_pathways_c <- full_join(combined_results_pathways_c, lm_contrasts_allFDR_contrasts)
combined_results_pathways_c$df <- NULL

#### c. Differential abundance phyla ----
Lynch_control_phyla_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/c_Lynch_control_phylafiltered.csv")
Lynch_control_phyla_filt <- as.data.frame(Lynch_control_phyla_filt)
rownames(Lynch_control_phyla_filt) <- Lynch_control_phyla_filt[, 1]
Lynch_control_phyla_filt <- Lynch_control_phyla_filt[, -1]
Lynch_control_phyla_filt_clr <- decostand(Lynch_control_phyla_filt, method = 'clr', pseudocount = min(Lynch_control_phyla_filt[Lynch_control_phyla_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bowel_Resection", "Bristol_score")

results_Lynch_phyla_c <- perform_DAA_taxa_Lynch(metadata = meta_Lynch_baseline, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = Lynch_control_phyla_filt_clr, 
                                                 totest = totest_Lynch_c,
                                                 covariates = covariates_daa)

combined_results_phyla_c <- bind_rows(results_Lynch_phyla_c)

plot_list_Lynch_phyla_c <- make_EMM_plots_phyla(totest_Lynch_c, results_Lynch_phyla_c, Lynch_control_phyla_filt_clr, meta_Lynch_baseline, covariates_daa)
#ggarrange(plotlist = plot_list_Lynch_genera_c, ncol = 2, nrow = 1)
#### merging dataframes from all Lynch ----
# species
DAA_species_Lynch <- c(plot_list_Lynch_species_a, plot_list_Lynch_species_b, plot_list_Lynch_species_c)
#ggarrange(plotlist = DAA_species_Lynch, ncol = 3, nrow = 4)

DAA_species_Lynch_results <- bind_rows(combined_results_species_a, combined_results_species_b, combined_results_species_c)
write.csv(DAA_species_Lynch_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/Lynch_species.csv")

# heatmap nominal species
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_species_Lynch_results %>% filter(`Pr(>|t|)` < 0.05)
significant_results$Bug <- gsub(".*s__", "s__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Feature, Estimate) %>%
  spread(key = Feature, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Feature, Estimate)

ggplot(heatmap_data_long, aes(x = Feature, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Nominally Significant Species",
       x = "Pheno",
       y = "Bug")

# strains
DAA_strains_Lynch <- c(plot_list_Lynch_strains_a, plot_list_Lynch_strains_b, plot_list_Lynch_strains_c)
#ggarrange(plotlist = DAA_strains_Lynch, ncol = 3, nrow = 4)

DAA_strains_Lynch_results <- bind_rows(combined_results_strains_a, combined_results_strains_b, combined_results_strains_c)
write.csv(DAA_strains_Lynch_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/Lynch_strains.csv")

# heatmap nominal strains
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_strains_Lynch_results %>% filter(`Pr(>|t|)` < 0.05)
significant_results$Bug <- gsub(".*t__", "t__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Feature, Estimate) %>%
  spread(key = Feature, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Feature, Estimate)

ggplot(heatmap_data_long, aes(x = Feature, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Nominally Significant Strains",
       x = "Pheno",
       y = "Bug")

# genera
DAA_genera_Lynch <- c(plot_list_Lynch_genera_a, plot_list_Lynch_genera_b, plot_list_Lynch_genera_c)
#ggarrange(plotlist = DAA_genera_Lynch, ncol = 3, nrow = 4)

DAA_genera_Lynch_results <- bind_rows(combined_results_genera_a, combined_results_genera_b, combined_results_genera_c)
write.csv(DAA_genera_Lynch_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/Lynch_genera.csv")

# heatmap nominal genera
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_genera_Lynch_results %>% filter(`Pr(>|t|)` < 0.05)
significant_results$Bug <- gsub(".*g__", "g__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Feature, Estimate) %>%
  spread(key = Feature, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Feature, Estimate)

ggplot(heatmap_data_long, aes(x = Feature, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Nominally Significant Genera",
       x = "Pheno",
       y = "Bug")

# pathways
DAA_pathways_Lynch <- c(plot_list_Lynch_pathways_a, plot_list_Lynch_pathways_b, plot_list_Lynch_pathways_c)
#ggarrange(plotlist = DAA_pathways_Lynch, ncol = 3, nrow = 4)

DAA_pathways_Lynch_results <- bind_rows(combined_results_pathways_a, combined_results_pathways_b, combined_results_pathways_c)
write.csv(DAA_pathways_Lynch_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/Lynch_pathways.csv")

# heatmap nominal pathways
# Step 1 pathways: Filter the dataframe to include only nominally significant results
significant_results <- DAA_pathways_Lynch_results %>% filter(`Pr(>|t|)` < 0.05)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Pathway) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Pathway")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Pathway)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Pathway, Feature, Estimate) %>%
  spread(key = Feature, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Pathway, Feature, Estimate)

ggplot(heatmap_data_long, aes(x = Feature, y = factor(Pathway, levels = unique(Pathway)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of Nominally Significant Pathways",
       x = "Pheno",
       y = "Pathway")




#### d. Differential abundance species ----
LynchGP_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_speciesfiltered.csv")
LynchGP_species_filt <- as.data.frame(LynchGP_species_filt)
rownames(LynchGP_species_filt) <- LynchGP_species_filt[, 1]
LynchGP_species_filt <- LynchGP_species_filt[, -1]
LynchGP_species_filt_clr <- decostand(LynchGP_species_filt, method = 'clr', pseudocount = min(LynchGP_species_filt[LynchGP_species_filt > 0])/2)

LynchGP_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_metadata.csv")
LynchGP_meta <- as.data.frame(LynchGP_meta)
rownames(LynchGP_meta) <- LynchGP_meta$Participant_ID
LynchGP_meta <- LynchGP_meta[, -1]

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_d <- perform_DAA_taxa_Lynch(metadata = LynchGP_meta, 
                                          ID = "Participant_ID", 
                                          CLR_transformed_data = LynchGP_species_filt_clr, 
                                          totest = totest_LynchGP_d,
                                          covariates = covariates_daa)
view(results_Lynch_d$cohort)
combined_results_species_d <- bind_rows(results_Lynch_d)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_d) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_d)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_d[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_species_d <- make_EMM_plots(totest_LynchGP_d, results_Lynch_d, LynchGP_species_filt_clr, LynchGP_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_species_d, ncol = 3, nrow = 1)

#### d. Differential abundance strains ----
LynchGP_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_strainsfiltered.csv")
LynchGP_strains_filt <- as.data.frame(LynchGP_strains_filt)
rownames(LynchGP_strains_filt) <- LynchGP_strains_filt[, 1]
LynchGP_strains_filt <- LynchGP_strains_filt[, -1]
LynchGP_strains_filt_clr <- decostand(LynchGP_strains_filt, method = 'clr', pseudocount = min(LynchGP_strains_filt[LynchGP_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_strains_d <- perform_DAA_taxa_Lynch(metadata = LynchGP_meta, 
                                                  ID = "Participant_ID", 
                                                  CLR_transformed_data = LynchGP_strains_filt_clr, 
                                                  totest = totest_LynchGP_d,
                                                  covariates = covariates_daa)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_d) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_strains_d)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_strains_d[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

combined_results_strains_d <- bind_rows(results_Lynch_strains_d)
FDR_sig <- results_Lynch_strains_d$cohort %>% filter(`FDR` < 0.05)
FDR_sig <- results_Lynch_strains_d$NohistoryCRC_cohort %>% filter(`FDR` < 0.05)

plot_list_Lynch_strains_d <- make_EMM_plots_strains(totest_LynchGP_d, results_Lynch_strains_d, LynchGP_strains_filt_clr, LynchGP_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_strains_d, ncol = 2, nrow = 1)

#### d. Differential abundance genera ----
LynchGP_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_generafiltered.csv")
LynchGP_genera_filt <- as.data.frame(LynchGP_genera_filt)
rownames(LynchGP_genera_filt) <- LynchGP_genera_filt[, 1]
LynchGP_genera_filt <- LynchGP_genera_filt[, -1]
LynchGP_genera_filt_clr <- decostand(LynchGP_genera_filt, method = 'clr', pseudocount = min(LynchGP_genera_filt[LynchGP_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_genera_d <- perform_DAA_taxa_Lynch(metadata = LynchGP_meta, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = LynchGP_genera_filt_clr, 
                                                 totest = totest_LynchGP_d,
                                                 covariates = covariates_daa)
test <- results_Lynch_genera_d$cohort
library(openxlsx)
write.xlsx(test, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/genera_LynchversusLifeLines.xlsx")

combined_results_genera_d <- bind_rows(results_Lynch_genera_d)
FDR_sig <- results_Lynch_genera_d$cohort %>% filter(`FDR` < 0.05)
FDR_sig <- results_Lynch_genera_d$NohistoryCRC_cohort %>% filter(`FDR` < 0.05)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_d) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_genera_d)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_genera_d[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_genera_d <- make_EMM_plots_genera(totest_LynchGP_d, results_Lynch_genera_d, LynchGP_genera_filt_clr, LynchGP_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_genera_d, ncol = 2, nrow = 1)

#### d. Differential abundance phyla ----
LynchGP_phyla_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_phylafiltered.csv")
LynchGP_phyla_filt <- as.data.frame(LynchGP_phyla_filt)
rownames(LynchGP_phyla_filt) <- LynchGP_phyla_filt[, 1]
LynchGP_phyla_filt <- LynchGP_phyla_filt[, -1]
LynchGP_phyla_filt_clr <- decostand(LynchGP_phyla_filt, method = 'clr', pseudocount = min(LynchGP_phyla_filt[LynchGP_phyla_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_phyla_d <- perform_DAA_taxa_Lynch(metadata = LynchGP_meta, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = LynchGP_phyla_filt_clr, 
                                                 totest = totest_LynchGP_d,
                                                 covariates = covariates_daa)
test <- results_Lynch_phyla_d$cohort
library(openxlsx)
write.xlsx(test, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/phyla_LynchversusLifeLines.xlsx")

combined_results_phyla_d <- bind_rows(results_Lynch_phyla_d)
FDR_sig <- results_Lynch_phyla_d$cohort %>% filter(`FDR` < 0.05)
FDR_sig <- results_Lynch_phyla_d$NohistoryCRC_cohort %>% filter(`FDR` < 0.05)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_d) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_phyla_d)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_phyla_d[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_phyla_d <- make_EMM_plots_phyla(totest_LynchGP_d, results_Lynch_phyla_d, LynchGP_phyla_filt_clr, LynchGP_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_phyla_d, ncol = 2, nrow = 1)

#### d. Differential abundance pathways ----
LynchGP_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/d_LynchGeneralPopulation_pathwaysfiltered.csv")
LynchGP_pathways_filt <- as.data.frame(LynchGP_pathways_filt)
rownames(LynchGP_pathways_filt) <- LynchGP_pathways_filt[, 1]
LynchGP_pathways_filt <- LynchGP_pathways_filt[, -1]
LynchGP_pathways_filt_clr <- decostand(LynchGP_pathways_filt, method = 'clr', pseudocount = min(LynchGP_pathways_filt[LynchGP_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_pathways_d <- perform_DAA_pwy(metadata = LynchGP_meta, 
                                            ID = "Participant_ID", 
                                            CLR_transformed_data = LynchGP_pathways_filt_clr, 
                                            totest = totest_LynchGP_d,
                                            covariates = covariates_daa)

combined_results_pathways_d <- bind_rows(results_Lynch_pathways_d)
FDR_sig <- results_Lynch_pathways_d$cohort %>% filter(`FDR` < 0.05)
FDR_sig <- results_Lynch_pathways_d$NohistoryCRC_cohort %>% filter(`FDR` < 0.05)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_d) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_pathways_d)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_pathways_d[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_pathways_d <- make_EMM_plots_pwy(totest_LynchGP_d, results_Lynch_pathways_d, LynchGP_pathways_filt_clr, LynchGP_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_pathways_d, ncol = 2, nrow = 1)


#### e. Differential abundance species ----
LynchLL_neoplasia_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_speciesfiltered.csv")
LynchLL_neoplasia_species_filt <- as.data.frame(LynchLL_neoplasia_species_filt)
rownames(LynchLL_neoplasia_species_filt) <- LynchLL_neoplasia_species_filt[, 1]
LynchLL_neoplasia_species_filt <- LynchLL_neoplasia_species_filt[, -1]
LynchLL_neoplasia_species_filt_clr <- decostand(LynchLL_neoplasia_species_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_species_filt[LynchLL_neoplasia_species_filt > 0])/2)

LynchLL_neoplasia_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_metadata.csv")
LynchLL_neoplasia_meta <- as.data.frame(LynchLL_neoplasia_meta)
rownames(LynchLL_neoplasia_meta) <- LynchLL_neoplasia_meta$Participant_ID
LynchLL_neoplasia_meta <- LynchLL_neoplasia_meta[, -1]

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_e <- perform_DAA_taxa_Lynch(metadata = LynchLL_neoplasia_meta, 
                                          ID = "Participant_ID", 
                                          CLR_transformed_data = LynchLL_neoplasia_species_filt_clr, 
                                          totest = totest_LynchGP_e,
                                          covariates = covariates_daa)
view(results_Lynch_e$LynchLL_Neoplasia)
combined_results_species_e <- bind_rows(results_Lynch_e)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_e) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_e)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_e[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_species_e <- make_EMM_plots(totest_LynchGP_e, results_Lynch_e, LynchLL_neoplasia_species_filt_clr, LynchLL_neoplasia_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_species_e, ncol = 2, nrow = 2)

#### e. Differential abundance strains ----
LynchLL_neoplasia_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_strainsfiltered.csv")
LynchLL_neoplasia_strains_filt <- as.data.frame(LynchLL_neoplasia_strains_filt)
rownames(LynchLL_neoplasia_strains_filt) <- LynchLL_neoplasia_strains_filt[, 1]
LynchLL_neoplasia_strains_filt <- LynchLL_neoplasia_strains_filt[, -1]
LynchLL_neoplasia_strains_filt_clr <- decostand(LynchLL_neoplasia_strains_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_strains_filt[LynchLL_neoplasia_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_strains_e <- perform_DAA_taxa_Lynch(metadata = LynchLL_neoplasia_meta, 
                                                  ID = "Participant_ID", 
                                                  CLR_transformed_data = LynchLL_neoplasia_strains_filt_clr, 
                                                  totest = totest_LynchGP_e,
                                                  covariates = covariates_daa)

combined_results_strains_e <- bind_rows(results_Lynch_strains_e)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_e) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_strains_e)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_strains_e[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_strains_e <- make_EMM_plots_strains(totest_LynchGP_e, results_Lynch_strains_e, LynchLL_neoplasia_strains_filt_clr, LynchLL_neoplasia_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_strains_e, ncol = 2, nrow = 2)

#### e. Differential abundance genera ----
LynchLL_neoplasia_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_generafiltered.csv")
LynchLL_neoplasia_genera_filt <- as.data.frame(LynchLL_neoplasia_genera_filt)
rownames(LynchLL_neoplasia_genera_filt) <- LynchLL_neoplasia_genera_filt[, 1]
LynchLL_neoplasia_genera_filt <- LynchLL_neoplasia_genera_filt[, -1]
LynchLL_neoplasia_genera_filt_clr <- decostand(LynchLL_neoplasia_genera_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_genera_filt[LynchLL_neoplasia_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_genera_e <- perform_DAA_taxa_Lynch(metadata = LynchLL_neoplasia_meta, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = LynchLL_neoplasia_genera_filt_clr, 
                                                 totest = totest_LynchGP_e,
                                                 covariates = covariates_daa)

combined_results_genera_e <- bind_rows(results_Lynch_genera_e)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_e) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_genera_e)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_genera_e[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_genera_e <- make_EMM_plots_genera(totest_LynchGP_e, results_Lynch_genera_e, LynchLL_neoplasia_genera_filt_clr, LynchLL_neoplasia_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_genera_e, ncol = 2, nrow = 2)

#### e. Differential abundance pathways ----
LynchLL_neoplasia_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_pathwaysfiltered.csv")
LynchLL_neoplasia_pathways_filt <- as.data.frame(LynchLL_neoplasia_pathways_filt)
rownames(LynchLL_neoplasia_pathways_filt) <- LynchLL_neoplasia_pathways_filt[, 1]
LynchLL_neoplasia_pathways_filt <- LynchLL_neoplasia_pathways_filt[, -1]
LynchLL_neoplasia_pathways_filt_clr <- decostand(LynchLL_neoplasia_pathways_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_pathways_filt[LynchLL_neoplasia_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_pathways_e <- perform_DAA_pwy(metadata = LynchLL_neoplasia_meta, 
                                            ID = "Participant_ID", 
                                            CLR_transformed_data = LynchLL_neoplasia_pathways_filt_clr, 
                                            totest = totest_LynchGP_e,
                                            covariates = covariates_daa)

combined_results_pathways_e <- bind_rows(results_Lynch_pathways_e)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_e) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_pathways_e)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_pathways_e[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_pathways_e <- make_EMM_plots_pwy(totest_LynchGP_e, results_Lynch_pathways_e, LynchLL_neoplasia_pathways_filt_clr, LynchLL_neoplasia_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_pathways_e, ncol = 2, nrow = 2)


#### e. Differential abundance phyla ----
LynchLL_neoplasia_phyla_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/e_LynchLifeLines_cases_phylafiltered.csv")
LynchLL_neoplasia_phyla_filt <- as.data.frame(LynchLL_neoplasia_phyla_filt)
rownames(LynchLL_neoplasia_phyla_filt) <- LynchLL_neoplasia_phyla_filt[, 1]
LynchLL_neoplasia_phyla_filt <- LynchLL_neoplasia_phyla_filt[, -1]
LynchLL_neoplasia_phyla_filt_clr <- decostand(LynchLL_neoplasia_phyla_filt, method = 'clr', pseudocount = min(LynchLL_neoplasia_phyla_filt[LynchLL_neoplasia_phyla_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_phyla_e <- perform_DAA_taxa_Lynch(metadata = LynchLL_neoplasia_meta, 
                                                ID = "Participant_ID", 
                                                CLR_transformed_data = LynchLL_neoplasia_phyla_filt_clr, 
                                                totest = totest_LynchGP_e,
                                                covariates = covariates_daa)

combined_results_phyla_e <- bind_rows(results_Lynch_phyla_e)
FDR_sig <- results_Lynch_phyla_e$LynchLL_Neoplasia %>% filter(`FDR` < 0.05)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_e) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_phyla_e)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_phyla_e[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_phyla_e <- make_EMM_plots_phyla(totest_LynchGP_e, results_Lynch_phyla_e, LynchLL_neoplasia_phyla_filt_clr, LynchLL_neoplasia_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_phyla_e, ncol = 2, nrow = 2)

#### f. Differential abundance species ----
LynchLL_controls_species_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_speciesfiltered.csv")
LynchLL_controls_species_filt <- as.data.frame(LynchLL_controls_species_filt)
rownames(LynchLL_controls_species_filt) <- LynchLL_controls_species_filt[, 1]
LynchLL_controls_species_filt <- LynchLL_controls_species_filt[, -1]
LynchLL_controls_species_filt_clr <- decostand(LynchLL_controls_species_filt, method = 'clr', pseudocount = min(LynchLL_controls_species_filt[LynchLL_controls_species_filt > 0])/2)

LynchLL_controls_meta <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_metadata.csv")
LynchLL_controls_meta <- as.data.frame(LynchLL_controls_meta)
rownames(LynchLL_controls_meta) <- LynchLL_controls_meta$Participant_ID
LynchLL_controls_meta <- LynchLL_controls_meta[, -1]

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_f <- perform_DAA_taxa_Lynch(metadata = LynchLL_controls_meta, 
                                          ID = "Participant_ID", 
                                          CLR_transformed_data = LynchLL_controls_species_filt_clr, 
                                          totest = totest_LynchGP_f,
                                          covariates = covariates_daa)
view(results_Lynch_f$LynchLL_Controls)
combined_results_species_f <- bind_rows(results_Lynch_f)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_f) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_f)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_f[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_species_f <- make_EMM_plots(totest_LynchGP_f, results_Lynch_f, LynchLL_controls_species_filt_clr, LynchLL_controls_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_species_f, ncol = 3, nrow = 1)

#### f. Differential abundance strains ----
LynchLL_controls_strains_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_strainsfiltered.csv")
LynchLL_controls_strains_filt <- as.data.frame(LynchLL_controls_strains_filt)
rownames(LynchLL_controls_strains_filt) <- LynchLL_controls_strains_filt[, 1]
LynchLL_controls_strains_filt <- LynchLL_controls_strains_filt[, -1]
LynchLL_controls_strains_filt_clr <- decostand(LynchLL_controls_strains_filt, method = 'clr', pseudocount = min(LynchLL_controls_strains_filt[LynchLL_controls_strains_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_strains_f <- perform_DAA_taxa_Lynch(metadata = LynchLL_controls_meta, 
                                                  ID = "Participant_ID", 
                                                  CLR_transformed_data = LynchLL_controls_strains_filt_clr, 
                                                  totest = totest_LynchGP_f,
                                                  covariates = covariates_daa)

combined_results_strains_f <- bind_rows(results_Lynch_strains_f)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_f) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_strains_f)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_strains_f[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_strains_f <- make_EMM_plots_strains(totest_LynchGP_f, results_Lynch_strains_f, LynchLL_controls_strains_filt_clr, LynchLL_controls_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_strains_f, ncol = 1, nrow = 1)

#### f. Differential abundance genera ----
LynchLL_controls_genera_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_generafiltered.csv")
LynchLL_controls_genera_filt <- as.data.frame(LynchLL_controls_genera_filt)
rownames(LynchLL_controls_genera_filt) <- LynchLL_controls_genera_filt[, 1]
LynchLL_controls_genera_filt <- LynchLL_controls_genera_filt[, -1]
LynchLL_controls_genera_filt_clr <- decostand(LynchLL_controls_genera_filt, method = 'clr', pseudocount = min(LynchLL_controls_genera_filt[LynchLL_controls_genera_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_genera_f <- perform_DAA_taxa_Lynch(metadata = LynchLL_controls_meta, 
                                                 ID = "Participant_ID", 
                                                 CLR_transformed_data = LynchLL_controls_genera_filt_clr, 
                                                 totest = totest_LynchGP_f,
                                                 covariates = covariates_daa)

combined_results_genera_f <- bind_rows(results_Lynch_genera_f)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_f) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_genera_f)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_genera_f[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_genera_f <- make_EMM_plots_genera(totest_LynchGP_f, results_Lynch_genera_f, LynchLL_controls_genera_filt_clr, LynchLL_controls_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_genera_f, ncol = 1, nrow = 1)

#### f. Differential abundance phyla ----
LynchLL_controls_phyla_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_phylafiltered.csv")
LynchLL_controls_phyla_filt <- as.data.frame(LynchLL_controls_phyla_filt)
rownames(LynchLL_controls_phyla_filt) <- LynchLL_controls_phyla_filt[, 1]
LynchLL_controls_phyla_filt <- LynchLL_controls_phyla_filt[, -1]
LynchLL_controls_phyla_filt_clr <- decostand(LynchLL_controls_phyla_filt, method = 'clr', pseudocount = min(LynchLL_controls_phyla_filt[LynchLL_controls_phyla_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_phyla_f <- perform_DAA_taxa_Lynch(metadata = LynchLL_controls_meta, 
                                                ID = "Participant_ID", 
                                                CLR_transformed_data = LynchLL_controls_phyla_filt_clr, 
                                                totest = totest_LynchGP_f,
                                                covariates = covariates_daa)

combined_results_phyla_f <- bind_rows(results_Lynch_phyla_f)
FDR_sig <- results_Lynch_phyla_f$LynchLL_Controls %>% filter(`FDR` < 0.05)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_f) {
  # Check if the variable exists in results_Lynch_d
  if (var %in% names(results_Lynch_phyla_f)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_phyla_f[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
summary_df

plot_list_Lynch_phyla_f <- make_EMM_plots_phyla(totest_LynchGP_f, results_Lynch_phyla_f, LynchLL_controls_phyla_filt_clr, LynchLL_controls_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_phyla_f, ncol = 2, nrow = 2)


#### f. Differential abundance pathways ----
LynchLL_controls_pathways_filt <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Dataframes_comparisons/f_LynchLifeLines_controls_pathwaysfiltered.csv")
LynchLL_controls_pathways_filt <- as.data.frame(LynchLL_controls_pathways_filt)
rownames(LynchLL_controls_pathways_filt) <- LynchLL_controls_pathways_filt[, 1]
LynchLL_controls_pathways_filt <- LynchLL_controls_pathways_filt[, -1]
LynchLL_controls_pathways_filt_clr <- decostand(LynchLL_controls_pathways_filt, method = 'clr', pseudocount = min(LynchLL_controls_pathways_filt[LynchLL_controls_pathways_filt > 0])/2)

covariates_daa <- c("Sex", "Age", "BMI", "Smoking", "Bristol_score")

results_Lynch_pathways_f <- perform_DAA_pwy(metadata = LynchLL_controls_meta, 
                                            ID = "Participant_ID", 
                                            CLR_transformed_data = LynchLL_controls_pathways_filt_clr, 
                                            totest = totest_LynchGP_f,
                                            covariates = covariates_daa)

combined_results_pathways_f <- bind_rows(results_Lynch_pathways_f)

summary_df <- data.frame(Variable = character(), Significant_Count = integer(), stringsAsFactors = FALSE)
for (var in totest_LynchGP_f) {
  # Check if the variable exists in results_Lynch_e
  if (var %in% names(results_Lynch_pathways_f)) {
    # Count the number of significant findings
    significant_count <- results_Lynch_pathways_f[[var]] %>% filter(FDR < 0.05) %>% nrow()
    summary_df <- summary_df %>% add_row(Variable = var, Significant_Count = significant_count)
  }
}
print(summary_df)

plot_list_Lynch_pathways_f <- make_EMM_plots_pwy(totest_LynchGP_f, results_Lynch_pathways_f, LynchLL_controls_pathways_filt_clr, LynchLL_controls_meta, covariates_daa)
ggarrange(plotlist = plot_list_Lynch_pathways_f, ncol = 2, nrow = 2)


#### merging dataframes from all Lynch ----
# species
DAA_species_LynchLifelines <- c(plot_list_Lynch_species_d, plot_list_Lynch_species_e, plot_list_Lynch_species_f)
ggarrange(plotlist = DAA_species_LynchLifelines, ncol = 3, nrow = 3)

DAA_species_LynchLifeLines_results <- bind_rows(combined_results_species_d, combined_results_species_e, combined_results_species_f)
write.csv(DAA_species_LynchLifeLines_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/LynchLifeLines_species.csv")

# heatmap FDR species
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_species_LynchLifeLines_results %>% filter(`FDR` < 0.05)
significant_results$Bug <- gsub(".*s__", "s__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Pheno, Estimate) %>%
  spread(key = Pheno, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Pheno, Estimate)

ggplot(heatmap_data_long, aes(x = Pheno, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of FDR Significant Species",
       x = "Pheno",
       y = "Bug")

# strains
DAA_strains_LynchLifeLines <- c(plot_list_Lynch_strains_d, plot_list_Lynch_strains_e, plot_list_Lynch_strains_f)
ggarrange(plotlist = DAA_strains_LynchLifeLines, ncol = 3, nrow = 3)

DAA_strains_LynchLifeLines_results <- bind_rows(combined_results_strains_d, combined_results_strains_e, combined_results_strains_f)
write.csv(DAA_strains_LynchLifeLines_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/LynchLifeLines_strains.csv")

# heatmap FDR strains
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_strains_LynchLifeLines_results %>% filter(`FDR` < 0.05)
significant_results$Bug <- gsub(".*t__", "t__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Pheno, Estimate) %>%
  spread(key = Pheno, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Pheno, Estimate)

ggplot(heatmap_data_long, aes(x = Pheno, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of FDR Significant Strains",
       x = "Pheno",
       y = "Bug")

# genera
DAA_genera_LynchLifeLines<- c(plot_list_Lynch_genera_d, plot_list_Lynch_genera_e, plot_list_Lynch_genera_f)
ggarrange(plotlist = DAA_genera_LynchLifeLines, ncol = 3, nrow = 3)

DAA_genera_LynchLifeLines_results <- bind_rows(combined_results_genera_d, combined_results_genera_e, combined_results_genera_f)
write.csv(DAA_genera_LynchLifeLines_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/LynchLifeLines_genera.csv")

# heatmap FDR genera
# Step 1 species: Filter the dataframe to include only nominally significant results
significant_results <- DAA_genera_LynchLifeLines_results %>% filter(`FDR` < 0.05)
significant_results$Bug <- gsub(".*g__", "g__", significant_results$Bug)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Bug) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Bug")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Bug)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Bug, Pheno, Estimate) %>%
  spread(key = Pheno, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Bug, Pheno, Estimate)

ggplot(heatmap_data_long, aes(x = Pheno, y = factor(Bug, levels = unique(Bug)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of FDR Significant Genera",
       x = "Pheno",
       y = "Bug")

# pathways
DAA_pathways_LynchLifeLines <- c(plot_list_Lynch_pathways_d, plot_list_Lynch_pathways_e, plot_list_Lynch_pathways_f)
ggarrange(plotlist = DAA_pathways_LynchLifeLines, ncol = 3, nrow = 3)

DAA_pathways_LynchLifeLines_results <- NULL
DAA_pathways_LynchLifeLines_results <- bind_rows(combined_results_pathways_d, combined_results_pathways_e, combined_results_pathways_f)
write.csv(DAA_pathways_LynchLifeLines_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/LynchLifeLines_pathways.csv")

# heatmap FDR pathways
# Step 1 pathways: Filter the dataframe to include only nominally significant results
significant_results <- DAA_pathways_LynchLifeLines_results %>% filter(`FDR` < 0.05)

# Step 2: Calculate the count of significant phenos for each Bug
bug_counts <- significant_results %>%
  group_by(Pathway) %>%
  summarise(significant_pheno_count = n())

significant_results <- significant_results %>%
  left_join(bug_counts, by = "Pathway")

# Step 3: Order the Bugs by the count of significant phenos
significant_results <- significant_results %>%
  arrange(desc(significant_pheno_count), Pathway)

# Step 4: Reshape the dataframe for the heatmap
heatmap_data <- significant_results %>%
  select(Pathway, Pheno, Estimate) %>%
  spread(key = Pheno, value = Estimate)

# Step 5: Create the heatmap
heatmap_data_long <- significant_results %>%
  select(Pathway, Pheno, Estimate)

ggplot(heatmap_data_long, aes(x = Pheno, y = factor(Pathway, levels = unique(Pathway)), fill = Estimate)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "blue", high = "red", na.value = "grey50", name = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Heatmap of FDR Significant Pathways",
       x = "Pheno",
       y = "Pathway")

# phyla
DAA_phyla_LynchLifeLines<- c(plot_list_Lynch_phyla_d)
ggarrange(plotlist = DAA_phyla_LynchLifeLines, ncol = 3, nrow = 1)

DAA_phyla_LynchLifeLines_results <- bind_rows(combined_results_phyla_d)
write.csv(DAA_phyla_LynchLifeLines_results, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Results/Differential_abundance/LynchLifeLines_phyla.csv")




