# LYNCH SYNDROME PROJECT, COLLABORATION AMSTERDAM GRONINGEN
# by Femke Prins, Q2 2024
# script for matching LifeLines samples and PALGA data

library(MatchIt)
setwd("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Raw_data")

#### Script for selecting data from the cluster ----
variables = read.table("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v27.csv", sep=",", header=T) %>% as_tibble()

variables_relevant <- variables %>% select(c("DAG3_sampleID", "ANTHRO.AGE", "ANTHRO.Sex", "ANTHRO.BMI",
                                             "EXP.SMOKING.Smoker.Now", "EXP.SMOKING.Smoker.stopped", 
                                             "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD", "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.UC", 
                                             "MED.DISEASES.Endocrine.DiabetesT2", "MED.DISEASES.Endocrine.Autoimmune.DiabetesT1",
                                             "MED.DISEASES.Cardiovascular.Hypertension", "MED.MEDS.Laxatives_osmotic_ATC_A06AD", "MED.MEDS.Laxatives_volume_increasing_ATC_A06AC",
                                             "EXP.DIET.Probiotics", "MED.MEDS.Vitamine_B_ATC_A11DA", "MED.MEDS.Vitamine_D_ATC_A11CC",
                                             "MED.MEDS.Vitamine_K_antagonists_ATC_B01AA", "MED.MEDS.Vitamine_supplements_ATC_A11A",
                                             "META.Antibiotics_3m", "MED.DISEASES.Cancer.Colon", "SOCIOECONOMIC.BUURT.Urbanicity",
                                             "MED.MEDS.PPIs_ATC_A02BC", "META.POOP.BristolMean", "META.POOP.COLLECTION_SEASON", "META.BATCH"))

table(variables_relevant$MED.DISEASES.Cancer.Colon)
write_tsv(variables_relevant, "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-rgacesa/PALGA/metadata_for_Femke/variables_for_matching_Lynch.tsv")

#### Loading metagenomic data DAG3----
prepCleanHumann <- function(inPath, dropUnintegrated = F, dropUnmapped = F, dropTaxonSpecific = T, presenceFilter = -1, minRelativeAbundance = -1, novogeneIdsClean = T) {
  inDF <- read.table(inPath,sep='\t',header=T,quote ='',comment.char = '')
  
  #fix pathway ID (these tend to be weird coming out of humann)
  rownames(inDF) <- inDF$X..Pathway
  inDF$X..Pathway <- NULL
  
  #drop "junk" from humann (unintegrated/unmapped data & taxon-specific pathways)
  if (dropTaxonSpecific) {
    inDF <- inDF[grep('\\|',rownames(inDF),invert = T),]
  }
  if (dropUnintegrated) {
    inDF <- inDF[grep('UNINTEGRATED',rownames(inDF),invert = T),]
  }
  if (dropUnmapped) {
    inDF <- inDF[grep('UNMAPPED',rownames(inDF),invert = T),]
  }
  rownames(inDF)[grep('PWY',rownames(inDF),invert = T)] <- paste0('PWY_',rownames(inDF)[grep('PWY',rownames(inDF),invert = T)]) # adds PWY to the rownames if not already present
  inDF <- as.data.frame(t.data.frame(inDF))
  
  # fix sample IDs, remove duplicates
  inDF$ID <- rownames(inDF)
  if (novogeneIdsClean) {
    print ('NOTE: doing cleaning of Novogene IDs')
    inDF$ID <- sub("_EKDN.*", "", inDF$ID)
  }
  if (sum(duplicated(inDF$ID) > 0)) {
    print(paste('WARNING: found ',sum(duplicated(inDF$ID) > 0),'duplicates, dropping them!'))
  }
  
  duplicates <- inDF$ID[duplicated(inDF$ID)]
  inDF <- inDF %>% filter(!ID %in% duplicates)
  rownames(inDF) <- inDF$ID
  inDF$ID <- NULL
  
  # make sure columns are actually numbers
  for (c in colnames(inDF)) {inDF[[c]] <- as.numeric(inDF[[c]])}
  
  # clean, rescale and save #not in this step, later per cohort
  #inDFt2 <- filterHumannDF(inDF,presPerc = presenceFilter,minMRelAb = minRelativeAbundance,minMedRelAb = -1,rescale = T,minSum = 1,verbose = T)
  inDFt2 <- inDF
  inDFt2 <- inDFt2[,colSums(inDFt2)!=0]
  inDFt2$ID <- row.names(inDFt2)
  print(paste('Done, returning ',nrow(inDFt2),'samples'))
  inDFt2
} #Ranko's function, adjusted a bit to fit this cohort

#pathways
inPath <- "DAG3_humann3.6_pathway_abundances_metacyc_2024Apr.txt" 
DAG3_pathways <- prepCleanHumann(inPath) #run once because takes a long time to load
rownames(DAG3_pathways) <- sub("_kneaddata.*", "", rownames(DAG3_pathways))
DAG3_pathways$ID <- NULL #8352 samples
row_sums <- rowSums(DAG3_pathways)
DAG3_pathways_cleaned <- DAG3_pathways[row_sums != 0, ] #8172 samples

#taxonomic data
my_DAG3 <- read.delim("DAG3_merged_metaphlan4_2024Apr.txt", header = FALSE) #will take a long time to load
tax <- my_DAG3 
tax <- tax[-1, ]
tax %>% t() %>% as.data.frame() -> tax2
colnames(tax2) <- tax2[1, ]
tax2 <- tax2[-1, ]
rownames(tax2) <- tax2$clade_name
tax2 <- tax2[, !(colnames(tax2) == "clade_name")]
tax3 <- as.data.frame(apply(tax2, 2, as.numeric))
rownames(tax3) <- rownames(tax2)
rownames(tax3) <- sub("_metaphlan", "", rownames(tax3)) #8354 obs, 9421 var

#samples with all layers available
tax4 <- tax3[rownames(tax3) %in% rownames(DAG3_pathways_cleaned), ] #8077 with available taxa data AND pathways

#### Preparing metadata "general population" ----
LifeLines_meta <- read_delim("variables_for_matching_Lynch2.tsv") #8229
reads_DAG3 <- read.csv("DAG3_kneaddata_merged.csv", header = TRUE)
reads_DAG3$Sample <- gsub("^(([^_]*_){2}[^_]*)_.*$", "\\1", reads_DAG3$Sample)
reads_DAG3 <- subset(reads_DAG3, select = c("Sample", "Final..p2."))
colnames(reads_DAG3) <- c("DAG3_sampleID", "reads")
LifeLines_meta$Bristol_score <- LifeLines_meta$META.POOP.BristolMean
LifeLines_meta$Bristol_score[is.na(LifeLines_meta$Bristol_score)] <- 4

LifeLines_meta_all <- merge(LifeLines_meta, reads_DAG3, by = "DAG3_sampleID") #7997
LifeLines_meta_gp <- merge(LifeLines_meta, reads_DAG3, by = "DAG3_sampleID") #7997
LifeLines_meta_gp <- LifeLines_meta_gp[LifeLines_meta_gp$DAG3_sampleID %in% rownames(tax4), ] #7819 GP with available taxa data AND pathways

#Excluding samples before matching (history of CRC in PALGA, self-reported CRC, self-reported IBD, antibiotic use <3m)
LifeLines_meta_gp <- LifeLines_meta_gp %>%
  filter(MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD == "N" & MED.DISEASES.Gastrointestinal.Autoimmune.IBD.UC == "N") #7728

PALGA <- read_excel("Copy of METADATA_PALGA_data_v7.xlsx", na = "NA")
PALGA_CRC <- PALGA %>% filter(Polyp.Neoplasia.T == "Adv.Neoplasia.CRC") #50 with coloncancer
PALGA_CRC <- PALGA_CRC %>% filter(META.deltaDate.ba == "Before.sampling") #24 with coloncancer
LifeLines_meta_gp <- LifeLines_meta_gp[!LifeLines_meta_gp$DAG3_sampleID %in% PALGA_CRC$DAG3_ID, ] #7708 participants left

LifeLines_meta_gp <- LifeLines_meta_gp %>% filter(MED.DISEASES.Cancer.Colon == "N") #7705

LifeLines_meta_gp <- LifeLines_meta_gp %>% filter(META.Antibiotics_3m == "N") #5896

#### Preparing metadata "controls" ----
# Exclude participants with self-reported endoscopies and in PALGA database
colonoscopies_DAG3 <- read_csv("PALGA_colonoscopies_DAG3.csv")
colonoscopies_DAG3_no <- colonoscopies_DAG3 %>% filter(`2b_q_1.gastrointestinal_endoscopy_adu_q_1_a` == 2) #3083
LifeLines_meta_controls <- LifeLines_meta_gp[LifeLines_meta_gp$DAG3_sampleID %in% colonoscopies_DAG3_no$DAG3_sampleID, ] #2203 that reported NO endoscopies and are in the metadata

PALGA_endo <- PALGA %>% filter(META.deltaDate.ba == "After.sampling") #706 with hit after sampling
LifeLines_meta_controls <- LifeLines_meta_controls[!LifeLines_meta_controls$DAG3_sampleID %in% PALGA_endo$DAG3_ID, ] #2142 not in PALGA after sampling and NO reported endoscopies

#### Merging metadata from both cohorts for matching controls----
Lynch_x <- meta_Lynch_baseline %>% filter(AllNeoplasia_controls == "Control") %>%
  select(Participant_ID, Sex, BMI, Age, Smoking)
Lynch_x <- mutate(Lynch_x, cohort = 'Lynch')
str(Lynch_x)

LifeLines_meta_controls <- LifeLines_meta_controls %>%
  mutate(Smoking = case_when(
    EXP.SMOKING.Smoker.stopped == "Y" ~ "ex_smoker",
    EXP.SMOKING.Smoker.Now == "Y" ~ "current_smoker",
    TRUE ~ "never_smoked"))

LifeLines_x <- select(LifeLines_meta_controls, c(DAG3_sampleID, ANTHRO.Sex, ANTHRO.AGE, Smoking, ANTHRO.BMI))
LifeLines_x <- mutate(LifeLines_x, cohort = 'LifeLines')
colnames(LifeLines_x) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "cohort")
LifeLines_x$Age <- as.numeric(LifeLines_x$Age)
str(LifeLines_x)

rbind(Lynch_x, LifeLines_x) -> df_matching
df_matching[sapply(df_matching, is.character)] <- lapply(df_matching[sapply(df_matching, is.character)], as.factor)
str(df_matching)

#### Matching with MatchIt controls ----
matchit(cohort ~ Age + Sex + BMI + Smoking , data = df_matching, method = "nearest", ratio = 6) -> matched
match.data(matched) -> df_matched
summary(matched, un = FALSE)
plot(matched, type = "density")

LL_participant_ids_control <- data.frame(Participant_ID = df_matched$Participant_ID[grep("^LL.*", df_matched$Participant_ID)])
#write.csv(LL_participant_ids_control, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Matched_IDs_816.csv" )

#### Re-assign neoplasia type (aligning with Lynch definitions) ----
PALGA <- read_excel("Copy of METADATA_PALGA_data_v7.xlsx", na = "NA")
palga_1 <- PALGA 

palga_1 <- palga_1 %>%
  mutate(new_neoplasia_type = case_when(
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.size.1cm == 1) ~ "Adv.Adenoma",
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.Dysplasia == "high") ~ "Adv.Adenoma",
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.Growth_type %in% c("TVA", "VA")) ~ "Adv.Adenoma",
    (Polyp.Histological_subtype %in% c("adenoma")) ~ "Adenoma",
    TRUE ~ NA))

palga_1 <- palga_1 %>%
  mutate(new_neoplasia_type = case_when(
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("HPP") & Polyp.size.1cm == 1) ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("TSA")) ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("SSP") & Polyp.Dysplasia == "yes") ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("SSP") & Polyp.size.1cm == 1) ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated")) ~ "Serrated",
    TRUE ~ new_neoplasia_type))

palga_1 <- palga_1 %>%
  mutate(new_neoplasia_type = case_when(
    (Lesion.malign.CRC == 1) ~ "CRC",
    TRUE ~ new_neoplasia_type))

#### Prepare palga data to only cases with hit within one year after stool sampling ----
palga_2 <- palga_1 %>% filter(DATASET == "PALGA_CASE") #filter only the rows with cases
palga_2 <- palga_2 %>% filter(META.deltaDate.ba == "After.sampling") #get the samples after fecal collection
palga_2 <- palga_2 %>% filter(META.deltaDate.days <= 365) #106 #get the samples with PALGA data within 1 year after fecal sampling
palga_2 <- palga_2 %>% filter(Polyp.Neoplasia.T != "Serrated") #83 #remove the serrated polyp cases (as these are not neoplasia in Lynch cohort)
unique(palga_2$DAG3_ID) #55 unique cases

#select cases that are not excluded in DAG3 (same exclusion criteria as Lynch)
palga_2 <- palga_2 %>% filter(DAG3_ID %in% LifeLines_meta_gp$DAG3_sampleID) #57
unique(palga_2$DAG3_ID) #34 unique cases

#Decide the most advance case for the IDs that have multiple biopsies
levels_order <- c("Adenoma", "Adv.Serrated", "Adv.Adenoma", "CRC")
palga_2$new_neoplasia_type <- factor(palga_2$new_neoplasia_type, levels = levels_order, ordered = TRUE)

palga_2 <- palga_2 %>% arrange(DAG3_ID, new_neoplasia_type) # Arrange the dataframe by ID and the advancement levels

# Identify the most advanced row for each unique ID
palga_2 <- palga_2 %>% group_by(DAG3_ID) %>%
  mutate(Most_Advanced = new_neoplasia_type == max(new_neoplasia_type)) %>% ungroup()

palga_3 <- palga_2 %>% filter(Most_Advanced == TRUE) %>% distinct(DAG3_ID, .keep_all = TRUE) #34 only keep the most advanced for each ID, and remove duplicates with same level of neoplasia
palga_3 <- palga_3 %>% filter(META.deltaDate.days >= 0) #remove one sample that is not taken after stool collection #54
table(palga_3$new_neoplasia_type)

#### Merging metadata from both cohorts for matching general population----
Lynch_x <- select(meta_Lynch_baseline, c(Participant_ID, Sex, BMI, Age, Smoking, AllNeoplasia_controls))
Lynch_x <- mutate(Lynch_x, cohort = 'Lynch')
str(Lynch_x)

LifeLines_meta_gp <- LifeLines_meta_gp %>%
  mutate(Smoking = case_when(
    EXP.SMOKING.Smoker.stopped == "Y" ~ "ex_smoker",
    EXP.SMOKING.Smoker.Now == "Y" ~ "current_smoker",
    TRUE ~ "never_smoked"))

LifeLines_x <- select(LifeLines_meta_gp, c(DAG3_sampleID, ANTHRO.Sex, ANTHRO.AGE, Smoking, ANTHRO.BMI))
LifeLines_x$AllNeoplasia_controls <- ifelse(LifeLines_x$DAG3_sampleID %in% palga_3$DAG3_ID, "Neoplasia", NA)
LifeLines_x$AllNeoplasia_controls <- ifelse(LifeLines_x$DAG3_sampleID %in% LL_participant_ids_control$Participant_ID, "Control", LifeLines_x$AllNeoplasia_controls)
LifeLines_x$AllNeoplasia_controls[is.na(LifeLines_x$AllNeoplasia_controls)] <- "Not_available"
table(LifeLines_x$AllNeoplasia_controls)

LifeLines_x <- mutate(LifeLines_x, cohort = 'LifeLines')
colnames(LifeLines_x) <- c("Participant_ID", "Sex", "Age", "Smoking", "BMI", "AllNeoplasia_controls", "cohort")
LifeLines_x$Age <- as.numeric(LifeLines_x$Age)
str(LifeLines_x)

rbind(Lynch_x, LifeLines_x) -> df_matching
df_matching[sapply(df_matching, is.character)] <- lapply(df_matching[sapply(df_matching, is.character)], as.factor)
str(df_matching)

#### Matching with MatchIt general population ----
matchit(cohort ~ Age + Sex + BMI + Smoking + AllNeoplasia_controls , data = df_matching, method = "optimal", ratio = 6) -> matched
match.data(matched) -> df_matched
summary(matched, un = FALSE)
plot(matched, type = "density")

LL_participant_ids_gp <- data.frame(Participant_ID = df_matched$Participant_ID[grep("^LL.*", df_matched$Participant_ID)])
#write.csv(LL_participant_ids_gp, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Matched_IDs_1188.csv" )

#### Check neoplasia in matched samples ----
palga_check <- PALGA %>% filter(DATASET == "PALGA_CASE") 
palga_check <- palga_check %>% filter(META.deltaDate.ba == "After.sampling") 
palga_check <- palga_check %>% filter(META.deltaDate.days <= 1825) 
table(palga_check$Polyp.Neoplasia.T)

palga_check2 <- palga_check %>% filter(DAG3_ID %in% LL_participant_ids_gp$Participant_ID)
table(palga_check2$Polyp.Neoplasia.T)

#### Get dataframes to work with general population ----
General_population_IDs <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Matched_IDs_1188.csv")

#get taxa data
General_population_metaphlan <- tax4[rownames(tax4) %in% General_population_IDs$Participant_ID,] #792 samples
General_population_metaphlan <- General_population_metaphlan[,colSums(General_population_metaphlan)>0] #from 9421 to 7140
saveRDS(General_population_metaphlan, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_metaphlan_complete.rds") 

#get meta data
General_population_metadata <- LifeLines_meta_gp[LifeLines_meta_gp$DAG3_sampleID %in% General_population_IDs$Participant_ID, ]
saveRDS(General_population_metadata, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_metadata_complete.rds") 

#get pathway data
General_population_pathways <- DAG3_pathways_cleaned[rownames(DAG3_pathways_cleaned) %in% General_population_IDs$Participant_ID,] #792 samples
General_population_pathways <- General_population_pathways[,colSums(General_population_pathways)>0] #from 647 to 574
saveRDS(General_population_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/GeneralPopulation_pathways_complete.rds") 

#### Get dataframes to work with controls ----
ControlsLL_IDs <- read_csv("~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Matched_IDs_816.csv")

#get taxa data
ControlsLL_metaphlan <- tax4[rownames(tax4) %in% ControlsLL_IDs$Participant_ID,] #816 samples
ControlsLL_metaphlan <- ControlsLL_metaphlan[,colSums(ControlsLL_metaphlan)>0] #from 9421 to 7014
saveRDS(ControlsLL_metaphlan, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_metaphlan_complete.rds") 

#get meta data
ControlsLL_metadata <- LifeLines_meta_controls[LifeLines_meta_controls$DAG3_sampleID %in% ControlsLL_IDs$Participant_ID, ]
saveRDS(ControlsLL_metadata, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_metadata_complete.rds") 

#get pathway data
ControlsLL_pathways <- DAG3_pathways_cleaned[rownames(DAG3_pathways_cleaned) %in% ControlsLL_IDs$Participant_ID,] #544 samples
ControlsLL_pathways <- ControlsLL_pathways[,colSums(ControlsLL_pathways)>0] #from 647 to 509
saveRDS(ControlsLL_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/Controls_LL_pathways_complete.rds") 

#### Get dataframes to work with PALGA cases ----
#get taxa data
PALGAcases_metaphlan <- tax4[rownames(tax4) %in% palga_3$DAG3_ID,] #34
PALGAcases_metaphlan <- PALGAcases_metaphlan[,colSums(PALGAcases_metaphlan)>0] #from 9421 to 3864
saveRDS(PALGAcases_metaphlan, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_metaphlan_complete.rds") 

#get meta data
PALGAcases_metadata <- LifeLines_meta_all %>% filter(DAG3_sampleID %in% palga_3$DAG3_ID) %>%
  left_join(palga_3[, c("DAG3_ID", "new_neoplasia_type")], by = c("DAG3_sampleID" = "DAG3_ID"))
saveRDS(PALGAcases_metadata, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_metadata_complete.rds") 

#get pathway data
PALGAcases_pathways <- DAG3_pathways_cleaned[rownames(DAG3_pathways_cleaned) %in% palga_3$DAG3_ID,] #34 samples
PALGAcases_pathways <- PALGAcases_pathways[,colSums(PALGAcases_pathways)>0] #from 647 to 465
saveRDS(PALGAcases_pathways, "~/Documents/MDPhD/Hfst_IBD_LynchSyndrome/Processed_data/PALGAcases_pathways_complete.rds") 

#Extra - Test Rose's definition overlap ----
palga_1 <- PALGA 

palga_1 <- palga_1 %>%
  mutate(neoplasia_type_Rose = case_when(
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.size.1cm == 1) ~ "Adv.Neoplasia",
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.Dysplasia == "high") ~ "Adv.Neoplasia",
    (Polyp.Histological_subtype %in% c("adenoma") & Polyp.Growth_type %in% c("VA")) ~ "Adv.Neoplasia",
    (Polyp.Histological_subtype %in% c("adenoma")) ~ "Neoplasia",
    TRUE ~ NA))

palga_1 <- palga_1 %>%
  mutate(neoplasia_type_Rose = case_when(
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("HPP")) ~ "Serrated",
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("SSP", "TSA") & Polyp.size.1cm == 1) ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated") & Polyp.Growth_type %in% c("SSP", "TSA") & Polyp.Dysplasia == "yes") ~ "Adv.Serrated",
    (Polyp.Histological_subtype %in% c("serrated")) ~ "Serrated",
    TRUE ~ neoplasia_type_Rose))

palga_1 <- palga_1 %>%
  mutate(neoplasia_type_Rose = case_when(
    (Lesion.malign.CRC == 1) ~ "Adv.Neoplasia.CRC",
    TRUE ~ neoplasia_type_Rose))

palga_1 <- palga_1 %>%
  mutate(overlap = Polyp.Neoplasia.T == neoplasia_type_Rose)






