####    Orkambi Kids Data Cleaning    ####

## Virginia Rossow
 

# install and load all needed packages
library(pacman)
pacman::p_load(tidyverse, phyloseq, magrittr, janitor, microbiome, lubridate, naniar, readxl)



# load data
sample_metadata_tbl <- read_delim("19_Sample_Metadata.csv", delim = NULL, col_names = TRUE, col_types = NULL)

table_decontam_tbl <- read_delim("19_ASV_table_decontam.csv", delim = NULL, col_names = TRUE, col_types = NULL) 

#clinical_data_tbl <- read_excel("ClinicalData_OrkambiKIDS.xlsx")
#clinical_data_long_tbl <- read_delim("OrkambiKIDSlong.csv", delim = NULL, col_names = TRUE, col_types = NULL)
#clinical_data_long_tbl_1 <- read_delim("clinical_data_orkambi.csv", delim = NULL, col_names = TRUE, col_types = NULL)
clinical_data_long_tbl <- read_delim("clinical_data_orkambi_1.csv", delim = NULL, col_names = TRUE, col_types = NULL)


lci_frc_tbl <- read_delim("LCI_FRC_new.csv", delim = NULL, col_names = TRUE, col_types = NULL)

percentiles_tbl <- read_delim("perzentile_long.csv", delim = NULL, col_names = TRUE, col_types = NULL)

cf_qpcr_tbl <- read_delim("CF_QPCR_results.csv", delim = NULL, col_names = TRUE, col_types = NULL)

samp_meta_sh15v8_tbl <- read_delim("Sample_MetadataSH15V8.csv", delim = NULL, col_names = TRUE, col_types = NULL)

asv_table_sh15v8_tbl <- read_delim("asv_table_SH15V8.csv", delim = NULL, col_names = TRUE, col_types = NULL)



## add SH15V8 in the metadata

samp_meta_sh15v8_tbl_clean <- samp_meta_sh15v8_tbl %>%
  select(-c(...1, raw_reads)) %>%
  mutate(project = as.character(project))

samp_meta_sh15v8_tbl_clean$project[samp_meta_sh15v8_tbl_clean$project == "FALSE"] <- "Orkambi"

sample_metadata_tbl_comp <- bind_rows(sample_metadata_tbl, samp_meta_sh15v8_tbl_clean)

#sample_metadata_tbl_comp %>% filter(X.SampleID == "IMP_40_V1")

sample_metadata_tbl_comp$Label[sample_metadata_tbl_comp$Label == "IMP_40_V1"]<-"SH_V8_St"
sample_metadata_tbl_comp$Label[sample_metadata_tbl_comp$Label == "IMP_40_V1_Thr"]<-"SH_V8_Thr"

sample_metadata_tbl_comp$X.SampleID[sample_metadata_tbl_comp$X.SampleID == "20IMP40V1"]<-"19SHV8St"
sample_metadata_tbl_comp$X.SampleID[sample_metadata_tbl_comp$X.SampleID == "20IMP40V1Thr"]<-"19SHV8Thr"

## add SH15V8 in the asv data

asv_table_sh15v8_tbl_clean <-asv_table_sh15v8_tbl %>%
  select(-c(...1))

names(asv_table_sh15v8_tbl_clean)[names(asv_table_sh15v8_tbl_clean) == "...2"] <- "...1" 

# remove columns that have only 0 as value in sh15v8
asv_cols <- colnames(asv_table_sh15v8_tbl_clean)

data_1 <- asv_table_sh15v8_tbl_clean[1, ] 
data_2 <- asv_table_sh15v8_tbl_clean[2, ] 

data_1 <- data_1 %>%
  replace_with_na_at(.vars = asv_cols ,condition = ~.x == 0)

data_2 <- data_2 %>%
  replace_with_na_at(.vars = asv_cols ,condition = ~.x == 0)

data_11 <- data_1[ , colSums(is.na(data_1)) == 0]
data_21 <- data_2[ , colSums(is.na(data_2)) == 0]

data_3 <- bind_rows(data_11, data_21)

any(duplicated(names(data_3))) # check for duplicates in column names

# bind sh15v8 to existing table_decontam
table_decontam_tbl_comp <- bind_rows(table_decontam_tbl, data_3)

table_decontam_tbl_comp$...1[table_decontam_tbl_comp$...1 == "20IMP40V1"]<-"19SHV8St"
table_decontam_tbl_comp$...1[table_decontam_tbl_comp$...1 == "20IMP40V1Thr"]<-"19SHV8Thr"

table_decontam_tbl_comp[is.na(table_decontam_tbl_comp)] <- 0

write.csv(table_decontam_tbl_comp, "table_decontam_clean.csv")


# sample_metadata_tbl #

sample_metadata_tbl_c <- sample_metadata_tbl_comp %>% clean_names()

sample_metadata_tbl_c$material[sample_metadata_tbl_c$material == "Sputum "]<-"Sputum"

sample_metadata_clean_tbl <- sample_metadata_tbl_c %>%
  mutate(material = as_factor(material)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(sample_or_control = as_factor(sample_or_control)) %>%
  select(!extraction_date) %>%
  select(-c(visit, age))



#sample_metadata_clean_tbl %>% filter(label == "SH_V8_St")

#sample_metadata_clean_tbl%>% purrr::map(levels)


## clinical_data_tbl #
# 
# clinical_data_tbl <- clinical_data_tbl %>% clean_names()
# 
# clinical_data_tbl$sputum_staphylococcus_aureus[clinical_data_tbl$sputum_staphylococcus_aureus == "1 ( das erste mal im Sputum)"]<-"1"
# 
# clinical_data_tbl <- clinical_data_tbl %>%
#   replace_with_na_at(.vars = c("sputum_lautropia_mirabilis","sputum_candida_spezies", "sputum_staphylococcus_aureus", 
#                                "sputum_haemophilus_influenzae", "sputum_serratia_marcencens", "sputum_moraxella_catarrhalis", 
#                                "sputum_streptococcus_pyogenes", "sputum_stenotrophomonas_maltophilia", "sputum_nocardia", 
#                                "sputum_aspergillus_fumigatus", "sputum_aspergillus_versicolor", "sputum_mould", 
#                                "sputum_nonfermenter_gruppe", "sputum_penicillinum_species", "sputum_escherichia_coli",
#                                "stool_consistency", "fat_stool"), condition = ~.x == 9999) %>%
#   replace_with_na_at(.vars = c("material_89" ,"sputum_lautropia_mirabilis","sputum_candida_spezies", "sputum_staphylococcus_aureus", 
#                                "sputum_haemophilus_influenzae", "sputum_serratia_marcencens", "sputum_moraxella_catarrhalis", 
#                                "sputum_streptococcus_pyogenes", "sputum_stenotrophomonas_maltophilia", "sputum_nocardia", 
#                                "sputum_aspergillus_fumigatus", "sputum_aspergillus_versicolor", "sputum_mould", 
#                                "sputum_nonfermenter_gruppe", "sputum_penicillinum_species", "sputum_escherichia_coli", "material_105", 
#                                "throat_moraxella_catarrhalis", "throat_stenotrophomonas_maltophilia_107", "throat_staphylococcus_aureus",
#                                "throat_streptococcus_pyogenes", "throat_candida_spezies", "throat_serratia_marcescens",
#                                "throat_stenotrophomonas_maltophilia_112", "throat_aspergillus_fumigatus", "throat_haemophilus_influenzae",
#                                "throat_lautropia_mirabilis", "throat_nocardia", "throat_aspergillus_versicolor",
#                                "throat_enterobacter_clocae", "throat_escherichia_coli"
#   ), condition = ~.x == "NA") 
# 
# # as.numeric
# clinical_data_clean_tbl <- clinical_data_tbl %>%
#   mutate(across (weight_kg:bmi, as.numeric)) %>%
#   mutate(across(sp_o2_percent:pankreasenzyme, as.numeric)) %>%
#   mutate(across (na_cl_6_percent:pancreatic_enzymes, ~ifelse(.==0, 0, 1))) %>%
#   mutate(across (probiotic:luffanest_homeopath, ~ifelse(.==0, 0, 1))) %>%
#   mutate(across (ne:ab15courses, as.numeric)) %>%
#   mutate(throat_stenotrophomonas_maltophilia_107 = as.numeric(throat_stenotrophomonas_maltophilia_107)) %>%
#   mutate(throat_stenotrophomonas_maltophilia_112 = as.numeric(throat_stenotrophomonas_maltophilia_112))
# 
# clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia <- apply(cbind(clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia_107,
#                                                                            clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia_112),
#                                                                      1,function(x)  ifelse(all(is.na(x)),NA,sum(x,na.rm=T)))
# 
# clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia [clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia == 2] <- 1
# clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia_107 <- NULL
# clinical_data_clean_tbl$throat_stenotrophomonas_maltophilia_112 <- NULL
# 
# cols_medication <- c("na_cl_6_percent", "na_cl_3_percent", "salbutamol", "salmeterol", "cortison_inhal", "cortison_nasenspray" , 
#                      "fenoterol", "ipatropium", "pulmozyme_dn_aase" , "dexamethason_infectodexakrupp_b_b", "aspecton_nasenspray", 
#                      "vitamines", "omeprazole10mg", "ursodeoxycholsaure", "pancreatic_enzymes", "probiotic", "movicol", 
#                      "monapax_homeopath", "fresubin_nutrini", "cetirizin_b_b", "acc", "contramutan_eupatorium_perfoliatum_echinacea",
#                      "angocin_tropaeoli_herbae_armoraciae_radicis", "naproxen", "luffanest_homeopath")
# 
# # as_factor
# clinical_data_clean_tbl <- clinical_data_clean_tbl %>%
#   mutate(visit = as_factor(visit)) %>%
#   mutate(patient_id = as_factor(patient_id), id = as_factor(id)) %>% 
#   mutate(stool_consistency = as_factor(stool_consistency), fat_stool = as_factor(fat_stool)) %>%
#   mutate(material_89 = as_factor(material_89)) %>%
#   mutate(material_105 = as_factor(material_105)) %>%
#   mutate(across (all_of(cols_medication), as_factor)) %>% 
#   mutate(across (probiotic:nasenbluten, as_factor)) %>%
#   mutate(across(sputum_lautropia_mirabilis:sputum_escherichia_coli, as_factor)) %>%
#   mutate(across(throat_moraxella_catarrhalis:throat_escherichia_coli, as_factor)) %>%
#   mutate(throat_stenotrophomonas_maltophilia = as_factor(throat_stenotrophomonas_maltophilia))
# 
# # change excel date to r date
# clinical_data_clean_tbl <- clinical_data_clean_tbl %>%
#   mutate(visit_date=as.numeric(visit_date), be_date=as.numeric(be_date),
#          throat_swab_date=as.numeric(throat_swab_date), stool_date=as.numeric(stool_date),
#          urine_date=as.numeric(urine_date))%>%
#   mutate(visit_date=excel_numeric_to_date(visit_date), be_date=excel_numeric_to_date(be_date), 
#          throat_swab_date=excel_numeric_to_date(throat_swab_date), stool_date=excel_numeric_to_date(stool_date),
#          urine_date=excel_numeric_to_date(urine_date))
# 
# 
# # check levels and data types:
# #clinical_data_clean_tbl%>% purrr::map(levels)


# clinical_data_long_tbl #

# a_rows <- clinical_data_long_tbl[159:160,] %>%
#   mutate(throat_stenotrophomonas_maltophilia_107 = as.numeric(throat_stenotrophomonas_maltophilia_107)) %>%
#   mutate(throat_stenotrophomonas_maltophilia_112 = as.numeric(throat_stenotrophomonas_maltophilia_112))
# 
# a_rows$throat_stenotrophomonas_maltophilia <- apply(cbind(a_rows$throat_stenotrophomonas_maltophilia_107,
#                                                           a_rows$throat_stenotrophomonas_maltophilia_112),
#                                                           1,function(x)  ifelse(all(is.na(x)),NA,sum(x,na.rm=T)))
# 
# a_rows$throat_stenotrophomonas_maltophilia [a_rows$throat_stenotrophomonas_maltophilia == 2] <- 1
# a_rows$throat_stenotrophomonas_maltophilia_107 <- NULL
# a_rows$throat_stenotrophomonas_maltophilia_112 <- NULL
# 
# a_rows %>%
#   mutate(throat_stenotrophomonas_maltophilia = as_factor(throat_stenotrophomonas_maltophilia)) %>%
#   relocate(throat_stenotrophomonas_maltophilia, .after = throat_escherichia_coli)

clinical_data_long_tbl$...1 <- NULL

# view(clinical_data_long_tbl_1 %>%
#   filter(label == "SA_V1_St"))
# 
# view(clinical_data_long_tbl_1 %>%
#   filter(label == "ML_V5_St"))

clinical_data_long_tbl <- clinical_data_long_tbl %>% clean_names()

# clinical_data_long_tbl$sweatchloride_mmol_l
# # correct faulty sweatchloride measurements
# clinical_data_long_tbl <- clinical_data_long_tbl %>%
#   mutate(sweatchloride_mmol_l = case_when(
#     ((visit == 1) & (id == 1)) ~ 106.0,
#     ((visit == 2) & (id == 1)) ~ 95.8,
#     ((visit == 5) & (id == 1)) ~ 82.5,
#     ((visit == 1) & (id == 2)) ~ 94.8,
#     ((visit == 2) & (id == 2)) ~ 78.8,
#     ((visit == 5) & (id == 2)) ~ 68.6,
#     ((visit == 1) & (id == 3)) ~ 91.3,
#     ((visit == 2) & (id == 3)) ~ 69.7,
#     ((visit == 5) & (id == 3)) ~ 79.1,
#     ((visit == 1) & (id == 4)) ~ 102.0,
#     ((visit == 2) & (id == 4)) ~ 70.0,
#     ((visit == 5) & (id == 4)) ~ 71.8,
#     ((visit == 1) & (id == 5)) ~ 105.0,
#     ((visit == 2) & (id == 5)) ~ 85.1,
#     ((visit == 5) & (id == 5)) ~ 77.8,
#     ((visit == 1) & (id == 6)) ~ 80.2,
#     ((visit == 2) & (id == 6)) ~ 70.6,
#     ((visit == 5) & (id == 6)) ~ 68.9,
#     ((visit == 1) & (id == 7)) ~ 103.0,
#     ((visit == 2) & (id == 7)) ~ 71.0,
#     ((visit == 5) & (id == 7)) ~ 76.2
#   ))

#write.csv(clinical_data_long_tbl, "/fast/AG_Forslund/Rebecca/CF_Scripts/clinical_data_orkambi_1.csv")

# 
# clinical_data_long_tbl_1 <- clinical_data_long_tbl_1 %>%
#   subset(label != "SA_V1_St") %>%
#   subset(label != "ML_V5_St")

# add a-rows to existing df
# clinical_data_long_tbl_1 <- rbind(clinical_data_long_tbl_1, a_rows)
# 
# clinical_data_long_tbl_1 <- clinical_data_long_tbl_1 %>%
#   relocate(throat_stenotrophomonas_maltophilia, .after = throat_escherichia_coli)

# remove mrsa from staph aureus
clinical_data_long_tbl$sputum_staphylococcus_aureus [clinical_data_long_tbl$sputum_staphylococcus_aureus == "2 MRSA"] <- "2"

# as.numeric
clinical_data_long_clean_tbl <- clinical_data_long_tbl %>%
  mutate(across (weight_kg:bmi, as.numeric)) %>%
  mutate(across(sp_o2_percent:pankreasenzyme, as.numeric)) %>%
  #mutate(across (na_cl_6_percent:pancreatic_enzymes, ~ifelse(.==0, 0, 1))) %>%
  #mutate(across (probiotic:luffanest_homeopath, ~ifelse(.==0, 0, 1))) %>%
  mutate(across (ne:ab15courses, as.numeric)) 
# %>%
#   mutate(ubelkeit_72 = as.numeric(ubelkeit_72)) %>%
#   mutate(ubelkeit_80 = as.numeric(ubelkeit_80))
# 
# clinical_data_long_clean_tbl_1$ubelkeit <- apply(cbind(clinical_data_long_clean_tbl_1$ubelkeit_72, clinical_data_long_clean_tbl_1$ubelkeit_80),
#                                           1,function(x)  ifelse(all(is.na(x)),NA,sum(x,na.rm=T)))
# 
# clinical_data_long_clean_tbl_1$ubelkeit_72 <- NULL
# clinical_data_long_clean_tbl_1$ubelkeit_80 <- NULL
# 
# clinical_data_long_clean_tbl_1 <- clinical_data_long_clean_tbl_1 %>%
#   relocate(ubelkeit, .after = nasenbluten) %>%
#   mutate(ubelkeit = (ceiling(ubelkeit/2))) %>%
#   mutate(ubelkeit = as_factor(ubelkeit))
# 

cols_medication <- c("na_cl_6_percent", "na_cl_3_percent", "salbutamol", "salmeterol", "cortison_inhal", "cortison_nasenspray" , 
                     "fenoterol", "ipatropium", "pulmozyme_dn_aase" , "dexamethason_infectodexakrupp_b_b", "aspecton_nasenspray", 
                     "vitamines", "omeprazole10mg", "ursodeoxycholsaure", "pancreatic_enzymes", "probiotic", "movicol", 
                     "monapax_homeopath", "fresubin_nutrini", "cetirizin_b_b", "acc", "contramutan_eupatorium_perfoliatum_echinacea",
                     "angocin_tropaeoli_herbae_armoraciae_radicis", "naproxen", "luffanest_homeopath")

# as_factor
clinical_data_long_clean_tbl <- clinical_data_long_clean_tbl %>%
  mutate(visit = as_factor(visit), material = as_factor(material)) %>%
  mutate(patient_id = as_factor(patient_id), id = as_factor(id)) %>% 
  mutate(stool_consistency = as_factor(stool_consistency), fat_stool = as_factor(fat_stool)) %>%
  #mutate(material_89 = as_factor(material_89)) %>%
  #mutate(material_105 = as_factor(material_105)) %>%
  mutate(across (all_of(cols_medication), as_factor)) %>% 
  mutate(across (probiotic:ubelkeit, as_factor)) %>%
  mutate(across(sputum_lautropia_mirabilis:sputum_escherichia_coli, as_factor)) %>%
  mutate(across(throat_moraxella_catarrhalis:throat_stenotrophomonas_maltophilia, as_factor))


# clinical_data_long_clean_tbl_1 <- clinical_data_long_clean_tbl_1 %>%
#   select(-c(material_89, material_105))


# # rename labels from another run
# clinical_data_long_clean_tbl_1 <- clinical_data_long_clean_tbl_1 %>%
#   mutate(label = case_when(id == 8 & label == "Run20" & material == "label_throat" ~ "SH_V8_Thr",
#                            id == 8 & label == "Run20" & material == "label_stool" ~ "SH_V8_St",
#                            TRUE ~ label))
# 
# # add x_sample_id
# clinical_data_long_clean_tbl_1 <- clinical_data_long_clean_tbl_1 %>%
#   mutate(x_sample_id = str_replace_all(label, "_", "")) %>%
#   mutate(x_sample_id = paste("19", x_sample_id)) %>%
#   mutate(x_sample_id = str_replace_all(x_sample_id, " ", ""))
# 
# #unique(clinical_data_long_clean_tbl_1$x_sample_id)
# 
# clinical_data_long_clean_tbl_1$label[clinical_data_long_clean_tbl_1$label == "MLV13St"]<-"ML_V13_St"
# 
# # view(clinical_data_long_clean_tbl_1 %>%
# #   filter(grepl("SA_V1_", label)))
# 
# clinical_data_long_clean_tbl_1$material[clinical_data_long_clean_tbl_1$label == "SA_V1_StA"] <- "label_stool"   



# lci_frc_tbl #

lci_frc_tbl <- lci_frc_tbl %>% clean_names()

lci_frc_clean_tbl <- lci_frc_tbl %>%
  mutate(id = as_factor(id)) %>%
  mutate(visit = as_factor(visit)) %>%
  mutate(neues_geraet = as_factor(neues_geraet)) %>%
  mutate(frc = as.numeric(sub(",", ".", frc, fixed = TRUE)))



# percentiles_tbl #

percentiles_tbl <- percentiles_tbl %>% clean_names()

cols_score <- c("gew_zscore", "groe_zscore", "bmi_zscore")

percentiles_clean_tbl <- percentiles_tbl %>%
  mutate(id = as_factor(id)) %>%
  mutate(patient_id = as_factor(patient_id)) %>%
  mutate(visit = as_factor(visit)) %>%
  mutate(gew_zscore = (str_replace(gew_zscore, ",", "."))) %>%
  mutate(groe_zscore = (str_replace(groe_zscore, ",", "."))) %>%
  mutate(bmi_zscore = (str_replace(bmi_zscore, ",", "."))) %>%
  mutate(across (all_of (c("gew_zscore", "groe_zscore", "bmi_zscore")), as.numeric))

# remove NA
percentiles_clean_tbl <- percentiles_clean_tbl %>%
  na.omit()



# cf_qpcr_tbl #

cf_qpcr_tbl <- cf_qpcr_tbl %>% clean_names()

cf_qpcr_clean_tbl <- cf_qpcr_tbl %>%
  replace_with_na_at(.vars = c("cq"), condition = ~.x == "Undetermined") %>%
  replace_with_na_at(.vars = c("cq_mean", "cq_standard_deviation", "quantity_mean", "quantity_standard_deviation"), condition = ~.x == "-")

cf_qpcr_clean_tbl$material[cf_qpcr_clean_tbl$material == "Sputum "]<-"Sputum"

cf_qpcr_clean_tbl <- cf_qpcr_clean_tbl %>%
  mutate(extraction_date = dmy(extraction_date)) %>%
  mutate(across(cq:cq_standard_deviation, as.numeric)) %>%
  mutate(across(quantity_mean:quantity_standard_deviation, as.numeric)) %>%
  mutate(material = as_factor(material)) 

cf_qpcr_clean_tbl <- cf_qpcr_clean_tbl %>%
  mutate(date_q_pcr = dmy(date_q_pcr))

# remove all rows that do not belong to our project
cf_qpcr_clean_tbl_or <- cf_qpcr_clean_tbl %>%
  filter(project == "Orkambi")


# transformations

# # to long format with only one lable column
# clinical_data <- clinical_data_clean_tbl %>%
#   pivot_longer(c(label_sputum, label_throat, label_stool))
# 
# # add age of participants to clinical data
# clinical_data <- clinical_data %>%
#   mutate(age = case_when(clinical_data$id == 1 ~ 10, clinical_data$id ==2 ~ 12, clinical_data$id ==3 ~ 10,
#                          clinical_data$id ==4 ~ 7, clinical_data$id ==5 ~ 11, clinical_data$id ==6 ~ 7, 
#                          clinical_data$id ==7 ~ 9, clinical_data$id ==8 ~ 5))
# 
# # was not possible to rename material = name and label = value otherwise, for whatever reason...
# # so here is a workaround
# clinical_data <- clinical_data %>%
#   mutate(material = name) %>%
#   mutate(material = as_factor(material)) %>%
#   mutate(label = as.character(value))%>%
#   filter(label!="NA")
# 
# # clinical data with correct col_names
# clinical_data1 <- clinical_data %>%
#   select(-c(name, value))

# write clinical_data_long_tbl
#write.csv (clinical_data_long_clean_tbl_1, "clinical_data_orkambi_1.csv")


# add percentiles_clean_tbl to clinical_data
clinical_data2 <- clinical_data_long_clean_tbl %>%
  left_join(percentiles_clean_tbl, by = c("patient_id", "visit", "id"))

# correct lables names
cf_qpcr_clean_tbl_or[cf_qpcr_clean_tbl_or$label == "MJ_V3_St ", ]$label <- "MJ_V3_St"
clinical_data2[clinical_data2$label == "MLV13St", ]$label <- "ML_V13_St"

# join clinical data with qpcr information
clinical_data3 <- clinical_data2 %>%
  full_join(cf_qpcr_clean_tbl_or, by = c("x_sample_id"), suffix = c("", "_qpcr"))

clinical_data3$project <- "Orkambi"


# clinical_data3[clinical_data3$label == "SH_V8_St", ]$x_sample_id <- "19SHV8St"
# clinical_data3[clinical_data3$label == "SH_V8_Thr", ]$x_sample_id <- "19SHV8Thr"

#anti_join(clinical_data2 ,cf_qpcr_clean_tbl_or, by="x_sample_id")


sample_metadata_clean_tbl_small <- sample_metadata_clean_tbl %>%
  select(-c(omnigene, qualitiy_control, exclusion, study_id, real_mec, band_cut_out)) %>%
  filter(project == "Orkambi")

# join clinical data with metadata
Metadata <- clinical_data3 %>%
  left_join(sample_metadata_clean_tbl_small, by=c("x_sample_id", "project")) 


# want to keep: sa12v1a und ml11v5
Metadata[159,1] <- as_factor(6)
Metadata[160,1] <- as_factor(4)
Metadata[159,2] <- "SA11"
Metadata[160,2] <- "ML11"
Metadata[159,4] <- as_factor(1)
Metadata[160,4] <- as_factor(5)

# add gender
Metadata <- Metadata %>%
  mutate(gender = case_when( # 0 ~ m, 1 ~ w
   id == 3 ~ 1,
   id == 1 ~ 1,
   id == 5 ~ 0,
   id == 2 ~ 0,
   id == 4 ~ 1,
   id == 6 ~ 1,
   id == 8 ~ 1,
   id == 7 ~ 0,
   TRUE ~ 0
  ))%>% 
  relocate(gender, .before = visit_date) %>%
  mutate(gender = as_factor(gender))


Metadata <- Metadata %>%
  mutate(base = case_when(
    visit == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  relocate(base, .before = visit_date) %>%
  mutate(base = as_factor(base))

## make column visit_2 for the small boxplots divided by timepoint
Metadata <- Metadata %>% 
  mutate(visit_2 = factor(case_when(
    visit == 1 ~ "baseline",
    visit == 2 | visit == 3 ~ "start",
    visit == 5 | visit == 6 ~ "mid", 
    visit == 8 | visit == 9 ~ "end",
    TRUE ~ NA_character_
  ), levels = c("baseline", "start", "mid", "end")))


Metadata$grouping.y[Metadata$grouping.y == "IMP40V1"]<-"SHV8St"
Metadata$grouping.y[Metadata$grouping.y == "IMP40V1Thr"]<-"SHV8Thr"

# export csv
write.csv (Metadata, "Metadata_complete.csv")



# remove unnecessary/duplicate columns
Metadata_clean <- Metadata %>%
  clean_names() %>%
  select(-c(sex, mother, location, comments_katja_x, comments_katja_y, material_qpcr, 
            dna_quant_ng_ul_y, total_ng_dna_x, label_y, label_qpcr, position_x, grouping_x)) 

# rename clolumns
names(Metadata_clean)[names(Metadata_clean) == "label_x"] <- "label"
names(Metadata_clean)[names(Metadata_clean) == "material_x"] <- "material_label"
names(Metadata_clean)[names(Metadata_clean) == "material_y"] <- "material"
names(Metadata_clean)[names(Metadata_clean) == "dna_quant_ng_ul_x"] <- "dna_quant_ng_ul"
names(Metadata_clean)[names(Metadata_clean) == "total_ng_dna_y"] <- "total_ng_dna"
names(Metadata_clean)[names(Metadata_clean) == "position_y"] <- "position"
names(Metadata_clean)[names(Metadata_clean) == "grouping_y"] <- "grouping"

# check if it worked 
nrow(Metadata_clean) - sum(is.na(Metadata_clean$x_sample_id)) == nrow(Metadata_clean)

Metadata_filt_2 <- Metadata_clean

Metadata_filt_2 <- Metadata_filt_2 %>%
  select(-c(180:186))

# want to keep: sa12v1a und ml11v5
Metadata_filt_2 <- Metadata_filt_2[-160,]
Metadata_filt_2 <- Metadata_filt_2[-159,]
Metadata_filt_2$x_sample_id[Metadata_filt_2$x_sample_id == "19MLV5StA"]<-"19MLV5St"
Metadata_filt_2$grouping[Metadata_filt_2$grouping == "SAV1St"]<-"SAV1StA"
Metadata_filt_2$label[Metadata_filt_2$label == "ML_V5_StA"]<-"ML_V5_St"

write.csv (Metadata_filt_2, "Metadata_work.csv")
saveRDS(Metadata_filt_2, "Metadata_filt_2.rds")

