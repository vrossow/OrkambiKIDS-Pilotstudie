
# Rebecca Knoll (based on Sarah Friedlmaier)
# adjusted by Virginia Rossow

# Create phyloseq object from ngless species file
# OrkambiKIDS



library(pacman)
pacman::p_load(tidyverse, dplyr, janitor, broom, phyloseq)

#load(url("https://github.com/AlessioMilanese/motus_taxonomy/blob/master/data/motus_taxonomy_2.6.1.Rdata?raw=true"))

metadata_tbl <- read_delim("Metadata_clean.csv", 
                           delim = NULL, col_names = TRUE, col_types = NULL)

wgs_species <- read_delim("orkambi_species_wg.csv", 
                          delim = NULL, col_names = TRUE, col_types = NULL)
wgs_species <- wgs_species %>%
  select(-c(1))


#### otu table

wgs_species_df <- as.data.frame(wgs_species)
# set rownames
rownames(wgs_species_df)<- wgs_species_df$species_names
wgs_species_df$species_names <- NULL
# omit NAs
wgs_species_df <- na.omit(wgs_species_df)
# Remove OTUs with 0 abundance
wgs_species_df <- wgs_species_df[rowSums(wgs_species_df) > 0,] 
32266+1305
# Extract samples with 0 reads
#wgs_species_removed <- wgs_species_df[,colSums(wgs_species_df) <= 0]
# Remove samples with 0 reads
wgs_species_df <- wgs_species_df[,colSums(wgs_species_df) > 0]
# visualize total read count per sample
hist(colSums(wgs_species_df))

wgs_species_df <- wgs_species_df %>%
  select(-c(mlv5st, sav1sta)) # the ones with more reads are kept in the df

# rename the one sample with a in the back -> better not, wrong data might end up together
#colnames(wgs_species_df)[colnames(wgs_species_df) == "mlv5sta"]<-"mlv5st"

#sort(colnames(wgs_species_df)) == sort(rownames(metadata_phy_df))
# shv8 aus dem imp table holen und hinzufugen bei wgs (IMP40)

#grep("mlv5st", colnames(wgs_species_df), value=TRUE) # in metadata schauen, wie die da genannt sind


rownames(wgs_species_df)  <- rownames(wgs_species_df) %>% 
  stringr::str_extract(pattern = "\\[.+\\]|^unassigned$") %>% #otu names within brackets
  stringr::str_remove_all(pattern = "\\[|\\]")

#### tax table

motus_wgs_adj <- wgs_species_df
# Extract matching OTUs from otu table
taxa_wgs <- motus2.6_taxonomy[motus2.6_taxonomy$mOTUs_ID %in% rownames(motus_wgs_adj), ]
# Adjust row names
rownames(taxa_wgs) <- taxa_wgs$mOTUs_ID
# Remove ID column
taxa_wgs$mOTUs_ID <- NULL
# Clear numbers from tax names
taxa_wgs[] <- lapply(taxa_wgs, gsub, pattern='[[:digit:]]+ ', replacement='') #brackets keep df format
head(taxa_wgs, n=10)


#### metadata preparation

metadata_tbl <- metadata_tbl %>% select(-c(...1))

metadata_tbl_interim <- metadata_tbl %>%
  select(-c(weight_kg, length_cm, bmi, visit_date, be_date, throat_swab_date, stool_date, urine_date, start_orkambi, sp_o2_percent, fev1,
            pred_fev1, fvc, pred_fvc, pred_lci, pp_lci, frc, pred_frc, filtratweight_g, filtratvolume_ml,
            kreon_average_i_eper_meal, sample_name, x1_calp_messung, x2_calp_messung, gew_zscore, groe_zscore, bmi_zscore)) %>%
  select(-c(label_qpcr, well, target_name, amp_status, task, cq_mean, cq_standard_deviation, quantity, quantity_mean,
            quantity_standard_deviation, material.x, label.x, baseline_start, grouping.x, material_qpcr, dna_quant_ng_ul.x,
            total_ng_dna.x, amount_stool_mg, dna_ng_per_mg_stool, grouping.y,
            intercept, r_squared, slope, efficiency, auto_threshold, auto_baseline, omit, date_q_pcr, position.x, unnamed_6,
            extraction_date, comments_katja.x, project, position.y, kit, sample_or_control, mec_eval, target_region, comments_katja.y,
            unfrozen, sex, mother, location, run_num, study_sample, interim_sample)) %>%
  mutate(label.y = gsub("_","",label.y)) %>%
  mutate(label.y = str_to_lower(label.y)) %>%
  mutate(material.y = str_to_lower(material.y)) %>%
  clean_names()

metadata_tbl_interim$sputum_staphylococcus_aureus 
# column if mrsa is there (0/1)
metadata_tbl_interim <- metadata_tbl_interim %>%
  mutate(mrsa_staph_aureus_sputum = case_when(
    str_detect(sputum_staphylococcus_aureus, "2 MRSA") == TRUE ~ "1",
    TRUE ~ "0"))

metadata_tbl_interim$sputum_staphylococcus_aureus[metadata_tbl_interim$sputum_staphylococcus_aureus == "2 MRSA"]<-"2"

# remove samples with less information
metadata_tbl_interim <- metadata_tbl_interim[metadata_tbl_interim$x_sample_id != "19MLV5St", ]
metadata_tbl_interim$label_y[metadata_tbl_interim$label_y == "mlv5st"]<-"mlv5sta"
metadata_tbl_interim <- metadata_tbl_interim[metadata_tbl_interim$x_sample_id != "19SAV1St", ]
metadata_tbl_interim$x_sample_id[metadata_tbl_interim$x_sample_id == "19SAV1StA"]<-"19SAV1St"


#  factor - like in orkambi_kids_clean.R
metadata_tbl_interim <- metadata_tbl_interim %>%
  mutate(mrsa_staph_aureus_sputum = as_factor(mrsa_staph_aureus_sputum)) %>%
  mutate(sputum_staphylococcus_aureus = as.numeric(sputum_staphylococcus_aureus)) %>%
  mutate(visit = as_factor(visit), patient_id = as_factor(patient_id), id = as_factor(id))
  #mutate(across(na_cl_6_percent:ubelkeit, as_factor)) %>%
  #mutate(across(sputum_lautropia_mirabilis:throat_stenotrophomonas_maltophilia, as_factor)) %>%
  #rename(label = label_y, material = material_y) %>%
  #rename(total_ng_dna = total_ng_dna_y, dna_quant_ng_ul = dna_quant_ng_ul_y) %>%
  #mutate(material = as_factor(material))

metadata_tbl_interim <- metadata_tbl_interim %>%
  mutate(label = label_y, material = material_y) %>%
  mutate(material = as_factor(material)) %>%
  mutate(tota_ng_dna = total_ng_dna_y, dna_quant_ng_ul = dna_quant_ng_ul_y)

metadata_tbl_interim$label_y <- NULL
metadata_tbl_interim$material_y <- NULL
metadata_tbl_interim$dna_quant_ng_ul_y <- NULL
metadata_tbl_interim$total_ng_dna_y <- NULL


metadata_tbl_interim$label[metadata_tbl_interim$label == "mjv3st " ]<-"mjv3st"


## make column visit_2 for the small boxplots divided by timepoint
metadata_tbl_interim <- metadata_tbl_interim %>% 
  mutate(visit_2 = factor(case_when(
    visit == 1 ~ "baseline",
    visit == 2 | visit == 3 ~ "start",
    visit == 5 | visit == 6 ~ "mid", 
    visit == 8 | visit == 9 ~ "end",
    TRUE ~ NA_character_
  ), levels = c("baseline", "start", "mid", "end")))

levels(metadata_tbl_interim$visit_2)

# add gender
metadata_tbl_interim <- metadata_tbl_interim %>%
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
  relocate(gender, .before = pp_fev1) %>%
  mutate(gender = as_factor(gender))


metadata_tbl_interim <- metadata_tbl_interim %>%
  mutate(base = case_when(
    visit == 1 ~ 1,
    TRUE ~ 0
  )) %>%
  relocate(base, .before = pp_fev1) %>%
  mutate(base = as_factor(base))




#### metadata

metadata_phy_df <- as.data.frame(metadata_tbl_interim)
rownames(metadata_phy_df) <- metadata_phy_df$label


#### create phyloseq object ####

# otu_table
otu_species_fil <- wgs_species_df [,colnames(wgs_species_df)%in%rownames(metadata_phy_df)]
otu_species <- otu_table(otu_species_fil, taxa_are_rows = TRUE)
otu_species <- otu_species[, gtools::mixedsort(colnames(otu_species_fil))]
# sample_metadata
sample_wgs <- sample_data(metadata_phy_df)
sample_wgs <- sample_wgs[gtools::mixedsort(rownames(sample_wgs)),]
# Confirm IDs
all(sample_names(sample_wgs) == sample_names(otu_species)) #same to merge
# tax table
tax_wgs <- tax_table(as.matrix(taxa_wgs))

### phyloseq
ps_wgs <- merge_phyloseq(otu_species, sample_wgs, tax_wgs)
#calculate total reads
sample_data(ps_wgs)$total_reads_wgs <- colSums(otu_table(ps_wgs)) # adds library_size to metadata

#### filter phyloseq object

# only visit 1-9
ps_wgs_1_9 <- subset_samples(ps_wgs, visit != 10) 
ps_wgs_1_9 <- subset_samples(ps_wgs_1_9, visit != 11)
ps_wgs_1_9 <- subset_samples(ps_wgs_1_9, visit != 12)
ps_wgs_1_9 <- subset_samples(ps_wgs_1_9, visit != 13)
sample_data(ps_wgs_1_9) [ ,1:5]

# cutoff at 50 reads per sample
ps_wgs_trimmed <- prune_samples(sample_sums(ps_wgs_1_9) >= 50, ps_wgs_1_9)

# relative abundance ps
ps_wgs_trimmed_relative <- transform_sample_counts(ps_wgs_trimmed , function(x) x/sum(x) )


## save as rds
# complete ps
saveRDS(ps_wgs, "wgs_orkambi.rds")
# subseted ps
saveRDS(ps_wgs_trimmed, "wgs_orkambi_trimmed.rds")
# subseted relative ps
saveRDS(ps_wgs_trimmed_relative, "wgs_orkambi_trimmed_relative.rds")
