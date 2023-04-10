# Virginia Rossow
# Whole Genome Data Overview

# includes preprocessing of OrkambiKIDS and IMP


library(pacman)
pacman::p_load(tidyverse, janitor)


###### read tax tables ######

phylum_df <- read_delim("mainz_wgs__motusv2.6_phylum.tsv",
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

class_df <- read_delim("mainz_wgs__motusv2.6_class.tsv",
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

order_df <- read_delim("mainz_wgs__motusv2.6_order.tsv",
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

family_df <- read_delim("mainz_wgs__motusv2.6_family.tsv",
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)

species_df <- read_delim("mainz_wgs__motusv2.6.tsv", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

###### make tables prettier ######

phylum_tbl <- as_tibble(phylum_df)
phylum_tbl <- phylum_tbl %>% rename_with(.cols = 1, ~"phylum_names") %>% clean_names()

class_tbl <- as_tibble(class_df)
class_tbl <- class_tbl %>% rename_with(.cols = 1, ~"class_names") %>% clean_names()

order_tbl <- as_tibble(order_df)
order_tbl <- order_tbl %>% rename_with(.cols = 1, ~"order_names") %>% clean_names()

family_tbl <- as_tibble(family_df)
family_tbl <- family_tbl %>% rename_with(.cols = 1, ~"family_names") %>% clean_names()

species_tbl <- as_tibble(species_df)
species_tbl <- species_tbl %>% rename_with(.cols = 1, ~"species_names") %>% clean_names()


# fixing the column names and splitting the datasets

# imp data phylum
p_imp_cols <- grep("imp", colnames(phylum_tbl), value=TRUE)

phylum_imp_tbl <- phylum_tbl %>%
  select(c(phylum_names, all_of(p_imp_cols)))

imp_phylum_cols <- replace(colnames(phylum_imp_tbl), colnames(phylum_imp_tbl) == 
          "preprocessed_sample_p04_d01_imp16v3thr_p04_d01_imp16v3thrfiltered_qc_pair", 
          "preprocessed_sample_p04_d01_imp16v3thr2_p04_d01_imp16v3thrfiltered_qc_pair")

#data.frame(table(imp_phylum_cols))

imp_phylum_cols <- substr(imp_phylum_cols, 29, 42)

# grep("imp16v3thr", imp_phylum_cols)
# imp_phylum_cols[123]

imp_phylum_cols <- replace(imp_phylum_cols, imp_phylum_cols == "imp_30_v7thr_p", "imp30v7thr_p04")
imp_phylum_cols <- gsub("_.*","",imp_phylum_cols)
imp_phylum_cols <- imp_phylum_cols[-c(1)]

colnames(phylum_imp_tbl) <- c("phylum_names", imp_phylum_cols)

# orkambi data phylum
phylum_orkambi_tbl <- phylum_tbl %>%
  select(-all_of(p_imp_cols)) %>%
  select(-c(77, 82, 162, 163)) # nk

orkambi_phylum_cols <- substr(colnames(phylum_orkambi_tbl), 29, 42)
orkambi_phylum_cols <- orkambi_phylum_cols[-c(1)]
orkambi_phylum_cols <- gsub("_.*","",orkambi_phylum_cols)

grep("mlv5st", orkambi_phylum_cols)
orkambi_phylum_cols <- replace(orkambi_phylum_cols, 51, "mlv5sta")
grep("sav1st", orkambi_phylum_cols)
orkambi_phylum_cols <- replace(orkambi_phylum_cols, 78, "sav1sta")

colnames(phylum_orkambi_tbl) <- c("phylum_names", orkambi_phylum_cols)

# phylum_orkambi_tbl$preprocessed_sample_p05_h02_nk_benzonase_treatment_p05_h02_nk_benzonase_treatmentfiltered_qc_pair
# phylum_orkambi_tbl$preprocessed_sample_p04_h12_nk31012022_p04_h12_nk31012022filtered_qc_pair
# phylum_orkambi_tbl$preprocessed_sample_p01_h03_nk100122_p01_h03_nk100122filtered_qc_pair

# imp data class
class_imp_tbl <- class_tbl %>%
  select(c(class_names, all_of(p_imp_cols)))

colnames(class_imp_tbl) <- c("class_names", imp_phylum_cols)

# orkambi class data

class_orkambi_tbl <- class_tbl %>%
  select(-all_of(p_imp_cols)) %>%
  select(-c(77, 82, 162, 163)) # nk

colnames(class_orkambi_tbl) <- c("class_names", orkambi_phylum_cols)

# imp data order
order_imp_tbl <- order_tbl %>%
  select(c(order_names, all_of(p_imp_cols)))

colnames(order_imp_tbl) <- c("order_names", imp_phylum_cols)

# orkambi order data

order_orkambi_tbl <- order_tbl %>%
  select(-all_of(p_imp_cols)) %>%
  select(-c(77, 82, 162, 163)) # nk

colnames(order_orkambi_tbl) <- c("order_names", orkambi_phylum_cols)

# imp data family
family_imp_tbl <- family_tbl %>%
  select(c(family_names, all_of(p_imp_cols)))

colnames(family_imp_tbl) <- c("family_names", imp_phylum_cols)

# orkambi family data

family_orkambi_tbl <- family_tbl %>%
  select(-all_of(p_imp_cols)) %>%
  select(-c(77, 82, 162, 163)) # nk

colnames(family_orkambi_tbl) <- c("family_names", orkambi_phylum_cols)

# imp data species
species_imp_tbl <- species_tbl %>%
  select(c(species_names, all_of(p_imp_cols)))

colnames(species_imp_tbl) <- c("species_names", imp_phylum_cols)

### imp40 is an orkambi sample as well, look for it in the imp species table
grep("imp40", colnames(species_imp_tbl), value=TRUE)

shv8_wgs_tbl <- species_imp_tbl %>%
  select(c(imp40v1thr, imp40v1st))

colnames(shv8_wgs_tbl) <- c("shv8thr", "shv8st")

# orkambi species data

species_orkambi_tbl <- species_tbl %>%
  select(-all_of(p_imp_cols)) %>%
  select(-c(77, 82, 162, 163)) # nk

colnames(species_orkambi_tbl) <- c("species_names", orkambi_phylum_cols)

species_orkambi_tbl <- bind_cols(species_orkambi_tbl, shv8_wgs_tbl)


######bring tables into shape
### exchange columns and rows (kind of flip, I guess)

# phylum imp
flipped_phylum_imp_tbl <- as_tibble(t(phylum_imp_tbl))

as.character(as.vector(flipped_phylum_imp_tbl[1,])) == phylum_names_orkambi

colnames(flipped_phylum_imp_tbl) <- phylum_names_orkambi
flipped_phylum_imp_tbl <- flipped_phylum_imp_tbl[-1,]

flipped_phylum_imp_tbl$sample_names <- imp_phylum_cols
flipped_phylum_imp_tbl <- flipped_phylum_imp_tbl %>% relocate(sample_names, .before = Acidobacteria)

# class imp
flipped_class_imp_tbl <- as_tibble(t(class_imp_tbl))

as.character(as.vector(flipped_class_imp_tbl[1,])) == class_names_orkambi

colnames(flipped_class_imp_tbl) <- class_names_orkambi
flipped_class_imp_tbl <- flipped_class_imp_tbl[-1,]

flipped_class_imp_tbl$sample_names <- imp_phylum_cols
flipped_class_imp_tbl <- flipped_class_imp_tbl %>% relocate(sample_names, .before = Acidimicrobiia)

# functions flips tables,
# wg_tbl = table with whole genome in fo of one level, sample_names_vec = sample names of said table
flip_tables <- function(wg_tbl, sample_names_vec) {

  flipped_tbl <- as_tibble(t(wg_tbl))

  colnames(flipped_tbl) <- pull(wg_tbl, 1)
  flipped_tbl <- flipped_tbl[-1,]
  flipped_tbl$sample_names <- sample_names_vec
  flipped_tbl <- flipped_tbl %>% relocate(sample_names)

  return(flipped_tbl)

}

# orkambi
flipped_phylum_orkambi_tbl <- flip_tables(phylum_orkambi_tbl, orkambi_phylum_cols)
flipped_class_orkambi_tbl <- flip_tables(class_orkambi_tbl, orkambi_phylum_cols)
flipped_order_orkambi_tbl <- flip_tables(order_orkambi_tbl, orkambi_phylum_cols)
flipped_family_orkambi_tbl <- flip_tables(family_orkambi_tbl, orkambi_phylum_cols)
flipped_species_orkambi_tbl <- flip_tables(species_orkambi_tbl, orkambi_phylum_cols)

# unique(colnames(class_imp_tbl))
# unique(imp_phylum_cols)
# data.frame(table(imp_phylum_cols)) # imp16v3thr  x  2

# imp
flipped_phylum_imp_tbl <- flip_tables(phylum_imp_tbl, imp_phylum_cols)
flipped_class_imp_tbl <- flip_tables(class_imp_tbl, imp_phylum_cols)
flipped_order_imp_tbl <- flip_tables(order_imp_tbl, imp_phylum_cols)
flipped_family_imp_tbl <- flip_tables(family_imp_tbl, imp_phylum_cols)
flipped_species_imp_tbl <- flip_tables(species_imp_tbl, imp_phylum_cols)

###### total reads 

# phylum orkambi
flipped_phylum_orkambi_tbl <- flipped_phylum_orkambi_tbl %>% mutate(across(c(2:146), as.numeric))
flipped_phylum_orkambi_tbl$total_reads <- rowSums(flipped_phylum_orkambi_tbl[, c(2:146)])

# class orkambi
flipped_class_orkambi_tbl <- flipped_class_orkambi_tbl %>% mutate(across(c(2:236), as.numeric))
flipped_class_orkambi_tbl$total_reads <- rowSums(flipped_class_orkambi_tbl[, c(2:236)])

# order orkambi
flipped_order_orkambi_tbl <- flipped_order_orkambi_tbl %>% mutate(across(c(2:404), as.numeric))
flipped_order_orkambi_tbl$total_reads <- rowSums(flipped_order_orkambi_tbl[, c(2:404)])

# family orkambi
flipped_family_orkambi_tbl <- flipped_family_orkambi_tbl %>% mutate(across(c(2:761), as.numeric))
flipped_family_orkambi_tbl$total_reads <- rowSums(flipped_family_orkambi_tbl[, c(2:761)])

# species orkambi
flipped_species_orkambi_tbl <- flipped_species_orkambi_tbl %>% mutate(across(c(2:33572), as.numeric))
flipped_species_orkambi_tbl$total_reads <- rowSums(flipped_species_orkambi_tbl[, c(2:33572)])

# phylum imp
flipped_phylum_imp_tbl <- flipped_phylum_imp_tbl %>% mutate(across(c(2:146), as.numeric))
flipped_phylum_imp_tbl$total_reads <- rowSums(flipped_phylum_imp_tbl[, c(2:146)])

# class imp
flipped_class_imp_tbl <- flipped_class_imp_tbl %>% mutate(across(c(2:236), as.numeric))
flipped_class_imp_tbl$total_reads <- rowSums(flipped_class_imp_tbl[, c(2:236)])

# order imp
flipped_order_imp_tbl <- flipped_order_imp_tbl %>% mutate(across(c(2:404), as.numeric))
flipped_order_imp_tbl$total_reads <- rowSums(flipped_order_imp_tbl[, c(2:404)])

# family imp
flipped_family_imp_tbl <- flipped_family_imp_tbl %>% mutate(across(c(2:761), as.numeric))
flipped_family_imp_tbl$total_reads <- rowSums(flipped_family_imp_tbl[, c(2:761)])

# species imp
flipped_species_imp_tbl <- flipped_species_imp_tbl %>% mutate(across(c(2:33572), as.numeric))
flipped_species_imp_tbl$total_reads <- rowSums(flipped_species_imp_tbl[, c(2:33572)])

flipped_phylum_imp_tbl <- flipped_phylum_imp_tbl[-125,]
flipped_phylum_imp_tbl$sample_names[flipped_phylum_imp_tbl$sample_names == "imp16v3thr2"]<-"imp16v3thr"

flipped_class_imp_tbl <- flipped_class_imp_tbl[-125,]
flipped_class_imp_tbl$sample_names[flipped_class_imp_tbl$sample_names == "imp16v3thr2"]<-"imp16v3thr"

flipped_order_imp_tbl <- flipped_order_imp_tbl[-125,]
flipped_order_imp_tbl$sample_names[flipped_order_imp_tbl$sample_names == "imp16v3thr2"]<-"imp16v3thr"

flipped_family_imp_tbl <- flipped_family_imp_tbl[-125,]
flipped_family_imp_tbl$sample_names[flipped_family_imp_tbl$sample_names == "imp16v3thr2"]<-"imp16v3thr"

flipped_species_imp_tbl <- flipped_species_imp_tbl[-125,]
flipped_species_imp_tbl$sample_names[flipped_species_imp_tbl$sample_names == "imp16v3thr2"]<-"imp16v3thr"

flipped_class_imp_tbl$sample_names

## add material

add_material <- function(wgs_tbl) {

  wgs_tbl <- wgs_tbl %>%
    mutate(material = case_when(
      str_detect(sample_names, "spu") == TRUE ~ "sputum",
      str_detect(sample_names, "thr") == TRUE ~ "throat",
      str_detect(sample_names, "st") == TRUE ~ "stool",
      TRUE ~ "fail"
    )) %>%
    mutate(material = as_factor(material))

  return(wgs_tbl)

}

# orkambi
flipped_phylum_orkambi_tbl <- add_material(flipped_phylum_orkambi_tbl)
flipped_class_orkambi_tbl <- add_material(flipped_class_orkambi_tbl)
flipped_order_orkambi_tbl <- add_material(flipped_order_orkambi_tbl)
flipped_family_orkambi_tbl <- add_material(flipped_family_orkambi_tbl)
flipped_species_orkambi_tbl <- add_material(flipped_species_orkambi_tbl)

# imp
flipped_phylum_imp_tbl <- add_material(flipped_phylum_imp_tbl)
flipped_class_imp_tbl <- add_material(flipped_class_imp_tbl)
flipped_order_imp_tbl <- add_material(flipped_order_imp_tbl)
flipped_family_imp_tbl <- add_material(flipped_family_imp_tbl)
flipped_species_imp_tbl <- add_material(flipped_species_imp_tbl)

### saving the flipped tables as .csv

write.csv (flipped_phylum_orkambi_tbl, "flipped_phylum_orkambi.csv")
write.csv (flipped_class_orkambi_tbl, "flipped_class_orkambi.csv")
write.csv (flipped_order_orkambi_tbl, "flipped_order_orkambi.csv")
write.csv (flipped_family_orkambi_tbl, "flipped_family_orkambi.csv")
write.csv (flipped_species_orkambi_tbl, "flipped_species_orkambi.csv")

write.csv (flipped_phylum_imp_tbl, "flipped_phylum_imp.csv")
write.csv (flipped_class_imp_tbl, "flipped_class_imp.csv")
write.csv (flipped_order_imp_tbl, "flipped_order_imp.csv")
write.csv (flipped_family_imp_tbl, "flipped_family_imp.csv")
write.csv (flipped_species_imp_tbl, "flipped_species_imp.csv")

# saving the tables as .csv

write.csv (phylum_orkambi_tbl, "orkambi_phylum_wg.csv")
write.csv (class_orkambi_tbl, "orkambi_class_wg.csv")
write.csv (order_orkambi_tbl, "orkambi_order_wg.csv")
write.csv (family_orkambi_tbl, "orkambi_family_wg.csv")
write.csv (species_orkambi_tbl, "orkambi_species_wg.csv")

write.csv (phylum_imp_tbl, "imp_phylum_wg.csv")
write.csv (class_imp_tbl, "imp_class_wg.csv")
write.csv (order_imp_tbl, "imp_order_wg.csv")
write.csv (family_imp_tbl, "imp_family_wg.csv")
write.csv (species_imp_tbl, "imp_species_wg.csv")


###### unique?
unique(orkambi_phylum_cols) == orkambi_phylum_cols

unique(phylum_orkambi_tbl$phylum_names)

phylum_orkambi_tbl %>% select(c(mlv5st, sav1st)) 

grep("mlv5st", colnames(phylum_tbl))
grep("sav1st", colnames(phylum_tbl))




