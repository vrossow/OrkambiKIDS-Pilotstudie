####    Orkambi Kids Phyloseq Creation    ####
## Virginia Rossow


getwd()

library(pacman)
pacman::p_load(tidyverse, phyloseq, magrittr, janitor)
pacman::p_load(microbiome, textshape, Biostrings)

Metadata_filt_2 <- readRDS("Metadata_filt_2.rds")
table_decontam_tbl_comp <- read_delim("table_decontam_clean.csv", delim = NULL, col_names = TRUE, col_types = NULL) 
taxonomy_decontam_tbl <- read_delim("19_ASV_taxonomy_decontam.csv", delim = NULL, col_names = TRUE, col_types = NULL) 


## tax_table
taxonomy_decontam_clean_tbl <- taxonomy_decontam_tbl %>% clean_names()
tax_clean_df <- as.data.frame(taxonomy_decontam_clean_tbl)
rownames(tax_clean_df) <- tax_clean_df$asv
tax_clean_df$asv <- NULL
tax_matrix <- as.matrix(tax_clean_df)

## otu_table
table_decontam_tbl_comp$...1 <- NULL
table_decontam_df <- as.data.frame(table_decontam_tbl_comp)
rownames(table_decontam_df) <- table_decontam_df$...2
table_decontam_df$...2 <- NULL
#req_rows <- Metadata_filt_2$x_sample_id
#table_decontam_small_df <- table_decontam_df[rownames(table_decontam_df)%in%req_rows, ]
#table_decontam_small_df[!complete.cases(table_decontam_small_df), ]


## samp
Metadata_filt_2_df <- as.data.frame(Metadata_filt_2)
rownames(Metadata_filt_2_df) = Metadata_filt_2_df$x_sample_id
Metadata_filt_2_df$x_sample_id <- NULL
#Metadata_filt_2 <- Metadata_clean[Metadata_clean$x_sample_id%in%rownames(table_decontam_df), ]
#rownames(Metadata_filt_2) = Metadata_filt_2$x_sample_id
#Metadata_filt_2_df <- as.data.frame(Metadata_filt_2)
#rownames(Metadata_filt_2_df) = Metadata_filt_2_df$x_sample_id
#Metadata_filt_2_df$x_sample_id <- NULL


#sum(nrow(unique(Metadata_filt_2$x_sample_id))) == sum(nrow(Metadata_filt_2$x_sample_id))

#### create phyloseq
otu <- otu_table(table_decontam_df, taxa_are_rows = FALSE)
samp <- sample_data(Metadata_filt_2_df)
tax <- tax_table(tax_matrix)

ps_orkambi <- phyloseq(otu, samp, tax) 
#sample_data(ps_orkambi)

# store the DNA sequences of our ASVs in the refseq slot of the phyloseq object
# and then rename our taxa to a short string, refseq(ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps_orkambi))
names(dna) <- taxa_names(ps_orkambi)
ps_orkambi <- merge_phyloseq(ps_orkambi, dna)
taxa_names(ps_orkambi) <- paste0("ASV", seq(ntaxa(ps_orkambi)))


### filtered phyloseqs

# only visit 1-9
ps_orkambi_trimmed <- subset_samples(ps_orkambi, visit != 10) 
ps_orkambi_trimmed <- subset_samples(ps_orkambi_trimmed, visit != 11)
ps_orkambi_trimmed <- subset_samples(ps_orkambi_trimmed, visit != 12)
ps_orkambi_trimmed <- subset_samples(ps_orkambi_trimmed, visit != 13)
#sample_data(ps_orkambi_trimmed) [ ,1:5]

# trim unused taxa away
ps_orkambi_trimmed <- prune_taxa(taxa_sums(ps_orkambi_trimmed) > 0, ps_orkambi_trimmed)

ps_orkambi_trimmed <- prune_samples(sample_sums(ps_orkambi_trimmed) <= 100000, ps_orkambi_trimmed)

# relative abundance ps
ps_orkambi_trimmed_relative <- transform_sample_counts(ps_orkambi_trimmed , function(x) x/sum(x) )


#### save as rds
# complete
saveRDS(ps_orkambi, "ps_orkambi.rds")
# trimmed 
saveRDS(ps_orkambi_trimmed, "ps_orkambi_trimmed.rds")
# trimmed relative abundance
saveRDS(ps_orkambi_trimmed_relative, "ps_orkambi_trimmed_relative.rds")


