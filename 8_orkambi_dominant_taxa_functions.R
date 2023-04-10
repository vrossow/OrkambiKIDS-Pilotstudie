# functions for dominant genus of orkambi project

# Virginia Rossow, derived from Rebecca Knoll


### functions for the 16S data

## find.top.taxa2 ##
#sourced from: https://github.com/joey711/phyloseq/issues/847

find.top.taxa2 <- function(x, taxa, num) {
  require(phyloseq)
  require(magrittr)
  
  top.taxa <- tax_glom(x, taxa)
  otu <-
    otu_table(top.taxa) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j1 <-
    apply(otu, 1, sort, index.return = T, decreasing = T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1, '[[', "ix") # select for index
  
  #k <- j[!duplicated(j)] # Replaced with unique() below
  l <- data.frame(unique(tax@.Data[unlist(j2),]))
  m <- data.frame(otu@.Data[, unique(unlist(j2))])
  #s <- as.name(taxa) # This isn't necessary!
  colnames(m) = l[, taxa]
  n <- apply(m, 1, sort, index.return = T, decreasing = T) %>%
    lapply('[[', "ix") %>% # Extract index
    lapply(head, n = num) # This to returns the top x tax
  
  p <- list()
  for (i in 1:length(n)) {
    p[[i]] <- colnames(m)[n[[i]]]
  }
  m$taxa <-
    p # replacing [,taxa], but now the new column is called "taxa" instead of the inputted taxonomic rank
  return(m)
}

## count genus ##
count.genus <- function(x, num){
  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired 
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = F)
  tax <- tax_table(x)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n
  
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names
  
  m <- data.frame(unlist(j2))
  
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    separate(rowname, c("SampleID", "ASV"))%>%
    dplyr::group_by(SampleID)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    #mutate(Main_taxa= Abundance>= num)%>%
    #dplyr::mutate(Type= case_when(Main_taxa== FALSE ~ "Satellites", TRUE ~ "Colonizers"))%>%
    arrange(SampleID, desc(genus))->m
  
  m$genus[is.na(m$genus)]<- "Unassigned" ##Change NA's into Unassigned 
  m$species<- NULL
  
  rm(otu, tax, j1, j2, n)
  return(m)
} # had to change Genus to genus and Species to species, because the orkambi 16s dataset is like that

## most abundant taxa ##
most_abundant_taxa <- function(ps_glom, ps_full) {
  
  top_sputum <- find.top.taxa2(ps_glom, "genus", 1)
  top_sputum$species<- NULL
  
  rslt <- top_sputum[, "taxa"]
  dd <- matrix(unlist(rslt), nrow=1)
  colnames(dd) <- rownames(top_sputum)
  top_sputum <- t(dd)
  
  top_sputum_df <- data.frame(x1 = row.names(top_sputum), top_sputum)%>%
    mutate(dominantGenus = top_sputum)
  top_sputum_df$top_sputum <- NULL
  
  # add the dominant genus to ps
  sample_data(ps_glom)$dominant_genus <- top_sputum_df$dominantGenus
  sample_data(ps_full)$dominant_genus <- top_sputum_df$dominantGenus
  
  ## pcoa
  plot_glom <- plot_ordination(ps_glom, ordinate(ps_glom, "MDS"), color = "visit", shape="dominant_genus") +
    ggtitle("Variation by dominant Genus and visit (glom)")
  plot_full <-plot_ordination(ps_full, ordinate(ps_full, "MDS"), color = "visit", shape="dominant_genus") +
    ggtitle("Variation by dominant Genus and visit (full)")
  plot_glom_dom <- plot_ordination(ps_glom, ordinate(ps_glom, "MDS"), color = "dominant_genus") +
    #ggtitle("Variation by dominant Genus (glom)")
    stat_ellipse() +
    theme_bw()
  
  
  return(list(plot_glom, plot_full, ps_glom, ps_full, plot_glom_dom))
  
}


## clustered heatmap ##
clustered_heatmap_prep <- function(ps_glom, ps_full) {
  
  #extract ASV and relative abundances from ou_table
  tmp<- as.data.frame(otu_table(ps_glom, "genus"))
  tmp <- t(tmp)
  tmp <- as.data.frame(tmp)
  
  # to later display columns in the right order I have to rename samples and sort them
  names(tmp) <- substring(names(tmp),3)
  tmp <- tmp %>% select(sort(names(tmp)))
  tmp1<- count.genus(ps_glom)
  
  tmp1%>%
    ungroup()%>%
    dplyr::select(c(genus, ASV))%>%
    unique()-> tmp1
  
  # tmp%>%
  #   rownames_to_column()%>%
  #   dplyr::rename(ASV= rowname)%>%
  #   left_join(tmp1, by="ASV")%>%
  #   unite(ASV_Genus, c("ASV", "genus"))%>%
  #   tibble::column_to_rownames(var = "ASV_Genus")->tmp
  # 
  tmp%>%
    rownames_to_column()%>%
    dplyr::rename(ASV= rowname)%>%
    left_join(tmp1, by="ASV")%>%
    unite(ASV_Genus, c("genus", "ASV"))%>%
    tibble::column_to_rownames(var = "ASV_Genus")->tmp
  tmp <- tmp[order(row.names(tmp)), ]
  
  ## Filter prevalence of taxa for heatmap
  #tmp<- tmp[apply(tmp[,-1], 1, function(x) !all(x==0)),]##Eliminate rows with zero
  tmp<- tmp[apply(tmp[,-1], 1, function(x) !all(x<2.5)),]# Eliminates taxa with overallabundance < 2.5%
  
  # #Create clusters
  # p_clust <- hclust(dist(t(tmp)), method = "ward.D") ##Dendogram
  # 
  # as.dendrogram(p_clust) %>%
  #   plot(horiz = T)
  # 
  # p_col <- cutree(tree = p_clust, k = 2)
  # p_col  <- data.frame(cluster = ifelse(test = p_col  == 1, yes = "cluster 1", no = "cluster 2"))
  # p_col$SampleID <- rownames(p_col)
  # 
  # #create annotation for columns
  # metadata_heatmpa <- as(sample_data(ps_glom),"data.frame")
  # #rename rownames and order them
  # row.names(metadata_heatmpa) <- substring(row.names(metadata_heatmpa),3)
  # metadata_heatmpa <- metadata_heatmpa[order(row.names(metadata_heatmpa)), ]
  # col_groups <- metadata_heatmpa %>%
  #   select(c(id,visit,dominant_genus)) ##Here It is possible to add the other characteristics
  
  return(tmp)
  
}

clustered_heatmap <- function(temp_object, map_o) {
  
  map_o <- pheatmap(temp_object, cluster_rows = F, cluster_cols = T,
                    #color = colorRampPalette(c("white","#832424FF"))(50), #"#3A3A98FF",
                    border_color = NA,
                    #annotation_col = col_groups, 
                    #annotation_colors = colour_groups,
                    show_rownames = T,
                    show_colnames = T,
                    fontsize = 16,
                    fontsize_row = 16,
                    fontsize_number = 16,
                    main= "Taxonomic composition of samples")
  
  return(map_o)
  
}

clustered_heatmap_ordered <- function(temp_object, map_o) {
  
  # pheatmap ordered by id and visit
  #just set cluster_cols to false to have it without clustering
  map_o <- pheatmap(temp_object, cluster_rows = F, cluster_cols = F, 
                    #color = colorRampPalette(c("white","#832424FF"))(50), #"#3A3A98FF",
                    border_color = NA,
                    #annotation_col = col_groups, 
                    #annotation_colors = colour_groups,
                    #gaps_col=c(4,9,12,15,20,24,28,33,35,39,41,45,49,53,58,63,67,72,76,80,84,87,89,93,95,96,97,98,102,105,110),
                    show_rownames = T,
                    show_colnames = T,
                    fontsize = 16,
                    fontsize_row = 16,
                    fontsize_number = 16#,
                   # main= "Taxonomic composition of samples, ordered by visit and id"
                   )
  
  return(map_o)
  
}

#### functions for wgs ####

find.top.taxa22 <- function(x, taxa, num) {
  require(phyloseq)
  require(magrittr)
  
  top.taxa <- tax_glom(x, taxa)
  otu <-
    t(otu_table(top.taxa)) # remove the transformation if using a merge_sample object
  tax <- tax_table(top.taxa)
  j1 <-
    apply(otu, 1, sort, index.return = T, decreasing = T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1, '[[', "ix") # select for index
  
  l <- data.frame(unique(tax@.Data[unlist(j2),]))
  
  m <- data.frame(otu@.Data[, unique(unlist(j2))])
  colnames(m) = l[, taxa]
  n <- apply(m, 1, sort, index.return = T, decreasing = T) %>%
    lapply('[[', "ix") %>% # Extract index
    lapply(head, n = num) # This to returns the top x tax
  
  p <- list()
  for (i in 1:length(n)) {
    p[[i]] <- colnames(m)[n[[i]]]
  }
  m$taxa <-
    p # replacing [,taxa], but now the new column is called "taxa" instead of the inputted taxonomic rank
  return(m)
}


## most abundant taxa ##
most_abundant_taxa_wgs <- function(ps_glom, ps_full) {
  
  top_sputum <- find.top.taxa22(ps_glom, "Genus", 1)
  top_sputum$species<- NULL
  
  rslt <- top_sputum[, "taxa"]
  dd <- matrix(unlist(rslt), nrow=1)
  colnames(dd) <- rownames(top_sputum)
  top_sputum <- t(dd)
  
  top_sputum_df <- data.frame(x1 = row.names(top_sputum), top_sputum)%>%
    mutate(dominantGenus = top_sputum)
  top_sputum_df$top_sputum <- NULL
  
  # add the dominant genus to ps
  sample_data(ps_glom)$dominant_genus <- top_sputum_df$dominantGenus
  sample_data(ps_full)$dominant_genus <- top_sputum_df$dominantGenus
  
  ## pcoa
  plot_glom <- plot_ordination(ps_glom, ordinate(ps_glom, "MDS"), color = "visit", shape="dominant_genus") +
    ggtitle("Variation by dominant Genus and visit (glom)")
  plot_full <-plot_ordination(ps_full, ordinate(ps_full, "MDS"), color = "visit", shape="dominant_genus") +
    ggtitle("Variation by dominant Genus and visit (full)")
  plot_glom_dom <- plot_ordination(ps_glom, ordinate(ps_glom, "MDS"), color = "dominant_genus") +
    #ggtitle("Variation by dominant Genus (glom)")
    stat_ellipse() +
    theme_bw()
  plot_glom_v2 <- plot_ordination(subset_samples(ps_glom, visit_2 != is.na(visit_2)), ordinate(subset_samples(ps_glom, visit_2 != is.na(visit_2)), "MDS"), color = "dominant_genus", shape = "visit_2") +
    ggtitle("Variation by dominant Genus and timepoint (glom)")
  
  return(list(plot_glom, plot_full, ps_glom, ps_full, plot_glom_dom, plot_glom_v2))
  
}

## count genus ##
count.genus.wgs <- function(x, num){

  require(phyloseq)
  require(magrittr)
  #x is a phyloseq object glomed to Genus
  #num is the threshold of Relative abundance desired 
  otu <- as(otu_table(x), "matrix")
  # transpose if necessary
  if(taxa_are_rows(x)){otu <- t(otu)}
  otu <- otu_table(otu, taxa_are_rows = T)
  tax <- tax_table(x)
  # Coerce to data.frame
  n <- as.data.frame(tax)
  n%>%
    rownames_to_column()%>%
    dplyr::rename(ASV = rowname)-> n
  
  j1 <- apply(otu,1,sort,index.return=T, decreasing=T) # modifying which.max to return a list of sorted index
  j2 <- lapply(j1,'[[',"x") # select for Names
  
  m <- data.frame(unlist(j2))
  
  m%>%
    rownames_to_column()%>%
    dplyr::filter(unlist.j2.!=0)%>%
    separate(rowname, c("SampleID", "ASV"), sep = "[.]")%>%
    dplyr::group_by(SampleID)%>%
    dplyr::rename(Abundance = unlist.j2.)%>%
    dplyr::mutate(Abundance = (Abundance/1E6)*100)%>%
    left_join(n, by="ASV")%>%
    #mutate(Main_taxa= Abundance>= num)%>%
    #dplyr::mutate(Type= case_when(Main_taxa== FALSE ~ "Satellites", TRUE ~ "Colonizers"))%>%
    arrange(SampleID, desc(Genus))->m
  
  m$Genus[is.na(m$Genus)]<- "Unassigned" ##Change NA's into Unassigned 
  m$Species<- NULL
  m$profiled <- NULL
  
  rm(otu, tax, j1, j2, n)
  return(m)
  
}

## clustered heatmap ##
clustered_heatmap_wgs_prep <- function(ps_glom, ps_full) {
  
  tmp<- as.data.frame(otu_table(ps_glom, "Genus"))
  #tmp <- t(tmp)
  #tmp <- as.data.frame(tmp)
  
  # to later display columns in the right order I have to rename samples and sort them
  #names(tmp) <- substring(names(tmp),6)
  #tmp <- tmp %>% select(sort(names(tmp)))
  
  tmp1<- count.genus.wgs(ps_glom)
  
  tmp1%>%
    ungroup()%>%
    dplyr::select(c(Genus, ASV))%>%
    unique()-> tmp1
  
  tmp%>%
    rownames_to_column()%>%
    dplyr::rename(ASV= rowname)%>%
    left_join(tmp1, by="ASV")%>%
    unite(ASV_Genus, c("Genus", "ASV"))%>%
    tibble::column_to_rownames(var = "ASV_Genus")->tmp
  
  ## Filter prevalence of taxa for heatmap
  tmp<- tmp[apply(tmp[,-1], 1, function(x) !all(x<2.5)),]# Eliminates taxa with overallabundance < 2.5%
  
  #Create clusters
  p_clust <- hclust(dist(t(tmp)), method = "ward.D") ##Dendogram
  
  as.dendrogram(p_clust) %>%
    plot(horiz = T)
  
  p_col <- cutree(tree = p_clust, k = 2)
  p_col  <- data.frame(cluster = ifelse(test = p_col  == 1, yes = "cluster 1", no = "cluster 2"))
  p_col$SampleID <- rownames(p_col)
  
  # for colours, i think
  
  # #create annotation for columns
  # metadata_heatmpa <- as(sample_data(ps_glom),"data.frame")
  # #rename rownames and order them
  # row.names(metadata_heatmpa) <- substring(row.names(metadata_heatmpa),3)
  # metadata_heatmpa <- metadata_heatmpa[order(row.names(metadata_heatmpa)), ]
  # col_groups <- metadata_heatmpa %>%
  #   select(c(id,visit,dominant_genus)) ##Here It is possible to add the other characteristics
  
  return(tmp)
  
}
