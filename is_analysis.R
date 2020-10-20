library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(reshape2)
library(gridExtra)
library(tidyr)
library(vegan)
library(cluster)
library(purrr)
set.seed(1)

#### Read data ####
# Metadata
metadata <- read.csv("db/SAMPLES/metadata.csv", stringsAsFactors = FALSE)
samples <- read.table("db/SAMPLES/samples.txt", stringsAsFactors = FALSE, header = FALSE)
metadata <- metadata[metadata$ID %in% samples$V1,]
metadata$Location_sampletype <- paste(metadata$sample_type, "-", metadata$Location)
metadata <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) 
metadata_summary <- metadata %>%
  group_by(Location, sample_type, timepoint) %>%
  summarise(n_total = n_distinct(ID))

# Save metadata for combined analysis
saveRDS(metadata, file = "data/metadata_TEs.RDS")

# Insertion site - ARG data
arg_is <- readRDS("data/filtered_arg_is.RDS")
arg_is <- arg_is[arg_is$ID %in% metadata$ID,]

# Insertion site - transposase data
trsps_is <- readRDS("data/filtered_trsps_is.RDS")
trsps_is <- trsps_is[trsps_is$ID %in% metadata$ID,]

# Contig data
no_contigs <- read.delim("data/no_contigs.txt", stringsAsFactors = FALSE, header = FALSE)
names(no_contigs) <- c("ID", "no_contigs")
no_contigs <- no_contigs[no_contigs$ID %in% metadata$ID,]

no_contigs_with_is <- read.delim("data/no_contigs_with_is.txt", stringsAsFactors = FALSE, header = FALSE)
names(no_contigs_with_is) <- c("ID", "no_contigs_is")
no_contigs_with_is <- no_contigs_with_is[no_contigs_with_is$ID %in% metadata$ID,]

# Contig data with ARGs
all_args <- read.csv("../../Phageome/vcarr_project/data/all_assemblies_card_90.csv", stringsAsFactors = FALSE)
all_args$ID <- gsub("_card.out", "", gsub(".*output_files/", "", all_args$filename))
all_args$name <- paste0(all_args$ID, "_", all_args$qseqid)

arg_no_is <- readRDS("data/filtered_arg_contig.RDS")

# GIT site colours
cols <- plasma(length(unique(metadata$sample_type)), end = 0.8)
names(cols) <- sort(unique(metadata$sample_type))

# Cohort 
cohort_cols <- c("grey", brewer.pal(9, "Blues")[c(5,7)], "seagreen3", brewer.pal(9, "YlOrRd")[c(3,5,7,9)], brewer.pal(9, "RdPu")[c(3,5,7,9)])
names(cohort_cols) <- sort(unique(metadata$Location_sampletype)) 

# Location
location_cols <- brewer.pal(11, "Spectral")[c(2,10,4,11,3)]
names(location_cols) <- sort(unique(metadata$Location))

#### ARG IS and Transposase IS catalogue ####
trsps_is_cat <- trsps_is[!duplicated(paste(trsps_is$target_name, trsps_is$itr_cluster)), c("target_name", "itr_cluster")]

#### Numbers ####
# Number of samples
length(unique(no_contigs_with_is$ID))

# Number of ITRs linked to transposase
length(unique(trsps_is$itr_cluster))

# Number/percentage contigs with sites linked to transposase -> putative insertion sites
length(unique(trsps_is$contig))
length(unique(trsps_is$contig))/sum(no_contigs_with_is$no_contigs)*100

# Number of types of transposase
trsps_types <- trsps_is %>%
  group_by(target_name) %>%
  summarise(n = n_distinct(query_name)) %>%
  arrange(desc(n))
trsps_types$target_name <- factor(trsps_types$target_name, levels = trsps_types$target_name[order(trsps_types$n, decreasing = TRUE)])

tiff("figures/insertion_site_clusters_transposase.tiff", width = 2200, height = 1400, res = 250)
ggplot(trsps_types, aes(target_name, n)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("Transposase protein family") +
  ylab("No. contigs with ISs") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 10),
        axis.title = element_text(size = 16))
dev.off()

# Number of unique insertion clusters-transposase
nrow(trsps_is_cat)

# How many insertion clusters-transposase are shared between non-related samples
trsps_dist <- trsps_is %>%
  filter(ID %in% metadata$ID[metadata$Visit_Number == 1]) %>%
  group_by(target_name, itr_cluster) %>%
  summarise(n = n_distinct(ID)) %>%
  arrange(desc(n))

# Number/percentage of unique transposase-IS pairs
trsps_uniq <- trsps_dist %>%
  filter(n == 1)
nrow(trsps_uniq)
nrow(trsps_uniq)/nrow(trsps_dist)*100

# Number of shared transposase-IS pairs
trsps_shared <- trsps_dist %>%
  filter(n > 1) %>%
  mutate(name = paste0("ITRc_", as.character(itr_cluster), " - ", target_name))
nrow(trsps_shared)
nrow(trsps_shared)/nrow(trsps_dist)*100

tiff("figures/shared_insertion_site_clusters_transposases.tiff", width = 1000, height = 1000, res = 200)
ggplot(trsps_shared, aes(n)) +
  geom_histogram(breaks = seq(2,max(trsps_shared$n),2)) +
  theme_classic() +
  xlab("No. samples") +
  ylab("No. unique ITR cluster-transposase pairs") +
  scale_x_continuous(breaks = seq(0, max(trsps_shared$n), 50)) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 14))
dev.off()

trsps_is_melt <- trsps_is %>%
  filter(ID %in% metadata$ID[metadata$Visit_Number == 1]) %>%
  mutate(name = paste0("ITRc_", as.character(itr_cluster), " - ", target_name)) %>%
  inner_join(trsps_shared) %>%
  ungroup() %>%
  select(ID, name) %>%
  unique()

trsps_is_mat <- dcast(trsps_is_melt, name ~ ID, length)

trsps_is_mat <- as.data.frame(trsps_is_mat)
row.names(trsps_is_mat) <- trsps_is_mat[,1]
trsps_is_mat <- trsps_is_mat[,-1]
trsps_is_mat[trsps_is_mat > 1] <- 1

# cluster_samples_nmds <- metaMDS(t(trsps_is_mat), distance = "jaccard", k = 2, trymax = 20)
# saveRDS(cluster_samples_nmds, file = "data/IS_nmds.RDS")
cluster_samples_nmds <- readRDS("data/IS_nmds.RDS")
df_cluster_samples_nmds <- as.data.frame(cluster_samples_nmds$points)

# Silhoette analysis of PAM (k-medoids)
avg_sil <- numeric(20)
for(k in 3:(length(avg_sil)+1)) {
  tmp <- silhouette(pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], k = k), df_cluster_samples_nmds[,c("MDS1", "MDS2")])
  avg_sil[k-1] <- mean(tmp[,3])
}

# Group by silhouette width
samples_clust <- pam(df_cluster_samples_nmds[,c("MDS1", "MDS2")], which.max(avg_sil)+1)
df_cluster_samples_nmds$cluster <- as.factor(samples_clust$cluster[row.names(df_cluster_samples_nmds)])
df_cluster_samples_nmds$ID <- row.names(df_cluster_samples_nmds)
df_cluster_samples_nmds$Sample.name <- as.character(sapply(df_cluster_samples_nmds$ID, function(x) metadata$Sample.name[metadata$ID == x]))
df_cluster_samples_nmds$Location <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$Location[metadata$ID == x])
df_cluster_samples_nmds$sample_type <- sapply(df_cluster_samples_nmds$ID, function(x) metadata$sample_type[metadata$ID == x])
df_cluster_samples_nmds$Location_sampletype <- paste(df_cluster_samples_nmds$sample_type, "-", df_cluster_samples_nmds$Location)

tiff("figures/nmds_clusters.tiff", width = 2000, height = 1000, res = 250)
ggplot(df_cluster_samples_nmds[!df_cluster_samples_nmds$ID %in% c("SRS050007", "SRS475922", "ERS2692276"),], aes(MDS1, MDS2, colour = Location_sampletype, shape = cluster)) +
  theme_classic() +
  geom_text(aes(label=cluster)) +
  scale_colour_manual("GIT Site - Region", values = cohort_cols, guide = guide_legend(override.aes = list(shape = 21, size = 4))) +
  xlab("NMDS 1") + ylab("NMDS 2") 
dev.off()

# Permanova
clusters_samples_kw <- apply(trsps_is_mat, 1, function(x) kruskal.test(x, df_cluster_samples_nmds$cluster)$p.value)
clusters_samples_kw = p.adjust(clusters_samples_kw, method = "bonferroni") #Correct p values
important_is <- clusters_samples_kw[clusters_samples_kw < 0.05]

# How many insertion clusters-transposase are shared between samples from each cohort
trsps_shared_git <- trsps_is %>%
  filter(ID %in% metadata$ID[metadata$Visit_Number == 1]) %>%
  ungroup() %>%
  group_by(target_name, itr_cluster) %>%
  mutate(n_total = n_distinct(ID)) %>%
  left_join(metadata) %>%
  ungroup() %>%
  group_by(Location, sample_type, Location_sampletype, target_name, itr_cluster, n_total) %>%
  summarise(n = n_distinct(ID)) %>%
  arrange(desc(n)) %>%
  mutate(name = paste0("ITRc_", as.character(itr_cluster), " - ", target_name)) %>%
  mutate(perc = n/n_total*100)

trsps_shared_git_prev <- trsps_shared_git %>%
  filter(name %in% names(important_is)) %>%
  arrange(desc(n_total))

# How many oral, gut and both
trsps_shared_site <- trsps_shared_git_prev %>%
  mutate(site = ifelse(sample_type == "stool", "gut", "oral")) 
trsps_oral <- unique(trsps_shared_site$name[trsps_shared_site$site == "oral"])
trsps_gut <- unique(trsps_shared_site$name[trsps_shared_site$site == "gut"])
trsps_both <- intersect(trsps_oral, trsps_gut)
trsps_oral <- trsps_oral[!trsps_oral %in% trsps_both]
trsps_gut <- trsps_gut[!trsps_gut %in% trsps_both]

tiff("figures/relative_abundance_is.tiff", height = 1000, width = 3500, res = 150)
ggplot(trsps_shared_git_prev, aes(reorder(name, n_total), perc, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  #facet_grid(sample_type ~ ., space = "free", scales = "free", switch = "both") +
  theme_classic() +
  xlab("IS (ITR cluster-Transposase)") + ylab("% samples with IS") +
  #ggtitle(unique_sample_type[i]) +
  scale_fill_manual("GIT Site - Region", values = cohort_cols) +
  theme(axis.title = element_text(size = 30),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# Save p-values
df_important_is <- data.frame(name = names(important_is), pvalue = important_is, 
                              group = sapply(names(important_is), function(x) names(which.max(table(as.numeric(df_cluster_samples_nmds$cluster)[as.logical(trsps_is_mat[row.names(trsps_is_mat) == x,])])))))
df_important_is <- df_important_is %>%
  arrange(pvalue)
write.csv(df_important_is, file = "data/important_insertion_sequences.csv", row.names = FALSE)

# Longitudinal ITR clusters-transposase
metadata_longus <- metadata %>% group_by(Location, sample_type, Sample.name) %>%
  mutate(timepoint = rank(as.numeric(Visit_Number))) %>%
  filter(Location == "USA") %>% 
  group_by(Sample.name, sample_type) %>%
  filter(any(timepoint == 3))

metadata_longus_summary <- metadata_longus %>%
  group_by(sample_type) %>%
  summarise(n = n_distinct(Sample.name))

is_long <- left_join(metadata_longus, trsps_is) %>%
  mutate(name = paste(target_name, " - ITRc_", as.character(itr_cluster)))

# Count average time points
cluster_no_tp <- is_long %>% group_by(sample_type, Sample.name, name) %>%
  summarise(no_tp = n_distinct(timepoint)) 

# Count all distinct clusters
cluster_all_summary <- is_long %>% group_by(sample_type, Sample.name) %>%
  summarise(total_is = n_distinct(name))

# Count clusters in at least x number of time points
cluster_count_tp <- map_df(.x = unique(cluster_no_tp$no_tp), .f = function(.x) {
  tmp <- data.frame(sample_type = rep(cluster_no_tp$sample_type[cluster_no_tp$no_tp == .x], .x),
                    Sample.name = rep(cluster_no_tp$Sample.name[cluster_no_tp$no_tp == .x], .x),
                    name= rep(cluster_no_tp$name[cluster_no_tp$no_tp == .x]),
                    no_tp = rep(cluster_no_tp$no_tp[cluster_no_tp$no_tp == .x], .x))
  tmp$tp <- c(1:.x)
  return(tmp)
})

cluster_count_tp_summary <- cluster_count_tp %>% group_by(sample_type, Sample.name, tp) %>%
  summarise(total_cluster_least = n_distinct(name))

cluster_one <- cluster_count_tp_summary %>%
  filter(tp == 1) %>%
  rename(total_cluster = total_cluster_least) %>%
  select(-tp)

cluster_count_tp_summary <- left_join(cluster_count_tp_summary, cluster_one, by = c("sample_type", "Sample.name")) %>%
  mutate(cluster_frac = total_cluster_least/total_cluster)

# Plot
tiff("figures/longitudinal_cluster_count.tiff", height = 800, width = 2000, res = 150)
set.seed(1)
ggplot(cluster_count_tp_summary, aes(as.factor(tp), cluster_frac, fill = as.factor(sample_type))) +
  geom_boxplot() +
  geom_jitter() +
  facet_grid(~ sample_type, scale = "free", space = "free") +
  theme_classic() + xlab("No. timepoints") + ylab("Proportion of ISs (ITR cluster/transposase pairs)") +
  theme(strip.background = element_blank(), strip.text.x = element_blank(), text = element_text(size = 16)) +
  scale_fill_manual("GIT Site", values = cols)
dev.off()

# Add ARG class and mechanism
card_aro_categories_index <- read.csv("db/CARD_DB/card-data/aro_categories_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_index <- read.csv("db/CARD_DB/card-data/aro_index.csv", stringsAsFactors = FALSE, sep = "\t")
card_metadata <- full_join(card_index, card_aro_categories_index, by = "Protein.Accession")
card_metadata <- card_metadata[!duplicated(card_metadata$ARO.Accession),]
arg_is$ARO.Accession <- paste0("ARO:", gsub("_.*", "", gsub(".*_ARO:", "", arg_is$sseqid)))
arg_is <- left_join(arg_is, card_metadata, by = "ARO.Accession")
arg_is$Drug.Class.alt <- arg_is$Drug.Class
arg_is$Drug.Class.alt[sapply(arg_is$Drug.Class.alt, function(x) str_count(x, ";")) > 1] <- "multidrug"
arg_is$Drug.Class.alt[arg_is$Resistance.Mechanism == "antibiotic efflux"] <- paste(arg_is$Drug.Class.alt[arg_is$Resistance.Mechanism == "antibiotic efflux"], "efflux")

# Filter ITR clusters that have been seen to link with transposase (because contigs aren't complete)
itrs <- data.frame(itr_cluster = unique(trsps_is$itr_cluster))
trsps_arg_is <- itrs %>%
  inner_join(arg_is, by = c("itr_cluster"="is_cluster")) 
  
# Number of unique contigs with ARG-ITR pairs
length(unique(trsps_arg_is$qseqid))

# Number of samples
length(unique(trsps_arg_is$ID))

# Number of unique ARG-ITR pairs
length(unique(paste(trsps_arg_is$ARO.Name, trsps_arg_is$is_cluster)))

# Add ARG metadata
trsps_arg_is$Drug.Class.alt <- trsps_arg_is$Drug.Class
trsps_arg_is$Drug.Class.alt[sapply(trsps_arg_is$Drug.Class.alt, function(x) str_count(x, ";")) > 2] <- "multidrug"
trsps_arg_is$Drug.Class.alt[trsps_arg_is$Resistance.Mechanism == "antibiotic efflux"] <- paste(trsps_arg_is$Drug.Class.alt[trsps_arg_is$Resistance.Mechanism == "antibiotic efflux"], "efflux")
trsps_arg_is$Drug.Class.alt[trsps_arg_is$Drug.Class.alt == "macrolide antibiotic;lincosamide antibiotic;streptogramin antibiotic"] <- "MLS antibiotic"

# Save for combined analysis
trsps_arg_is_save <- trsps_arg_is %>% 
  ungroup() %>%
  select(ID, qseqid, ARO.Name, AMR.Gene.Family, Drug.Class.alt, Resistance.Mechanism, itr_cluster)
saveRDS(trsps_arg_is_save, file = "data/arg_TEs.RDS")

# Use ARG family for multiple ARGs
trsps_arg_is <- trsps_arg_is %>% group_by(ID, qseqid, qstart, qend, itr_cluster, evalue, Drug.Class.alt, AMR.Gene.Family, Resistance.Mechanism) %>%
   summarise_all(paste, collapse = ",")
trsps_arg_is$AMR.Gene.Family <- sapply(trsps_arg_is$AMR.Gene.Family, function(x) paste(unique(strsplit(x, ",")[[1]]), collapse = ","))
trsps_arg_is$ARO.Name[sapply(trsps_arg_is$ARO.Name, function(x) str_count(x, ",")) > 0] <- paste(trsps_arg_is$AMR.Gene.Family[sapply(trsps_arg_is$ARO.Name, function(x) str_count(x, ",")) > 0], "family")

# Enrichment of ARGs on contigs with insertion sites
arg_is_counts <- trsps_arg_is %>%
  left_join(metadata, by = "ID") %>%
  group_by(Location, sample_type, ID) %>%
  summarise(n_contigs_is_args = n_distinct(qseqid)) %>%
  left_join(no_contigs_with_is, by = "ID") %>%
  mutate(perc_contigs_is = n_contigs_is_args/no_contigs_is*100)

arg_no_is_counts <- arg_no_is %>%
  left_join(metadata, by = "ID") %>%
  group_by(Location, sample_type, ID) %>%
  summarise(n_contigs_no_is_args = n_distinct(qseqid)) %>%
  left_join(no_contigs, by = "ID") %>%
  left_join(no_contigs_with_is, by = "ID") %>%
  mutate(no_contigs_no_is = no_contigs - no_contigs_is) %>%
  mutate(perc_contigs_no_is = n_contigs_no_is_args/no_contigs_no_is*100)

arg_contig_counts <- inner_join(arg_is_counts, arg_no_is_counts)
arg_contig_counts$n_contigs_is_args[is.na(arg_contig_counts$n_contigs_is_args)] <- 0

# Wilcoxon Signed-Rank Test
arg_contig_counts_melt <- melt(arg_contig_counts[, c("ID", "perc_contigs_no_is", "perc_contigs_is")])
wilcox_res <- wilcox.test(value ~ variable, arg_contig_counts_melt, paired = TRUE)
wilcox_res$p.value

tiff("figures/arg_is_enrichment.tiff", width = 1000, height = 1000, res = 150)
ggplot(arg_contig_counts_melt, aes(variable, value)) +
  geom_boxplot() +
  theme_classic() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(labels = c("Contigs without ITR cluster", "Contigs with ITR cluster")) +
  xlab("") +
  ylab("% contigs with ARG") +
  annotate("text", x = 1.5, y = 1, label = "*", size = 8)
dev.off()

# Summarise associations by ARG
metadata <- metadata %>%
  group_by(Location, sample_type, Visit_Number) %>%
  mutate(n_total = n_distinct(ID))

trsps_arg_is_summary <- trsps_arg_is %>%
  left_join(metadata, by = "ID") %>%
  filter(Visit_Number == 1) %>%
  group_by(Location_sampletype, Location, sample_type, sseqid, ARO.Name, Drug.Class.alt) %>%
  summarise(n = n_distinct(ID)) %>%
  mutate(labels = paste0(Drug.Class.alt, " - ", ARO.Name)) 

trsps_arg_is_summary$labels <- factor(trsps_arg_is_summary$labels, levels = rev(levels(factor(trsps_arg_is_summary$labels))))

tiff("figures/samples_arg_is.tiff", height = 5500, width = 4000, res = 200)
ggplot(trsps_arg_is_summary, aes(labels, n, fill = Location_sampletype)) +
  geom_bar(stat = "identity") +
  #facet_grid(Drug.Class.alt ~ ., space = "free", scales = "free", switch = "both") +
  theme_bw() +
  coord_flip() +
  xlab("ARG Class - ARG") + ylab("No. samples") +
  #ggtitle(unique_sample_type[i]) +
  scale_fill_manual("GIT Site -\nRegion", values = cohort_cols, labels = names(cohort_cols)) +
  #scale_y_continuous(breaks = seq(0, max(trsps_arg_is_summary$n_total), 1)) +
  theme(strip.text.y = element_text(angle = 180, size = 8),
        axis.title = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
dev.off()

# Number of ARGs
length(unique(trsps_arg_is_summary$ARO.Name))

# Number of classes
length(unique(trsps_arg_is_summary$Drug.Class.alt))

# Number of ARGs and classes in gut
arg_oral <- unique(trsps_arg_is_summary$ARO.Name[trsps_arg_is_summary$sample_type != "stool"])
arg_gut <- unique(trsps_arg_is_summary$ARO.Name[trsps_arg_is_summary$sample_type == "stool"])
arg_both <- intersect(arg_oral, arg_gut)
arg_oral <- arg_oral[!arg_oral %in% arg_both]
arg_gut <- arg_gut[!arg_gut %in% arg_both]

arg_class_oral <- unique(trsps_arg_is_summary$Drug.Class.alt[trsps_arg_is_summary$sample_type != "stool"])
arg_class_gut <- unique(trsps_arg_is_summary$Drug.Class.alt[trsps_arg_is_summary$sample_type == "stool"])
arg_class_both <- intersect(arg_class_oral, arg_class_gut)
arg_class_oral <- arg_class_oral[!arg_class_oral %in% arg_class_both]
arg_class_gut <- arg_class_gut[!arg_class_gut %in% arg_class_both]

# Longitudinal ARGs
arg_long <- trsps_arg_is %>% left_join(metadata, by = "ID")
arg_long <- arg_long[arg_long$Sample.name %in% arg_long$Sample.name[arg_long$Visit_Number > 1],]  
arg_long <- arg_long %>% 
  group_by(sample_type, Sample.name, ARO.Name, itr_cluster) %>%
  summarise(n = n_distinct(Visit_Number)) %>%
  group_by(sample_type, ARO.Name, itr_cluster, n) %>%
  summarise(freq = n_distinct(Sample.name)) %>%
  group_by(sample_type, ARO.Name, n) %>%
  summarise(total = sum(freq)) %>%
  group_by(sample_type) %>%
  filter(ARO.Name %in% ARO.Name[n > 1])

tiff("figures/longitudinal_args.tiff", width = 3000, height =3000, res = 250)
ggplot(arg_long, aes(as.character(n), total, fill = (as.character(n)))) +
  geom_bar(stat = "identity") +
  scale_fill_manual("No. timepoints", values = brewer.pal(11, "Spectral")[c(2,10,4)]) +
  facet_grid(sample_type + ARO.Name ~ ., scale = "free", space = "free", switch = "both") +
  theme_bw() +
  coord_flip() +
  theme(strip.text.y = element_text(angle = 180, size = 8),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab("No. individuals") + xlab("") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()  

# Number of ARG-ITRs in at least 2 timepoints
length(unique(paste(arg_long$itr_cluster, arg_long$ARO.Name)))

# Number of ARG-ITRs in USA samples
arg_us <- trsps_arg_is %>% 
  left_join(metadata, by = "ID") %>%
  filter(Sample.name %in% Sample.name[Visit_Number > 1])
length(unique(paste(arg_us$itr_cluster, arg_us$ARO.Name)))

# # Prevalence of drug classes
# trsps_arg_is %>%
#   left_join(metadata, by = "ID") %>%
#   filter(Visit_Number == 1) %>%
#   group_by(Location_sampletype, Location, sample_type, sseqid, ARO.Name, Drug.Class.alt, itr_cluster) %>%
#   summarise(n_cohort = n_distinct(ID)) %>%
#   ungroup() %>%
#   group_by(sample_type, sseqid, ARO.Name, Drug.Class.alt, itr_cluster) %>%
#   mutate(n_total = sum(n_cohort)) %>%
#   group_by(sample_type) %>%
#   filter(n_total > 1) 
# 
# # Summary of number of contigs per sample for each category
# no_contigs_is_trsps <- trsps_is %>%
#   group_by(ID) %>%
#   summarise(no_contigs_is_trsps = n_distinct(contig))
# 
# no_contigs_is_args <- arg_is %>%
#   group_by(ID)%>%
#   summarise(no_contigs_is_args = n_distinct(qseqid))
# 
# no_contigs_is_trsps_args <- trsps_arg_is %>%
#   group_by(ID) %>%
#   summarise(no_contigs_is_trsps_args = n_distinct(qseqid))
# 
# contig_dist <- no_contigs %>% left_join(no_contigs_with_is, by = "ID") %>%
#   left_join(no_contigs_is_trsps, by = "ID") %>%
#   left_join(no_contigs_is_args, by = "ID") %>%
#   left_join(no_contigs_is_trsps_args, by = "ID") %>%
#   melt() %>%
#   mutate(value = replace(value, is.na(value), 0)) %>%
#   left_join(metadata, by = "ID") %>%
#   mutate(log_count = log10(value+1))
# 
# tiff("figures/contig_counts.tiff", width = 3000, height = 1500, res = 200)
# ggplot(contig_dist, aes(variable, log_count, fill = sample_type)) +
#   geom_boxplot() +
#   theme_classic() +
#   facet_grid(~Location, scale = "free", space = "free") +
#   scale_fill_manual("GIT Site", values = cols) +
#   ylab("log10(no. contigs+1)") + xlab("") +
#   theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
#   scale_x_discrete(labels = c("All contigs", "Contigs w/ insertion site(s)", 
#                               "Contigs w/ insertion site(s) &\n transposase(s)", 
#                               "Contigs w/ insertion site(s) &\n ARG(s)", 
#                               "Contigs w/ insertion sites(s),\n transposase(s) & ARG(s)"))
# dev.off()  
# 
# 
# 
# # Percentage of contigs with ARGs associated to ITR per ARG
# arg_no_is$ARO.Accession <- paste0("ARO:", gsub("_.*", "", gsub(".*_ARO:", "", arg_no_is$sseqid)))
# arg_no_is <- left_join(arg_no_is, card_metadata, by = "ARO.Accession")
# arg_no_is$Drug.Class.alt <- arg_no_is$Drug.Class
# arg_no_is$Drug.Class.alt[sapply(arg_no_is$Drug.Class.alt, function(x) str_count(x, ";")) > 1] <- "multidrug"
# arg_no_is$Drug.Class.alt[arg_no_is$Resistance.Mechanism == "antibiotic efflux"] <- paste(arg_no_is$Drug.Class.alt[arg_no_is$Resistance.Mechanism == "antibiotic efflux"], "efflux")
# arg_no_is$itr_cluster <- NA
# 
# arg_all <- rbind(arg_is, arg_no_is) %>%
#   left_join(metadata, by = "ID")
# 
# arg_with_is_dist <- arg_all %>% group_by(ID, Location, sample_type, ARO.Name) %>%
#   mutate(total_args = n_distinct(qseqid)) %>%
#   filter(!is.na(itr_cluster)) %>%
#   filter(Visit_Number == 1) %>%
#   group_by(ID, Location, sample_type, Location_sampletype, ARO.Name, Drug.Class.alt, Resistance.Mechanism, total_args) %>%
#   summarise(n = n_distinct(qseqid)) %>%
#   mutate(perc_args = n/total_args*100) %>%
#   filter(total_args >= 10) %>%
#   mutate(labels = paste0(ARO.Name, "\n(", Drug.Class.alt, ")"))
# 
# tiff("figures/percentage_args_with_is.tiff", width = 2500, height = 3500, res = 200)
# ggplot(arg_with_is_dist, aes(ARO.Name, perc_args, fill = sample_type)) +
#   geom_boxplot() +
#   geom_point() +
#   theme_classic() +
#   xlab("") + ylab("% ARGs that are associated to a potential insertion site") +
#   scale_fill_manual("GIT Site", values = cols) +
#   coord_flip() +
#   facet_grid(labels + Location_sampletype ~ ., space = "free", scale = "free", switch = "both") +
#   theme(strip.text.y = element_text(angle = 180, size = 6),
#         axis.title = element_text(size = 16),
#         axis.text = element_text(size = 12),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# dev.off()
# 

# #### Prevalence of ARG-carrying phages vs DDD
# antibiotic_use <- read.csv("db/resistanceMap_use_190221.csv", stringsAsFactors = FALSE)
# 
# arg_is_use <- trsps_arg_is %>%
#   left_join(metadata, by = "ID") %>%
#   group_by(Location, sample_type, Sample.name, Location_sampletype) %>%
#   ungroup() %>%
#   select(Location, sample_type, Location_sampletype, ID, ARO.Name, Drug.Class, Resistance.Mechanism) %>%
#   mutate(ARO.Name = c(strsplit(ARO.Name, ","))) %>%
#   unnest(ARO.Name) %>%
#   ungroup() %>%
#   mutate(Drug.Class = c(strsplit(Drug.Class, ";"))) %>%
#   unnest(Drug.Class) %>%
#   group_by(Location, sample_type, Drug.Class) %>%
#   summarise(n = n_distinct(ID)) %>%
#   inner_join(metadata_summary[metadata_summary$timepoint == 1,]) %>%
#   mutate(perc = n/n_total*100) %>%
#   inner_join(antibiotic_use, c("Location", "Drug.Class"="CARD.Class"))
# 
# # Add linear model
# arg_is_use_lm <- arg_is_use %>%
#   group_by(Location, sample_type) %>% 
#   do(mod = lm(perc ~ DDD.Per.1000.Pop, data = .))
# 
# # Plot
# class_cols <- c("black", brewer.pal(length(unique(arg_is_use$Drug.Class))-1, "Paired"))
# names(class_cols) <- unique(arg_is_use$Drug.Class)
# tiff("figures/arg_with_is_antibiotic_use.tiff", width = 4000, height = 700, res = 200)
# ggplot(arg_is_use, aes(DDD.Per.1000.Pop, perc, colour = Drug.Class)) +
#   geom_point(shape = 1) +
#   facet_grid(~Location + sample_type) +
#   theme_bw() +
#   scale_colour_manual("ARG Class", values = class_cols) +
#   ylab("% samples") + xlab("Defined Daily Dose Per 1000 Individuals")
# dev.off()
