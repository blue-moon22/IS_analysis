library(dplyr)

#### Functions ####
getBestBlastHit <- function(blast_result, qseqid_id){
  blast_result_tmp <- blast_result[blast_result$qseqid == qseqid_id,]
  remove_ids <- c()
  for(i in 1:(nrow(blast_result_tmp)-1)){
    
    for(j in (i+1):nrow(blast_result_tmp)){
      
      if(blast_result_tmp$qend[i] >= blast_result_tmp$qstart[j] & blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qstart[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qend[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[i] >= blast_result_tmp$qend[j] & blast_result_tmp$qstart[i] <= blast_result_tmp$qstart[j]){
        overlap_prop <- (blast_result_tmp$qend[i] - blast_result_tmp$qstart[i]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else if(blast_result_tmp$qend[j] >= blast_result_tmp$qend[i] & blast_result_tmp$qstart[j] <= blast_result_tmp$qstart[i]){
        overlap_prop <- (blast_result_tmp$qend[j] - blast_result_tmp$qstart[j]) / min(blast_result_tmp$length[i], blast_result_tmp$length[j])
      } else {
        overlap_prop <- 0
      }
      
      if(overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] != blast_result_tmp[j, "evalue"]){
        ind <- which.max(blast_result_tmp$evalue[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] & 
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] & 
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))
        
      } else if (overlap_prop > 0.2 & blast_result_tmp[i, "evalue"] == blast_result_tmp[j, "evalue"] & blast_result_tmp[i, "pident"] != blast_result_tmp[j, "pident"]){
        ind <- which.min(blast_result_tmp$pident[c(i,j)])
        remove_id <- which(blast_result$qseqid == blast_result_tmp$qseqid[c(i,j)[ind]] & 
                             blast_result$sseqid == blast_result_tmp$sseqid[c(i,j)[ind]] &
                             blast_result$qstart == blast_result_tmp$qstart[c(i,j)[ind]] & 
                             blast_result$qend == blast_result_tmp$qend[c(i,j)[ind]])
        remove_ids <- unique(c(remove_ids, remove_id))
      }
    }
  }
  return(remove_ids)
}

filterBySmallestEvalue <- function(alignment_results) {
  alignment_results_filtered <- alignment_results %>% 
    group_by(query_name) %>%
    filter(evalue == min(evalue)) %>%
    ungroup() %>%
    group_by(query_name) %>%
    filter(evalue_domain == min(evalue_domain))
  return(alignment_results_filtered)
}

#### Read samples ####
samples <- read.table("db/SAMPLES/samples.txt", stringsAsFactors = FALSE, header = FALSE)

#### Transposase ####
trsps_files <- list.files("data/TRANSPOSASE_IS", pattern = "*.out", full.names = TRUE)
trsps <- data.frame()
for (i in 1:length(trsps_files)) {
  info = file.info(trsps_files[i])
  if (info$size != 0) {
    tmp <- read.table(trsps_files[i], header = FALSE, stringsAsFactors = FALSE) %>%
      unique() %>%
      as.data.frame()
    tmp$ID <- gsub("_.*", "", gsub("data/TRANSPOSASE_IS/", "", trsps_files[i]))
    trsps <- rbind(trsps, tmp)
  }
}

names(trsps) <- c("query_name", "query_accession", "target_name", "target_accession", "evalue", 
                     "score", "bias", "evalue_domain", "score_domain", "bias_domain", "exp", "reg", 
                     "clu", "ov", "env", "dom", "rep", "inc", "ID")
trsps <- trsps[trsps$ID %in% samples$V1,]
trsps$contig <- sub("_[^_]+$", "", trsps$query_name)
trsps_filter <- filterBySmallestEvalue(trsps)

# Read transposase ITR file
trsps_itr_files <- list.files("data/TRANSPOSASE_ITRS", pattern = "*.tab", full.names = TRUE)
trsps_itr <- data.frame()
for (i in 1:length(trsps_itr_files)) {
  info = file.info(trsps_itr_files[i])
  if (info$size != 0) {
    tmp <- read.delim(trsps_itr_files[i], stringsAsFactors = FALSE)
    tmp$ID <- gsub("_.*", "", gsub("data/TRANSPOSASE_ITRS/", "", trsps_itr_files[i]))
    trsps_itr <- rbind(trsps_itr, tmp)
  }
}

trsps_is_filter <- inner_join(trsps_filter, trsps_itr, by = c("query_name"="transposase_seq", "contig", "ID"))

# Save data
saveRDS(trsps_is_filter, "data/filtered_trsps_is.RDS")

#### Insertion sites with ARGs ####
# Read blast out files
arg_is_files <- list.files("data/ARG_IS", pattern = "*.out", full.names = TRUE)
arg_is <- data.frame()
for (i in 1:length(arg_is_files)) {
  info = file.info(arg_is_files[i])
  if (info$size != 0) {
    tmp <- read.delim(arg_is_files[i], header = FALSE, stringsAsFactors = FALSE)
    tmp$ID <- gsub("_.*", "", gsub("data/ARG_IS/", "", arg_is_files[i]))
    arg_is <- rbind(arg_is, tmp)
  }
}

names(arg_is) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "is_cluster", "ID")
arg_is <- arg_is[arg_is$ID %in% samples$V1,]

# Filter
arg_is_filter_pident <- arg_is[arg_is$pident > 90,]

# Get top hits on exact regions
arg_is_filter <- arg_is_filter_pident %>%
  group_by(qseqid, qstart, qend) %>%
  filter(evalue == min(evalue)) %>%
  filter(bitscore == max(bitscore))

# Get best ARG hits
arg_is_filter_ident <- unique(arg_is_filter[arg_is_filter$qseqid %in% arg_is_filter$qseqid[duplicated(paste(arg_is_filter$qseqid, arg_is_filter$qstart, arg_is_filter$qend))],])
arg_is_filter_nonident <- arg_is_filter[!(arg_is_filter$qseqid %in% arg_is_filter_ident$qseqid),]
arg_is_filter_dup <- arg_is_filter_nonident[arg_is_filter_nonident$qseqid %in% arg_is_filter_nonident$qseqid[duplicated(arg_is_filter_nonident$qseqid)],]
arg_is_filter_sig <- arg_is_filter_nonident[!(arg_is_filter_nonident$qseqid %in% arg_is_filter_dup$qseqid),]
arg_is_remove_ids <- unlist(lapply(unique(arg_is_filter_dup$qseqid), function(x) getBestBlastHit(arg_is_filter_dup, x)))
arg_is_filter_dup <- arg_is_filter_dup[-arg_is_remove_ids,]
arg_is_filter <- rbind(arg_is_filter_ident, arg_is_filter_dup, arg_is_filter_sig)

# Add IS clusters back
arg_is_filter_hits <- arg_is_filter %>%
  select(-is_cluster) %>%
  inner_join(unique(arg_is_filter_pident[,c("qseqid", "is_cluster")]))

# Save data
saveRDS(unique(arg_is_filter_hits), "data/filtered_arg_is.RDS")

#### Contigs with ARGs ####
# Read blast out files
arg_contig_files <- list.files("data/ARG_CONTIGS", pattern = "*.out", full.names = TRUE)
arg_contig <- data.frame()
for (i in 1:length(arg_contig_files)) {
  info = file.info(arg_contig_files[i])
  if (info$size != 0) {
    tmp <- read.delim(arg_contig_files[i], header = FALSE, stringsAsFactors = FALSE)
    tmp$ID <- gsub("_.*", "", gsub("data/ARG_CONTIGS/", "", arg_contig_files[i]))
    arg_contig <- rbind(arg_contig, tmp)
  }
}

names(arg_contig) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "ID")
arg_contig <- arg_contig[arg_contig$ID %in% samples$V1,]

# Filter
arg_contig_filter <- arg_contig[arg_contig$pident > 90,]

# Remove contigs that have insertion sites
arg_contig_filter <- arg_contig_filter[!(paste(arg_contig_filter$qseqid, arg_contig_filter$ID) %in% paste(arg_is$qseqid, arg_is$ID)),]

# Get top hits on exact regions
arg_contig_filter <- arg_contig_filter %>%
  group_by(qseqid, qstart, qend) %>%
  filter(evalue == min(evalue)) %>%
  filter(bitscore == max(bitscore))

# Get best ARG hits
arg_filter_ident <- arg_contig_filter[arg_contig_filter$qseqid %in% arg_contig_filter$qseqid[duplicated(paste(arg_contig_filter$qseqid, arg_contig_filter$qstart, arg_contig_filter$qend))],]
arg_filter_nonident <- arg_contig_filter[!(arg_contig_filter$qseqid %in% arg_filter_ident$qseqid),]
arg_filter_dup <- arg_filter_nonident[arg_filter_nonident$qseqid %in% arg_filter_nonident$qseqid[duplicated(arg_filter_nonident$qseqid)],]
arg_filter_sig <- arg_filter_nonident[!(arg_filter_nonident$qseqid %in% arg_filter_dup$qseqid),]
remove_ids <- unlist(lapply(unique(arg_filter_dup$qseqid), function(x) getBestBlastHit(arg_filter_dup, x)))
arg_filter_dup <- arg_filter_dup[-remove_ids,]
arg_contig_filter <- rbind(arg_filter_ident, arg_filter_dup, arg_filter_sig)

# Save data
saveRDS(arg_contig_filter, "data/filtered_arg_contig.RDS")
