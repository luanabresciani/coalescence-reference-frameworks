```r
# Balanced Donor-Sum (BDS) simulation
# This script adds BDS-simulated microbial community profiles to a phyloseq object.
#
# Input:
#   1. A phyloseq object containing the original OTU table and taxonomy
#   2. A folder containing BDS simulation files in tab-delimited format
#
# Output:
#   A phyloseq object containing the original samples plus BDS-simulated samples

library(phyloseq)

# Read and align one BDS simulation file to the OTU table in the phyloseq object
make_bds_full <- function(physeq, file) {
  
  sim_df <- read.table(
    file,
    header = TRUE,
    sep = "\t",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(sim_df)[1] <- "OTU"
  rownames(sim_df) <- sim_df$OTU
  sim_df$OTU <- NULL
  
  # Extract simulation tag from filename, e.g. avg1, avg2, etc.
  fname <- basename(file)
  avg_tag <- sub("^.*(avg[0-9]+).*$", "\\1", fname)
  
  old_names <- colnames(sim_df)
  new_sample_names <- paste0(avg_tag, "_", old_names)
  colnames(sim_df) <- new_sample_names
  
  sim_mat <- as.matrix(sim_df)
  
  otu_orig <- as(otu_table(physeq), "matrix")
  
  if (!taxa_are_rows(otu_table(physeq))) {
    otu_orig <- t(otu_orig)
  }
  
  taxa_orig <- rownames(otu_orig)
  sim_samples <- colnames(sim_mat)
  
  sim_full <- matrix(
    0,
    nrow = length(taxa_orig),
    ncol = ncol(sim_mat),
    dimnames = list(taxa_orig, sim_samples)
  )
  
  common_otus <- intersect(rownames(sim_mat), taxa_orig)
  
  if (length(common_otus) == 0) {
    stop("No common OTUs between simulation file and phyloseq OTUs: ", fname)
  }
  
  sim_full[common_otus, ] <- sim_mat[common_otus, ]
  
  return(sim_full)
}


# Parse metadata from BDS simulated sample names
parse_bds_metadata <- function(sample_names) {
  
  parts <- strsplit(sample_names, "_")
  
  avg_tag <- vapply(parts, `[`, "", 1)
  pair    <- vapply(parts, `[`, "", 3)
  div_raw <- vapply(parts, `[`, "", 4)
  rep     <- vapply(parts, `[`, "", 5)
  
  avg_num <- as.numeric(sub("^avg", "", avg_tag))
  div <- as.numeric(sub("^d", "", div_raw))
  
  data.frame(
    ID        = sample_names,
    Soil      = paste0("BDS_avg", avg_num, "_", pair),
    Diversity = div,
    Time      = "T0",
    Rep       = rep,
    row.names = sample_names,
    stringsAsFactors = FALSE
  )
}


# Add multiple BDS simulation files to a phyloseq object
add_bds_to_physeq <- function(physeq, files) {
  
  sim_list <- lapply(files, make_bds_full, physeq = physeq)
  sim_all <- do.call(cbind, sim_list)
  
  otu_orig <- as(otu_table(physeq), "matrix")
  taxa_rows <- taxa_are_rows(otu_table(physeq))
  
  if (!taxa_rows) {
    otu_orig <- t(otu_orig)
  }
  
  otu_new_mat <- cbind(otu_orig, sim_all)
  
  if (taxa_rows) {
    otu_new <- otu_table(otu_new_mat, taxa_are_rows = TRUE)
  } else {
    otu_new <- otu_table(t(otu_new_mat), taxa_are_rows = FALSE)
  }
  
  sam_orig <- as.data.frame(sample_data(physeq))
  orig_cols <- colnames(sam_orig)
  
  new_sample_names <- colnames(sim_all)
  sam_meta <- parse_bds_metadata(new_sample_names)
  
  for (nm in setdiff(orig_cols, colnames(sam_meta))) {
    sam_meta[[nm]] <- NA
  }
  
  sam_meta <- sam_meta[, orig_cols]
  
  sam_combined <- rbind(sam_orig, sam_meta)
  sam_new <- sample_data(sam_combined)
  
  physeq_new <- phyloseq(
    otu_new,
    tax_table(physeq),
    sam_new
  )
  
  return(physeq_new)
}


# Example usage:
#
# load("path/to/phyloseq_object.Rdata")
#
# bds_path <- "path/to/BDS_simulation_files/"
#
# bds_files <- list.files(
#   bds_path,
#   pattern = "\\.txt$",
#   full.names = TRUE
# )
#
# physeq_with_bds <- add_bds_to_physeq(
#   physeq = physeq_original,
#   files = bds_files
# )
#
# save(
#   physeq_with_bds,
#   file = "physeq_with_BDS_simulations.Rdata"
# )
```
