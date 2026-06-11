```r
# Balanced Averaging Null Model (BANM)
#
# This script implements the BANM reference framework used to generate
# simulated microbial community coalescence profiles.
#
# The input should be a taxon-by-sample abundance table containing two donor
# communities. The first n_reps columns should correspond to Donor X and the
# next n_reps columns should correspond to Donor Y.
#
# The function computes bidirectional randomized donor mixtures:
#   1. Donor X focal replicates mixed with randomly sampled Donor Y replicates
#   2. Donor Y focal replicates mixed with randomly sampled Donor X replicates
#
# Sampling is conducted with replacement across N permutations.
# The final BANM profile is the mean of both directional mixtures.

compute_BANM <- function(otu_table_pair,
                         N = 1000,
                         n_reps = 6,
                         donor_1_prefix = "mixed_avg_X",
                         donor_2_prefix = "mixed_avg_Y",
                         final_prefix = "BANM_rep",
                         seed = NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  otu_table_pair <- as.data.frame(otu_table_pair)
  
  if (ncol(otu_table_pair) != n_reps * 2) {
    stop(
      "Input table must contain ",
      n_reps * 2,
      " columns: ",
      n_reps,
      " replicates from each donor."
    )
  }
  
  mixed_avg_1 <- matrix(0, nrow = nrow(otu_table_pair), ncol = n_reps)
  mixed_avg_2 <- matrix(0, nrow = nrow(otu_table_pair), ncol = n_reps)
  
  rownames(mixed_avg_1) <- rownames(otu_table_pair)
  rownames(mixed_avg_2) <- rownames(otu_table_pair)
  
  colnames(mixed_avg_1) <- paste0(donor_1_prefix, "_", seq_len(n_reps))
  colnames(mixed_avg_2) <- paste0(donor_2_prefix, "_", seq_len(n_reps))
  
  # Direction 1: Donor X -> Donor Y
  for (x_rep in seq_len(n_reps)) {
    
    permuted_sums_matrix <- matrix(
      0,
      nrow = nrow(otu_table_pair),
      ncol = N
    )
    
    for (perm in seq_len(N)) {
      
      sampled_y <- sample(
        (n_reps + 1):(n_reps * 2),
        n_reps,
        replace = TRUE
      )
      
      sums_matrix <- matrix(
        0,
        nrow = nrow(otu_table_pair),
        ncol = n_reps
      )
      
      for (i in seq_len(n_reps)) {
        y_rep <- sampled_y[i]
        sums_matrix[, i] <- otu_table_pair[, x_rep] + otu_table_pair[, y_rep]
      }
      
      permuted_sums_matrix[, perm] <- rowMeans(sums_matrix)
    }
    
    mixed_avg_1[, x_rep] <- rowMeans(permuted_sums_matrix)
  }
  
  # Direction 2: Donor Y -> Donor X
  for (y_rep in seq_len(n_reps)) {
    
    permuted_sums_matrix <- matrix(
      0,
      nrow = nrow(otu_table_pair),
      ncol = N
    )
    
    for (perm in seq_len(N)) {
      
      sampled_x <- sample(
        seq_len(n_reps),
        n_reps,
        replace = TRUE
      )
      
      sums_matrix <- matrix(
        0,
        nrow = nrow(otu_table_pair),
        ncol = n_reps
      )
      
      for (i in seq_len(n_reps)) {
        x_rep <- sampled_x[i]
        sums_matrix[, i] <- otu_table_pair[, n_reps + y_rep] + otu_table_pair[, x_rep]
      }
      
      permuted_sums_matrix[, perm] <- rowMeans(sums_matrix)
    }
    
    mixed_avg_2[, y_rep] <- rowMeans(permuted_sums_matrix)
  }
  
  final_mixed <- (mixed_avg_1 + mixed_avg_2) / 2
  
  final_mixed <- as.data.frame(final_mixed)
  colnames(final_mixed) <- paste0(final_prefix, "_", seq_len(n_reps))
  rownames(final_mixed) <- rownames(otu_table_pair)
  
  return(final_mixed)
}


convert_to_integer_counts <- function(df) {
  df[] <- lapply(df, function(x) as.integer(round(x)))
  return(df)
}


# Example usage:
#
# source("BANM/BANM_simulation.R")
#
# banm_AB_d1 <- compute_BANM(
#   otu_table_pair = simul_otu_ab_d1_t1,
#   N = 1000,
#   n_reps = 6,
#   donor_1_prefix = "mixed_avg_A",
#   donor_2_prefix = "mixed_avg_B",
#   final_prefix = "BANM_AB_d1",
#   seed = 123
# )
#
# banm_AB_d1 <- convert_to_integer_counts(banm_AB_d1)
```
