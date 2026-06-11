```r
# Create example phyloseq object
#
# Input files:
#   data/otu_table_example.csv
#   data/metadata_example.csv
#   data/tax_table_example.csv
#
# Output:
#   example_phyloseq_object.rds

library(phyloseq)

# Read input files
otu_table_example <- read.csv(
  "data/otu_table_example.csv",
  row.names = 1,
  check.names = FALSE
)

metadata_example <- read.csv(
  "data/metadata_example.csv",
  row.names = 1,
  check.names = FALSE
)

tax_table_example <- read.csv(
  "data/tax_table_example.csv",
  row.names = 1,
  check.names = FALSE
)

# Convert to matrices
otu_mat <- as.matrix(otu_table_example)
tax_mat <- as.matrix(tax_table_example)

# Make sure sample names match between OTU table and metadata
common_samples <- intersect(colnames(otu_mat), rownames(metadata_example))

otu_mat <- otu_mat[, common_samples, drop = FALSE]
metadata_example <- metadata_example[common_samples, , drop = FALSE]

# Make sure taxa names match between OTU table and taxonomy table
common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))

otu_mat <- otu_mat[common_taxa, , drop = FALSE]
tax_mat <- tax_mat[common_taxa, , drop = FALSE]

# Create phyloseq components
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
SAM <- sample_data(metadata_example)
TAX <- tax_table(tax_mat)

# Create phyloseq object
physeq_example <- phyloseq(OTU, SAM, TAX)

# Check object
physeq_example

# Save object
saveRDS(
  physeq_example,
  file = "data/example_phyloseq_object.rds"
)
```
