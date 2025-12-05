#!/usr/bin/env Rscript
#
# Create phyloseq object from abundance and taxonomy tables
#

suppressPackageStartupMessages({
  library(phyloseq)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-a", "--abundance-table"), type="character",
              help="Abundance table TSV file (taxa as rows, samples as columns)"),
  make_option(c("-t", "--taxonomy-table"), type="character",
              help="Taxonomy table TSV file (taxa as rows, taxonomic ranks as columns)"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Sample metadata TSV file (samples as rows) [optional]"),
  make_option(c("-o", "--output"), type="character",
              help="Output RDS file (phyloseq object)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$`abundance-table`) || is.null(opt$`taxonomy-table`) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

cat("Loading data...\n")

# Load abundance table (OTU table)
cat(sprintf("  Reading: %s\n", opt$`abundance-table`))
otu_data <- read.table(opt$`abundance-table`, 
                       header=TRUE, 
                       sep="\t", 
                       row.names=1, 
                       check.names=FALSE,
                       quote="",
                       comment.char="")
otu_mat <- as.matrix(otu_data)
cat(sprintf("  Loaded %d taxa × %d samples\n", nrow(otu_mat), ncol(otu_mat)))

# Load taxonomy table
cat(sprintf("  Reading: %s\n", opt$`taxonomy-table`))
tax_data <- read.table(opt$`taxonomy-table`, 
                       header=TRUE, 
                       sep="\t", 
                       row.names=1, 
                       check.names=FALSE,
                       quote="",
                       comment.char="")
tax_mat <- as.matrix(tax_data)
cat(sprintf("  Loaded %d taxa × %d taxonomic ranks\n", nrow(tax_mat), ncol(tax_mat)))

# Verify matching taxa
if (!all(rownames(otu_mat) == rownames(tax_mat))) {
  warning("Taxa names don't match between abundance and taxonomy tables!")
  common_taxa <- intersect(rownames(otu_mat), rownames(tax_mat))
  cat(sprintf("  Using %d common taxa\n", length(common_taxa)))
  otu_mat <- otu_mat[common_taxa, ]
  tax_mat <- tax_mat[common_taxa, ]
}

# Create phyloseq components
OTU <- otu_table(otu_mat, taxa_are_rows=TRUE)
TAX <- tax_table(tax_mat)

# Load metadata if provided
if (!is.null(opt$metadata)) {
  cat(sprintf("  Reading metadata: %s\n", opt$metadata))
  sample_data_df <- read.table(opt$metadata, 
                                header=TRUE, 
                                sep="\t", 
                                row.names=1, 
                                check.names=FALSE,
                                quote="",
                                comment.char="")
  
  # Verify matching samples
  if (!all(colnames(otu_mat) %in% rownames(sample_data_df))) {
    warning("Not all samples have metadata!")
  }
  
  SAM <- sample_data(sample_data_df)
  ps <- phyloseq(OTU, TAX, SAM)
} else {
  ps <- phyloseq(OTU, TAX)
}

# Summary
cat("\nPhyloseq object created:\n")
cat(sprintf("  Taxa: %d\n", ntaxa(ps)))
cat(sprintf("  Samples: %d\n", nsamples(ps)))
cat(sprintf("  Taxonomic ranks: %s\n", paste(rank_names(ps), collapse=", ")))
if (!is.null(opt$metadata)) {
  cat(sprintf("  Sample variables: %s\n", paste(sample_variables(ps), collapse=", ")))
}

# Save phyloseq object
cat(sprintf("\nSaving: %s\n", opt$output))
saveRDS(ps, opt$output)

cat("✓ Done!\n")