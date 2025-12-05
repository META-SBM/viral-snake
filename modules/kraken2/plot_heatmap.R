#!/usr/bin/env Rscript
#
# Create filtered abundance heatmap from phyloseq object
# Filters taxa by prevalence and creates log-transformed heatmap
#

suppressPackageStartupMessages({
  library(phyloseq)
  library(microbiome)
  library(ComplexHeatmap)
  library(optparse)
  library(circlize)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-p", "--phyloseq"), type="character",
              help="Phyloseq RDS file"),  # CHANGED: single input now
  make_option(c("-o", "--output"), type="character",
              help="Output PDF file"),
  make_option(c("--prevalence"), type="numeric", default=0.1,
              help="Prevalence threshold (0-1) [default: 0.1 = 10%% of samples]"),
  make_option(c("-d", "--detection"), type="numeric", default=0,
              help="Detection threshold (min reads) [default: 0]"),
  make_option(c("--domain"), type="character", default=NULL,
              help="Filter by domain (e.g., 'Viruses', 'Bacteria', 'Archaea') [default: all domains]")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$phyloseq) || is.null(opt$output)) {  # CHANGED
  print_help(opt_parser)
  stop("Missing required arguments", call.=FALSE)
}

cat("Loading phyloseq object...\n")

# Load phyloseq object
ps <- readRDS(opt$phyloseq)  # CHANGED: just load it

cat(sprintf("Initial phyloseq: %d taxa, %d samples\n", 
            ntaxa(ps), nsamples(ps)))

# Filter by prevalence using microbiome package
cat(sprintf("Filtering: prevalence >= %.1f%%, detection >= %d\n", 
            opt$prevalence * 100, opt$detection))

ps_filtered <- core(ps, 
                    detection = opt$detection, 
                    prevalence = opt$prevalence)

cat(sprintf("Filtered phyloseq: %d taxa, %d samples\n", 
            ntaxa(ps_filtered), nsamples(ps_filtered)))

# Filter by domain if specified
if (!is.null(opt$domain)) {
  cat(sprintf("Filtering by domain: %s\n", opt$domain))
  
  # Get taxonomy table
  tax_table_temp <- as.data.frame(tax_table(ps_filtered))
  
  # Find taxa matching the domain (case-insensitive)
  domain_match <- grepl(opt$domain, tax_table_temp$Domain, ignore.case = TRUE)
  
  if (sum(domain_match) == 0) {
    stop(sprintf("No taxa found for domain '%s'! Available domains: %s", 
                 opt$domain, 
                 paste(unique(tax_table_temp$Domain), collapse=", ")), 
         call.=FALSE)
  }
  
  # Subset phyloseq object
  taxa_to_keep <- rownames(tax_table_temp)[domain_match]
  ps_filtered <- prune_taxa(taxa_to_keep, ps_filtered)
  
  cat(sprintf("After domain filter: %d taxa\n", ntaxa(ps_filtered)))
}

if (ntaxa(ps_filtered) == 0) {
  stop("No taxa passed filtering! Try reducing prevalence threshold.", call.=FALSE)
}

# Extract abundance matrix and log transform
abundance_mat <- as.matrix(otu_table(ps_filtered))
log_abundance <- log10(abundance_mat + 1)

# Get taxonomy table for filtered taxa
tax_filtered <- as.data.frame(tax_table(ps_filtered))

# Create row labels: tax_id + Species name
tax_ids <- rownames(log_abundance)
species_names <- tax_filtered$Species
row_labels <- paste(tax_ids, species_names)

# Get Domain information for annotation (handle NAs)
domains <- tax_filtered$Domain
domains[is.na(domains)] <- "Unknown"  # Replace NA with "Unknown"
domains <- as.character(domains)      # Ensure character vector
domains <- trimws(domains)            # Remove whitespace

# Create color palette for domains
unique_domains <- unique(domains)
domain_colors <- setNames(
  rainbow(length(unique_domains), s=0.7, v=0.8),
  unique_domains
)

# Verify all domains have colors
if (!all(domains %in% names(domain_colors))) {
  stop("Some domains don't have corresponding colors!", call.=FALSE)
}

cat("Creating heatmap...\n")

# Calculate dimensions based on fixed cell size
cell_width <- 4    # mm per cell (smaller)
cell_height <- 3   # mm per cell (smaller)
n_samples <- ncol(log_abundance)
n_taxa <- nrow(log_abundance)

# Calculate PDF dimensions
pdf_width <- (n_samples * cell_width / 25.4) + 5  # add 5" for labels/legend
pdf_height <- (n_taxa * cell_height / 25.4) + 3   # add 3" for title/legend

cat(sprintf("PDF dimensions: %.1f x %.1f inches\n", pdf_width, pdf_height))

# Create row annotation for Domain
row_ha <- rowAnnotation(
  Domain = domains,
  col = list(Domain = domain_colors),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Domain = list(title = "Domain", 
                  title_gp = gpar(fontsize = 10, fontface = "bold"),
                  labels_gp = gpar(fontsize = 9))
  ),
  width = unit(5, "mm")
)

# Create heatmap with fixed cell size
pdf(opt$output, width=pdf_width, height=pdf_height)

ht <- Heatmap(
  log_abundance,
  name = "log10(counts+1)",
  
  # Clustering with gaps between clusters
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_km = 4,
  column_km = 7,
  
  # Fixed cell dimensions
  width = unit(n_samples * cell_width, "mm"),
  height = unit(n_taxa * cell_height, "mm"),
  
  # Row labels (taxa with species names)
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 7),
  row_names_max_width = unit(12, "cm"),
  
  # Column labels (sample names)
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  
  # Color scheme
  col = colorRamp2(
    c(0, max(log_abundance)/2, max(log_abundance)),
    c("blue", "white", "red")
  ),
  
  # Row annotation (Domain)
  right_annotation = row_ha,
  
  # Options
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = sprintf("Taxa abundance (n=%d taxa, %d samples)", 
                        n_taxa, n_samples),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  
  # Show dendrograms
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  
  row_dend_width = unit(3, "cm"),
  column_dend_height = unit(3, "cm"),
  # Legend parameters
  heatmap_legend_param = list(
    title = "Abundance",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    legend_height = unit(4, "cm")
  )
)

draw(ht)
dev.off()

cat(sprintf("✓ Heatmap saved: %s\n", opt$output))