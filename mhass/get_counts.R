#!/usr/bin/env Rscript
# get_counts.R
# Purpose: Simulate Titan ASV counts using metaSPARSim with genome-aware distribution logic

suppressMessages({
  # Check for required packages
  packages_needed <- c("optparse", "BiocManager")
  packages_to_install <- packages_needed[!sapply(packages_needed, requireNamespace, quietly = TRUE)]
  
  if (length(packages_to_install) > 0) {
    message("Installing required packages: ", paste(packages_to_install, collapse = ", "))
    install.packages(packages_to_install, repos = "https://cloud.r-project.org/", quiet = TRUE)
  }
  
  # Check for Biostrings
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    message("Installing Biostrings from Bioconductor")
    BiocManager::install("Biostrings", quiet = TRUE)
  }
  
  # Check for metaSPARSim
  if (!requireNamespace("metaSPARSim", quietly = TRUE)) {
    message("Installing metaSPARSim from GitLab")
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes", repos = "https://cloud.r-project.org/", quiet = TRUE)
    }
    
    # Try to install from GitLab
    tryCatch({
      remotes::install_gitlab("sysbiobig/metaSPARSim", quiet = TRUE)
    }, error = function(e) {
      # Try alternative URL if the first fails
      tryCatch({
        remotes::install_git("https://gitlab.com/sysbiobig/metasparsim.git", quiet = TRUE)
      }, error = function(e2) {
        message("Failed to install metaSPARSim. Please install manually and try again.")
        message("Error details: ", e2$message)
        quit(status = 1)
      })
    })
  }
  
  # Load required libraries
  library(metaSPARSim)
  library(Biostrings)
  library(optparse)
})

# CLI options
option_list <- list(
  make_option(c("-f", "--fasta"), type="character", help="Input FASTA file of ASVs"),
  make_option(c("-n", "--samples"), type="integer", default=10),
  make_option(c("-r", "--reads"), type="integer", default=10000),
  make_option(c("-d", "--dispersion"), type="double", default=0.1),
  make_option(c("-G", "--genome-map"), type="character", help="TSV: asv_id<TAB>genome_id"),
  make_option(c("--genome-distribution"), type="character", default="uniform",
              help="One of: uniform, lognormal, powerlaw, or empirical:<file.tsv>"),
  make_option(c("-o", "--out"), type="character", default="sim_counts.tsv"),
  make_option(c("-m", "--meta"), type="character", default="sim_meta.tsv")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Validate required parameters
if (is.null(opt$fasta)) {
  stop("Required parameter --fasta is missing")
}
if (is.null(opt$`genome-map`)) {
  stop("Required parameter --genome-map is missing")
}

# Set seed for reproducibility
set.seed(42)

# Load ASVs
message("Loading ASV sequences...")
asv_seqs <- readDNAStringSet(opt$fasta)
asv_ids <- names(asv_seqs)
n_samples <- opt$samples
lib_sizes <- rep(opt$reads, n_samples)
message(sprintf("=> Loaded %d ASVs from %s <=", length(asv_ids), opt$fasta))

# Load genome map
message("Loading genome mapping...")
genome_map <- read.delim(opt$`genome-map`, header=TRUE, stringsAsFactors=FALSE)
stopifnot(all(c("asvid", "genomeid") %in% colnames(genome_map)))
genome_map <- genome_map[genome_map$asvid %in% asv_ids, ]
asv_to_genome <- setNames(genome_map$genomeid, genome_map$asvid)
genome_ids <- unique(genome_map$genomeid)
n_genomes <- length(genome_ids)
message(sprintf("=> Found %d genomes in mapping <=", n_genomes))

# Compute genome proportions
message("Computing genome proportions using distribution: ", opt$`genome-distribution`)
dist_spec <- opt$`genome-distribution`
if (startsWith(dist_spec, "empirical:")) {
  file_path <- sub("empirical:", "", dist_spec)
  if (!file.exists(file_path)) {
    stop("Empirical abundance file not found: ", file_path)
  }
  gdf <- read.delim(file_path, header=FALSE)
  colnames(gdf) <- c("genomeid", "proportion")
  genome_props <- setNames(gdf$proportion, gdf$genomeid)
  genome_props <- genome_props[names(genome_props) %in% genome_ids]
  if (length(genome_props) < n_genomes) {
    warning("Not all genomes found in empirical abundance file. Missing genomes will have zero abundance.")
  }
  genome_props <- genome_props / sum(genome_props)
} else if (dist_spec == "uniform") {
  genome_props <- setNames(rep(1 / n_genomes, n_genomes), genome_ids)
} else if (dist_spec == "lognormal") {
  vals <- rlnorm(n_genomes, meanlog=0, sdlog=1)
  genome_props <- setNames(vals / sum(vals), genome_ids)
} else if (dist_spec == "powerlaw") {
  alpha <- 1.5
  vals <- (1:n_genomes)^(-alpha)
  genome_props <- setNames(vals / sum(vals), genome_ids)
} else {
  stop("Invalid --genome-distribution. Use uniform, lognormal, powerlaw, or empirical:<file>")
}

# Build ASV-level mu
message("Building ASV-level intensity parameters...")
asv_counts_per_genome <- table(genome_map$genomeid)
mu_vec <- sapply(asv_to_genome, function(gid) {
  prop <- genome_props[as.character(gid)]
  n_asvs <- asv_counts_per_genome[as.character(gid)]
  if (is.na(prop) || is.na(n_asvs)) return(0)
  prop / n_asvs
})
mu_vec <- mu_vec / sum(mu_vec)  # Normalize again just in case

# Final prep for metaSPARSim
feature_names <- paste0("Feature_", seq_along(mu_vec))
names(mu_vec) <- feature_names
phi_vec <- rep(opt$dispersion, length(mu_vec))
names(phi_vec) <- feature_names

dataset_parameters <- list(G1 = list(
  intensity = mu_vec,
  variability = phi_vec,
  lib_size = lib_sizes
))

message("=> Simulating counts via metaSPARSim...")
result <- metaSPARSim(dataset_parameters)

counts <- result$counts
rownames(counts) <- names(asv_to_genome)
colnames(counts) <- paste0("Sample", seq_len(ncol(counts)))

# Output
counts_df <- data.frame(ASVID = rownames(counts), counts, check.names = FALSE)
write.table(counts_df, file = opt$out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(data.frame(Sample=colnames(counts), LibrarySize=colSums(counts)),
          file = opt$meta, row.names = FALSE, sep = "\t", quote = FALSE)
message("=> Simulation complete. <=")
message(paste("- Count matrix:", opt$out))
message(paste("- Metadata:", opt$meta))
