# ============================================================
# Purpose: Process paired-end 16S V4 amplicon reads with DADA2,
# generate a non-chimeric ASV table, assign taxonomy with SILVA,
# and build a phyloseq object for downstream analyses.
# ============================================================

# =========================
# 1. Load libraries
# =========================
library(dada2)
library(ggplot2)
library(phyloseq)
library(dplyr)
library(reshape2)
library(tibble)
library(RColorBrewer)
library(randomcoloR)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)

# =========================
# 2. Define input paths
# =========================
# Update these paths as needed for your system/project structure.
# Recommended structure:
#   data/raw/       -> paired-end FASTQ files
#   data/metadata/  -> sample metadata
#   ref/            -> SILVA reference database
#   output/         -> generated tables/files

path <- "data/raw"
metadata_file <- "data/metadata/metadata.txt"
silva_ref <- "ref/silva_nr99_v138.2_toGenus_trainset.fa.gz"
output_dir <- "output"
filt_path <- file.path(output_dir, "filtered")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(filt_path, recursive = TRUE, showWarnings = FALSE)

# Confirm that the FASTQ files are visible
list.files(path)

# =========================
# 3. Identify input FASTQ files
# =========================
fnFs <- sort(list.files(path, pattern = "_R1.fastq.gz$", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq.gz$", full.names = TRUE))

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No FASTQ files found. Check the 'path' value and file naming pattern.")
}
if (length(fnFs) != length(fnRs)) {
  stop("Forward and reverse read files are not matched in number.")
}

# Sample name is everything before the first underscore
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# =========================
# 4. Create filtered file names
# =========================
filtFs <- file.path(filt_path, paste0(sample.names, "_1_filt.fq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_2_filt.fq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# =========================
# 5. Quality filtering and trimming
# =========================
# Filtering choices:
#   maxN = 0      -> discard reads with ambiguous bases
#   truncQ = 10   -> remove low-quality tails while retaining reads
#   trimLeft = 24 and trimRight = 24 -> remove primer/adaptor-associated bases
#   minLen = 50   -> discard very short reads after trimming
# On Windows, multithread should generally be FALSE.
out <- filterAndTrim(
  fnFs, filtFs, fnRs, filtRs,
  maxN = 0,
  truncQ = 10,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = FALSE,
  trimLeft = 24,
  trimRight = 24,
  minLen = 50
)

write.table(out, file = file.path(output_dir, "filtered_information.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# =========================
# 6. Learn error rates
# =========================
# Error models are estimated separately for forward and reverse reads.
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

# =========================
# 7. Dereplication
# =========================
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

# =========================
# 8. Denoising / ASV inference
# =========================
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# =========================
# 9. Merge paired reads
# =========================
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# =========================
# 10. Construct sequence table
# =========================
seqtab <- makeSequenceTable(mergers)
print(dim(seqtab))

# Inspect distribution of sequence lengths
ASVseqlengthtab <- table(nchar(getSequences(seqtab)))
write.table(ASVseqlengthtab,
            file = file.path(output_dir, "ASV_sequence_length_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# =========================
# 11. Remove chimeras
# =========================
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = FALSE,
  verbose = TRUE
)

print(dim(seqtab.nochim))
print(sum(seqtab.nochim) / sum(seqtab))

write.table(seqtab.nochim,
            file = file.path(output_dir, "seqtab_nochim.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# =========================
# 12. Track reads through the pipeline
# =========================
getN <- function(x) sum(getUniques(x))

track <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(mergers, getN),
  rowSums(seqtab),
  rowSums(seqtab.nochim)
)

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names

write.table(track,
            file = file.path(output_dir, "dada_read_stats.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# Save the non-chimeric ASV table as an R object for downstream use
saveRDS(seqtab.nochim, file = file.path(output_dir, "seqtab.nochim.rds"))

# Optional reload step
seqtab.nochim <- readRDS(file.path(output_dir, "seqtab.nochim.rds"))

# Transpose ASV table for export if needed
seqtab.nochim_t <- t(seqtab.nochim)
write.table(seqtab.nochim_t,
            file = file.path(output_dir, "silva_ASVs_table_t.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# =========================
# 13. Assign taxonomy
# =========================
# Taxonomic assignment is performed against SILVA v138.2.
# Resolution is expected to be strongest at higher ranks and often limited at genus/species level for the V4 region.
taxa <- assignTaxonomy(seqtab.nochim, silva_ref, multithread = FALSE)

# Fill missing ranks conservatively using the next highest assigned rank
taxon <- as.data.frame(taxa, stringsAsFactors = FALSE)
taxon$Phylum[is.na(taxon$Phylum)] <- taxon$Kingdom[is.na(taxon$Phylum)]
taxon$Class[is.na(taxon$Class)] <- taxon$Phylum[is.na(taxon$Class)]
taxon$Order[is.na(taxon$Order)] <- taxon$Class[is.na(taxon$Order)]
taxon$Family[is.na(taxon$Family)] <- taxon$Order[is.na(taxon$Family)]
taxon$Genus[is.na(taxon$Genus)] <- taxon$Family[is.na(taxon$Genus)]

write.table(taxon,
            file = file.path(output_dir, "silva_taxa_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(seqtab.nochim,
            file = file.path(output_dir, "silva_ASVs_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# =========================
# 14. Create phyloseq object
# =========================
otu <- read.table(
  file.path(output_dir, "silva_ASVs_table_t.txt"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)
otu <- t(otu)

taxon <- read.table(
  file.path(output_dir, "silva_taxa_table.txt"),
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

samples <- read.table(
  metadata_file,
  sep = "\t", header = TRUE, row.names = 1, check.names = FALSE
)

OTU <- otu_table(otu, taxa_are_rows = FALSE)
TAX <- tax_table(as.matrix(taxon))
sampledata <- sample_data(samples)

ps <- phyloseq(OTU, sampledata, TAX)
print(ps)

# Remove non-target sequences
ps <- subset_taxa(ps, Family != "Mitochondria")
ps <- subset_taxa(ps, Order != "Chloroplast")
ps <- subset_taxa(ps, Kingdom != "Eukaryota")
ps <- subset_taxa(ps, !is.na(Kingdom))
print(ps)

# Export cleaned OTU and taxonomy tables
otu_clean <- as(otu_table(ps), "matrix")
taxon_clean <- as(tax_table(ps), "matrix")

write.table(otu_clean,
            file = file.path(output_dir, "silva_nochloronomito_otu_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(taxon_clean,
            file = file.path(output_dir, "silva_nochloronomito_taxa_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# Summarize retained taxonomic diversity
length(get_taxa_unique(ps, "Phylum"))
length(get_taxa_unique(ps, "Class"))
length(get_taxa_unique(ps, "Order"))
length(get_taxa_unique(ps, "Family"))
length(get_taxa_unique(ps, "Genus"))

# =========================
# 15. Filter low-abundance ASVs
# =========================
# Retain ASVs with mean abundance > 5 reads across samples.
psmin5 <- filter_taxa(ps, function(x) mean(x) > 5, TRUE)

ntaxa(psmin5)
print(psmin5)

length(get_taxa_unique(psmin5, "Phylum"))
length(get_taxa_unique(psmin5, "Class"))
length(get_taxa_unique(psmin5, "Order"))
length(get_taxa_unique(psmin5, "Family"))
length(get_taxa_unique(psmin5, "Genus"))

otu_min5 <- as(otu_table(psmin5), "matrix")
taxon_min5 <- as(tax_table(psmin5), "matrix")

write.table(otu_min5,
            file = file.path(output_dir, "psmin5_ASV_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(taxon_min5,
            file = file.path(output_dir, "psmin5_taxa_table.txt"),
            sep = "\t", quote = FALSE, col.names = NA)
