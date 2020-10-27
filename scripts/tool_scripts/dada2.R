# Support code for sequence-filtering using dada2 as in the analysis workflow described in
# Callahan et al. (2016), DADA2: High-resolution sample inference from Illumina amplicon data.
# https://doi.org/10.1038/nmeth.3869

library(dada2)
library(ShortRead)
library(ggplot2)


# ===== Sequence-filtering methods =====

# Parameters of filter methods:
#  maxN      After truncation, sequences with more than maxN Ns will be discarded. Note that dada currently does not allow Ns. (default: 0)
#  maxEE     After truncation, reads with higher than maxEE "expected errors" will be discarded. Expected errors are
#            calculated from the nominal definition of the quality score: EE = sum(10^(-Q/10)) (default: inf)
#  truncQ    Truncate reads at the first instance of a quality score less than or equal to truncQ. (default: 2)
#  truncLen  Truncate reads after truncLen bases. Reads shorter than this are discarded. (default: 0 = deactivated)
#  trimLeft  The number of nucleotides to remove from the start of each read. If both truncLen and trimLeft are provided,
#            filtered reads will have length truncLen-trimLeft. (default: 0)

filter_single <- function(forward_reads, max_n, max_ee, trunc_q, trunc_len, trim_left, plot_profile = FALSE) {

    cat("=== Filtering parameters ===", fill = TRUE)
    cat("forward reads:", forward_reads, fill = TRUE)
    cat("max_n:", max_n, fill = TRUE)
    cat("max_ee:", max_ee, fill = TRUE)
    cat("trunc_q:", trunc_q, fill = TRUE)
    cat("trunc_len:", trunc_len, fill = TRUE)
    cat("trim_left:", trim_left, fill = TRUE)

    if (plot_profile) {

        plotQualityProfile(forward_reads)
        ggsave(file = gsub(".fastq.*", ".pdf", forward_reads), device = "pdf")

    }

    # filter & trim
    filtered_forward <- gsub(".fastq", "_sf.fastq", forward_reads)
    fastqFilter(forward_reads, filtered_forward, maxN = max_n, maxEE = max_ee, truncQ = trunc_q, truncLen = trunc_len,
                trimLeft = trim_left, compress = TRUE, verbose = TRUE)

}

filter_paired <- function(forward_reads, reverse_reads, max_n, max_ee, trunc_q, trunc_len, trim_left, plot_profile = FALSE) {

    cat("=== Filtering parameters ===", fill = TRUE)
    cat("forward reads:", forward_reads, fill = TRUE)
    cat("reverse reads:", reverse_reads, fill = TRUE)
    cat("max_n:", max_n, fill = TRUE)
    cat("max_ee:", max_ee, fill = TRUE)
    cat("trunc_q:", trunc_q, fill = TRUE)
    cat("trunc_len:", trunc_len, fill = TRUE)
    cat("trim_left:", trim_left, fill = TRUE)

    if (plot_profile) {

        plotQualityProfile(forward_reads)
        ggsave(file = gsub(".fastq.*", ".pdf", forward_reads), device = "pdf")
        plotQualityProfile(reverse_reads)
        ggsave(file = gsub(".fastq.*", ".pdf", reverse_reads), device = "pdf")

    }

    # filter & trim
    filtered_forward <- gsub(".fastq", "_pf.fastq", forward_reads)
    filtered_reverse <- gsub(".fastq", "_pf.fastq", reverse_reads)
    fastqPairedFilter(c(forward_reads, reverse_reads), c(filtered_forward, filtered_reverse), maxN = max_n, maxEE = max_ee,
                      truncQ = trunc_q, truncLen = trunc_len, trimLeft = trim_left, compress = TRUE, verbose = TRUE)

}


# ===== Clustering methods =====

dada2_cluster_single <- function(forward_reads, otus_file) {

  derep_forward <- derepFastq(forward_reads, verbose = TRUE)

  dada_forward <- dada(derep_forward, err = inflateErr(tperr1, 3), selfConsist = TRUE, OMEGA_A = 1e-40)

  bim_forward <- isBimeraDenovo(dada_forward, allowOneOff = TRUE, verbose = TRUE)

  clustering <- dada_forward$clustering
  clustering_nc <- dada_forward$clustering[!bim_forward,]

  write.csv(clustering[,c("sequence", "abundance")], paste0(otus_file, "_dada2.txt"), row.names = TRUE)
  write.csv(clustering_nc[,c("sequence", "abundance")], paste0(otus_file, "_dada2_nc.txt"), row.names = TRUE)


  # OTU output
  deflines_fw <- readLines(forward_reads)
  deflines_fw <- deflines_fw[seq(1, length(deflines_fw), 4)]
  deflines_fw <- sapply(deflines_fw, function (x) strsplit(x, split = " ")[[1]][1], USE.NAMES = FALSE)
  deflines_fw <- sapply(deflines_fw, function (x) substr(x, 2, nchar(x)), USE.NAMES = FALSE)

  forward_partition <- dada_forward$map[derep_forward$map]
  partition_membership <- data.frame(sid = deflines_fw, forward = forward_partition, stringsAsFactors = F)

  get_otu <- function(pf) {
    member_ids <- partition_membership[which(partition_membership$forward == pf), "sid"]
    paste0(member_ids, "_1", collapse = " ")
  }


  rng <- as.numeric(row.names(clustering))#1:length(clustering$sequence)
  rng_nc <- as.numeric(row.names(clustering_nc))#1:length(clustering_nc$sequence)
  otus <- mapply(get_otu, rng, USE.NAMES = FALSE)
  otus_nc <- mapply(get_otu, rng_nc, USE.NAMES = FALSE)

  write(paste0(otus, collapse = "\n"), paste0(otus_file, "_dada2_otus.txt"))
  write(paste0(otus_nc, collapse = "\n"), paste0(otus_file, "_dada2_nc_otus.txt"))

}

dada2_cluster_paired <- function(forward_reads, reverse_reads, otus_file) {

  derep_forward <- derepFastq(forward_reads, verbose = TRUE)
  derep_reverse <- derepFastq(reverse_reads, verbose = TRUE)

  dada_forward <- dada(derep_forward, err = inflateErr(tperr1, 3), selfConsist = TRUE, OMEGA_A = 1e-40)
  dada_reverse <- dada(derep_reverse, err = inflateErr(tperr1, 3), selfConsist = TRUE, OMEGA_A = 1e-40)

  bim_forward <- isBimeraDenovo(dada_forward, allowOneOff = TRUE, verbose = TRUE)
  bim_reverse <- isBimeraDenovo(dada_reverse, allowOneOff = TRUE, verbose = TRUE)

  merger <- mergePairs(dada_forward, derep_forward, dada_reverse, derep_reverse, verbose = TRUE)
  merger_nc <- merger[!bim_forward[merger$forward] & !bim_reverse[merger$reverse],]

  write.csv(merger[,c("sequence", "abundance")], paste0(otus_file, "_dada2.txt"), row.names = TRUE)
  write.csv(merger_nc[,c("sequence", "abundance")], paste0(otus_file, "_dada2_nc.txt"), row.names = TRUE)


  # OTU output
  deflines_fw <- readLines(forward_reads)
  deflines_fw <- deflines_fw[seq(1, length(deflines_fw), 4)]
  deflines_fw <- sapply(deflines_fw, function (x) strsplit(x, split = " ")[[1]][1], USE.NAMES = FALSE)
  deflines_fw <- sapply(deflines_fw, function (x) substr(x, 2, nchar(x)), USE.NAMES = FALSE)

  forward_partition <- dada_forward$map[derep_forward$map]
  reverse_partition <- dada_reverse$map[derep_reverse$map]
  partition_membership <- data.frame(sid = deflines_fw, forward = forward_partition, reverse = reverse_partition, stringsAsFactors = F)


  get_otu <- function(pf, pr) {
    member_ids <- partition_membership[which(partition_membership$forward == pf & partition_membership$reverse == pr), "sid"]
    paste0(member_ids, "_1", collapse = " ")
  }

  otus <- mapply(get_otu, merger$forward, merger$reverse, USE.NAMES = FALSE)
  otus_nc <- mapply(get_otu, merger_nc$forward, merger_nc$reverse, USE.NAMES = FALSE)

  write(paste0(otus, collapse = "\n"), paste0(otus_file, "_dada2_otus.txt"))
  write(paste0(otus_nc, collapse = "\n"), paste0(otus_file, "_dada2_nc_otus.txt"))

}