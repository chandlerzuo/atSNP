library(MotifAnalysis)
library(Rcpp)

motif_file <- "/p/keles/ENCODE-CHARGE/volume1/ENCODE-Motifs/encode_motifs_for_fimo.txt"

motif_library <- LoadMotifLibrary(motif_file)

snpInfo <- LoadSNPData("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed")

if(FALSE) {
  snpInfo$sequence_matrix <- snpInfo$sequence_matrix[, 1:100]
  snpInfo$a1 <- snpInfo$a1[1:100]
  snpInfo$a2 <- snpInfo$a2[1:100]
}

system.time(results <- ComputeMotifScore(motif_library, snpInfo, ncores = 20))

save(results, file = "/p/keles/ENCODE-CHARGE/volume2/SNP/my_motif_analysis.Rda")

table(results$on_odds[results$log_lik_ratio < 0] > 0)

table(results$off_odds[results$log_lik_ratio > 0] > 0)

## compare the results with fimo

## check for matches on a1

match_a1 <- results$match_a1
colnames(match_a1) <- motif_library$motif
rownames(match_a1) <- colnames(snpInfo$sequence_matrix)

fimo_a1 <- read.table("/p/keles/ENCODE-CHARGE/volume2/SNP/fasta/snp_ref_window_30_a1.fasta_fimo_filtered.txt", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "")
fimo_a1 <- fimo_a1[fimo_a1$snpid %in% colnames(snpInfo$sequence_matrix) & fimo_a1$pattern_name %in% motif_library$motif, ]
ids <- which(fimo_a1$strand == "-")
match_pos <- - fimo_a1$SNP + fimo_a1$start + 30
match_rev <- fimo_a1$SNP + 30 - fimo_a1$stop
match_pos[ids] <- - match_rev[ids]

my_match_pos <- match_a1[cbind(fimo_a1$snpid, fimo_a1$pattern_name)]

table(abs(my_match_pos - match_pos) > 0)

## check for matches on a2

match_a2 <- results$match_a2
colnames(match_a2) <- motif_library$motif
rownames(match_a2) <- colnames(snpInfo$sequence_matrix)

fimo_a2 <- read.table("/p/keles/ENCODE-CHARGE/volume2/SNP/fasta/snp_ref_window_30_a2.fasta_fimo_filtered.txt", header = TRUE,
                     stringsAsFactors = FALSE, comment.char = "")
fimo_a2 <- fimo_a2[fimo_a2$snpid %in% colnames(snpInfo$sequence_matrix) & fimo_a2$pattern_name %in% motif_library$motif, ]
ids <- which(fimo_a2$strand == "-")
match_pos <- - fimo_a2$SNP + fimo_a2$start + 30
match_rev <- fimo_a2$SNP + 30 - fimo_a2$stop
match_pos[ids] <- - match_rev[ids]

my_match_pos <- match_a2[cbind(fimo_a2$snpid, fimo_a2$pattern_name)]

table(abs(my_match_pos - match_pos) > 0)
