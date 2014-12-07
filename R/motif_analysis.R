#' @name LoadMotifLibrary
#' @title Load position weight matrices.
#' @description Load the file for position weight matrices for motifs.
#' @param filename A file containing MEME format: http://meme.nbcr.net/meme/doc/meme-format.html.
#' @param tag A string that marks the description line of the position weight matrix.
#' @param skiprows Number of description lines before each position weight matrix.
#' @param skipcols Number of columns to be skipped in the position weight matrix. 
#' @param transpose If TRUE (default), then the position weight matrix should have 4 columns. Otherwise, it should have 4 rows.
#' @param field The index of the field in the description line, seperated by space, that indicates the motif name.
#' @param sep The string seperator to separate each lines of the matrix. Default: " ".
#' @details This function reads the formatted file containing motif information and convert them into a list of position weight matrices. The list of arguments should provide enough flexibility of importing a varying number of formats. Som eexamples are the following:
#' For MEME format, the suggested arguments are: tag = 'Motif', skiprows = 2, skipcols = 0, transpose = FALSE, field = 2, sep = " ";
#' For motif files from JOHNSON lab (i.e. http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt), the suggested arguments are: tag = '/NAME', skiprows = 1, skipcols = 0, transpose = FALSE, field = 2, sep = "\t";
#' For JASPAR pfm matrices (i.e. http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt), the suggested arguments are: tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1, sep = "\t";
#' For the TRANSFAC library provided by UCF bioinformatics groups (i.e. http://gibbs.biomed.ucf.edu/PreDREM/download/nonredundantmotif.transfac), the suggested arguments are: tag = "DE", skiprows = 1, skipcols = 1, transpose = FALSE, field = 2, sep = "\t".
#' @return A list object of position weight matrices.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' pwms <- LoadMotifLibrary("/p/keles/ENCODE-CHARGE/volume1/ENCODE-Motifs/encode_motifs_for_fimo.txt")
#' pwms <- LoadMotifLibrary("http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt", tag = "/NAME", skiprows = 1, skipcols = 0, transpose = FALSE, field = 2, sep = "\t")
#' pwms <- LoadMotifLibrary("http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt", tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1, sep = "\t")
#' pwms <- LoadMotifLibrary("http://gibbs.biomed.ucf.edu/PreDREM/download/nonredundantmotif.transfac", tag = "DE", skiprows = 1, skipcols = 1, transpose = FALSE, field = 2, sep = "\t")
#' }
#' @useDynLib atSNP
#' @export
LoadMotifLibrary <- function(filename, tag = "MOTIF", transpose = FALSE, field = 2, sep = " ", skipcols = 0, skiprows = 2) {
  lines <- readLines(filename)
  motifLineNums <- grep(tag, lines)
  if(length(strsplit(lines[motifLineNums[1]], " ")[[1]]) > 1) {
    motifnames <-
      sapply(strsplit(lines[motifLineNums], " "), function(x) x[field])
  } else {
    motifnames <-
      sapply(strsplit(lines[motifLineNums], "\t"), function(x) x[field])
  }
  allmats <- as.list(seq_along(motifnames))
  
  for(matrixId in seq_along(motifLineNums)) {
    motifLineNum <- motifLineNums[matrixId] + skiprows
    if(!transpose) {
      pwm <- NULL
      nrows <- 0
      while(length(tmp <- strsplit(lines[nrows + motifLineNum], split = sep)[[1]]) >= 4 + skipcols) {
        tmp <- as.numeric(tmp[skipcols + seq(4)])
        if(sum(is.na(tmp)) == 0) {
          pwm <- rbind(pwm, tmp)
          nrows <- nrows + 1
        }
      }
    } else {
      nrows <- 4
      pwm <-
        matrix(as.numeric(unlist(strsplit(lines[seq(nrows) + motifLineNum - 1], split = sep))), ncol = 4)
      if(skipcols > 0) {
        pwm <- pwm[-seq(skipcols), ]
      }
    }
    pwm <- pwm / apply(pwm, 1, sum)
    pwm <- t(apply(pwm, 1,
                   function(x) {
                     x[x < 1e-10] <- 1e-10 / (1 - 1e-10 * sum(x < 1e-10)) * sum(x)
                     return(x / sum(x))
                   }))
    rownames(pwm) <- NULL
    allmats[[matrixId]] <- pwm
  }
  names(allmats) <- motifnames
  return(allmats)
}


#' @name LoadSNPData
#' @title Load the SNP information and code the genome sequences around the SNP locations.
#' @description Load the SNP data.
#' @param filename A table containing the SNP information. Must contain at least four columns:
#' \tabular{ll}{
#' chr \tab chromosome.\cr
#' snp \tab The nucleotide position of the SNP.\cr
#' snpid \tab The names of the SNPs.\cr
#' a1 \tab The deoxyribose at the reference genome.\cr
#' a2 \tab The deoxyribose at the SNP genome.\cr
#' }
#' @param genome.lib A string of the library name for the genome version. Default: "BSgenome.Hsapiens.UCSC.hg19".
#' @param half.window.size An integer for the half window size around the SNP within which the motifs are matched. Default: 30.
#' @param ... Other parameters passed to 'read.table'.
#' @details TODO.
#' @return A list object containing the following components:
#' \tabular{ll}{
#' sequence_matrix \tab A list of integer vectors representing the deroxyribose sequence around each SNP.\cr
#' a1 \tab An integer vector for the deroxyribose at the SNP location on the reference genome.\cr
#' a2 \tab An integer vector for the deroxyribose at the SNP location on the SNP genome.\cr
#' }
#' The results are coded as: "A"-1, "C"-2, "G"-3, "T"-4.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{LoadSNPData("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed")}
#' @useDynLib atSNP
#' @export
LoadSNPData <- function(filename, genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                        half.window.size = 30, ...) {
  ## load the corresponding genome version
  library(package = genome.lib, character.only = TRUE)
  tbl <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, ...)
  ## check if the input file has the required information
  if(sum(!c("snp", "chr", "a1", "a2", "snpid") %in% names(tbl)) > 0) {
    stop("Error: 'filename' must be a table containing 'snp' and 'chr' columns.")
  }
  seqvec <- getSeq(Hsapiens,
                   as.character(tbl$chr),
                   start=tbl$snp - half.window.size,
                   end=tbl$snp + half.window.size,
                   as.character=TRUE)
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  sequences <- sapply(seqvec, function(x) codes[strsplit(x, "")[[1]]])
  colnames(sequences) <- tbl$snpid
  rownames(sequences) <- NULL
  a1 <- codes[tbl$a1]
  a2 <- codes[tbl$a2]
  names(a1) <- names(a2) <- NULL
  transition <- .Call("transition_matrix", sequences, package = "atSNP")
  prior <- apply(transition, 1, sum)
  prior <- prior / sum(prior)
  transition <- transition / apply(transition, 1, sum)
  names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  a1.ref.base.id <- which(a1 == sequences[half.window.size + 1, ])
  a2.ref.base.id <- which(a2 == sequences[half.window.size + 1, ])
  ## store SNPs that have the same base in SNP and REF alleles only once
  a2.ref.base.id <- a2.ref.base.id[!a2.ref.base.id %in% a1.ref.base.id]
  sequences <- sequences[, c(a1.ref.base.id, a2.ref.base.id)]
  ref.base <- c(a1[a1.ref.base.id], a2[a2.ref.base.id])
  snp.base <- c(a2[a1.ref.base.id], a1[a2.ref.base.id])
  return(list(
              sequence_matrix= sequences,
              ref_base = ref.base,
              snp_base = snp.base,
              transition = transition,
              prior = prior
              ))
}

#' @name ComputeMotifScore
#' @title Compute the scores for SNP effects on motifs.
#' @description Compute the log-likelihood scores for motifs.
#' @param motif.lib A list object with the output format of function 'LoadMotifLibrary'.
#' @param snp.info A list object with the output format of function 'LoadSNPData'.
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @details TODO.
#' @return A list of two data.table's. Field 'snp.tbl' contains:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleobase sequence.\cr
#' snp_seq \tab SNP allele nucleobase sequence.\cr
#' ref_seq_rev \tab Reference allele nucleobase sequence on the reverse strand.\cr
#' snp_seq_rev \tab SNP allele nucleobase sequence on the reverse strand.\cr}
#' Field 'motif.score' contains:
#' \tabular{cc}{
#' motif \tab Name of the motif.\cr
#' motif_len \tab Length of the motif.\cr
#' ref_start, ref_end, ref_strand \tab Location of the best matching subsequence on the reference allele.\cr
#' snp_start, snp_end, snp_strand \tab Location of the best matching subsequence on the SNP allele.\cr
#' log_lik_ref \tab Log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab Log-likelihood score for the SNP allele.\cr
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' log_enhance_odds \tab Difference in log-likelihood ratio between SNP allele and reference allele based on the best matching subsequence on the reference allele.\cr
#' log_reduce_odds \tab Difference in log-likelihood ratio between reference allele and SNP allele based on the best matching subsequence on the SNP allele.\cr
#' }
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' ComputeMotifScore(motif_library, snpInfo, ncores = 2)
#' @useDynLib atSNP
#' @import data.table doMC
#' @export
ComputeMotifScore <- function(motif.lib, snp.info, ncores = 1) {
  ## check arguments
  if(sum(!unlist(sapply(motif.lib, is.matrix))) > 0 | sum(unlist(sapply(motif.lib, ncol)) != 4) > 0) {
    stop("Error: 'motif.lib' must be a list of numeric matrices each with 4 columns.")
  }
  if(sum(!c("sequence_matrix", "snp_base", "ref_base") %in% names(snp.info)) > 0) {
    stop("Error: 'snp.info' must contain three components: 'ref_base', 'snp_base', 'sequence_matrix'.")
  }
  if(ncol(snp.info$sequence_matrix) != length(snp.info$ref_base) | length(snp.info$ref_base) != length(snp.info$snp_base)) {
    stop("Error: the number of columns of 'snp.info$sequence_matrix', the length of 'snp.info$ref_base' and the length of 'snp.info$snp_base' must be the same.")
  }
  if(sum(sort(unique(c(c(snp.info$sequence_matrix), snp.info$ref_base, snp.info$snp_base))) != seq(4)) > 0) {
    stop("Error: 'snp.info$sequence_matrix', 'snp.info$ref_base', 'snp.info$snp_base' can only contain entries in 1, 2, 3, 4.")
  }
  if(nrow(snp.info$sequence_matrix) / 2 == as.integer(nrow(snp.info$sequence_matrix) / 2)) {
    stop("Error: 'snp.info$sequence_matrix' must have an odd number of rows so that the central row refers to the SNP nucleotide.")
  }

  motifs <- names(motif.lib)
  snpids <- colnames(snp.info$sequence_matrix)
  nsnps <- ncol(snp.info$sequence_matrix)
  nmotifs <- length(motif.lib)
  
  motif_tbl <- data.table(motif = motifs,
                          motif_len = sapply(motif.lib, nrow))

  ncores <- min(c(ncores, length(snp.info$ref_base)))
  registerDoMC(ncores)

  motif_score_par <- foreach(i = seq(ncores)) %dopar% {
    k <- as.integer(length(snp.info$ref_base) / ncores)
    if(i < ncores) {
      ids <- seq(k) + k * (i - 1)
    } else {
      ids <- (k * (ncores - 1) + 1):length(snp.info$ref_base) 
    }
    this.snp.info <- list(sequence_matrix = snp.info$sequence_matrix[, ids],
                          ref_base = snp.info$ref_base[ids], snp_base = snp.info$snp_base[ids])
    motif.scores <- .Call("motif_score", motif.lib, this.snp.info, package = "atSNP")
    for(i in seq_along(motif.scores)) {
      rownames(motif.scores[[i]]) <- snpids[ids]
      colnames(motif.scores[[i]]) <- motifs
    }
    
    len_seq <- nrow(snp.info$sequence_matrix)
    strand_ref <- (motif.scores$match_ref_base > 0)
    ref_start <- motif.scores$match_ref_base
    ref_start[!strand_ref] <- len_seq + motif.scores$match_ref_base[!strand_ref] + 1
    
    strand_snp <- (motif.scores$match_snp_base > 0)
    snp_start <- motif.scores$match_snp_base
    snp_start[!strand_snp] <- len_seq + motif.scores$match_snp_base[!strand_snp] + 1
    
    motif_score_tbl <- data.table(snpid = rep(snpids[ids], nmotifs),
                                  motif = rep(motifs, each = length(ids)),
                                  log_lik_ref = c(motif.scores$log_lik_ref),
                                  log_lik_snp = c(motif.scores$log_lik_snp),
                                  log_lik_ratio = c(motif.scores$log_lik_ratio),
                                  log_enhance_odds = c(motif.scores$log_enhance_odds),
                                  log_reduce_odds = c(motif.scores$log_reduce_odds),
                                  ref_start = c(ref_start),
                                  snp_start = c(snp_start),
                                  ref_strand = c("-", "+")[1 + as.integer(strand_ref)],
                                  snp_strand = c("-", "+")[1 + as.integer(strand_snp)]
                                  )
    setkey(motif_score_tbl, motif)
    setkey(motif_tbl, motif)
    motif_score_tbl <- motif_tbl[motif_score_tbl]
    motif_score_tbl[ref_strand == "-", ref_start := ref_start - motif_len + 1]
    motif_score_tbl[, ref_end := ref_start + motif_len - 1]
    motif_score_tbl[snp_strand == "-", snp_start := snp_start - motif_len + 1]
    motif_score_tbl[, snp_end := snp_start + motif_len - 1]

    motif_score_tbl
    
  }

  motif.scores <- motif_score_par[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      motif.scores <- rbind(motif.scores, motif_score_par[[i]])
    }
  }

  ## sequences on the reference genome
  ref_seqs <- apply(snp.info$sequence_matrix, 2, function(x) paste(c("A", "C", "G", "T")[x], collapse = ""))
  ref_seqs_rev <- apply(snp.info$sequence_matrix, 2, function(x) paste(c("A", "C", "G", "T")[5 - rev(x)], collapse = ""))
  ## sequences on the snp allele
  id1 <- seq(as.integer(nrow(snp.info$sequence_matrix) / 2))
  id2 <- id1 + (nrow(snp.info$sequence_matrix) + 1) / 2
  snp_seqs <- apply(rbind(snp.info$sequence_matrix, snp.info$snp_base), 2,
                    function(x)
                    paste(c("A", "C", "G", "T")[x[c(id1, length(x), id2)]],
                          collapse = "")
                    )
  snp_seqs_rev <- apply(rbind(snp.info$sequence_matrix, snp.info$snp_base), 2,
                    function(x)
                    paste(c("A", "C", "G", "T")[5 - rev(x[c(id1, length(x), id2)])],
                          collapse = "")
                    )

  snp_tbl <- data.table(snpid = snpids,
                        ref_seq = ref_seqs,
                        snp_seq = snp_seqs,
                        ref_seq_rev = ref_seqs_rev,
                        snp_seq_rev = snp_seqs_rev)
  
  return(list(snp.tbl = snp_tbl, motif.scores = motif.scores))
}

#' @name MatchSubsequence
#' @title Compute the matching subsequence.
#' @description TODO.
#' @param snp.tbl A data.table with the following information:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleobase sequence.\cr
#' snp_seq \tab SNP allele nucleobase sequence.\cr
#' ref_seq_rev \tab Reference allele nucleobase sequence on the reverse strand.\cr
#' snp_seq_rev \tab SNP allele nucleobase sequence on the reverse strand.\cr}
#' @param motif.scores A data.table with the following information:
#' \tabular{cc}{
#' motif \tab Name of the motif.\cr
#' motif_len \tab Length of the motif.\cr
#' ref_start, ref_end, ref_strand \tab Location of the best matching subsequence on the reference allele.\cr
#' snp_start, snp_end, snp_strand \tab Location of the best matching subsequence on the SNP allele.\cr
#' log_lik_ref \tab Log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab Log-likelihood score for the SNP allele.\cr
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' log_enhance_odds \tab Difference in log-likelihood ratio between SNP allele and reference allele based on the best matching subsequence on the reference allele.\cr
#' log_reduce_odds \tab Difference in log-likelihood ratio between reference allele and SNP allele based on the best matching subsequence on the SNP allele.\cr
#' }
#' @param motif.lib A list of the position weight matrices for the motifs.
#' @param snpids A subset of snpids to compute the subsequences. Default: NULL, when all snps are computed.
#' @param motifs A subset of motifs to compute the subsequences. Default: NULL, when all motifs are computed.
#' @param ncores The number of cores used for parallel computing.
#' @return A data.table containing all columns in both 'snp.tbl' and 'motif.scores'. In addition, the following columns are added:
#' \tabular{ll}{
#' ref_match_seq \tab Best matching subsequence on the reference allele.\cr
#' snp_match_seq \tab Best matching subsequence on the SNP allele.\cr
#' ref_seq_snp_match \tab Subsequence on the reference allele corresponding to the best matching location on the SNP allele.\cr
#' snp_seq_ref_match \tab Subsequence on the SNP allele corresponding to the best matching location on the reference allele.\cr
#' }
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' MatchSubsequence(motif_scores$snp.tbl, motif_scores$motif.scores, motif_library)
#' @useDynLib atSNP
#' @import data.table doMC
#' @export
MatchSubsequence <- function(snp.tbl, motif.scores, motif.lib, snpids = NULL, motifs = NULL, ncores = 2) {
  if(is.null(snpids)) {
    snpids <- unique(snp.tbl$snpid)
  }
  if(is.null(motifs)) {
    motifs <- unique(motif.scores$motif)
  }
  if(sum(! motifs %in% names(motif.lib)) > 0) {
    stop("Error: some motifs are not included in 'motif.lib'.")
  }
  if(sum(! snpids %in% motif.scores$snpid) > 0) {
    stop("Error: some snpids are not included in 'motif.scores'.")
  }
  snpids <- unique(snpids)
  motifs <- unique(motifs)
  motif.scores <- motif.scores[snpid %in% snpids & motif %in% motifs, ]
  snp.tbl <- snp.tbl[snpid %in% snpids, ]
  snp.tbl[, len_seq := nchar(ref_seq)]

  ## get the IUPAC subsequence for the motifs
  motif.tbl <- data.table(
    motif = motifs,
    IUPAC = sapply(motif.lib[motifs],
      function(x) GetIUPACSequence(x, prob = 0.25))
  )
  setkey(motif.tbl, motif)

  ncores <- min(c(ncores, length(snpids)))
  registerDoMC(ncores)

  motif_score_par <- foreach(i = seq(ncores)) %dopar% {
    k <- as.integer(length(snpids) / ncores)
    if(i < ncores) {
      ids <- seq(k) + k * (i - 1)
    } else {
      ids <- (k * (ncores - 1) + 1):length(snpids) 
    }
    motif.scores <- motif.scores[snpid %in% snpids[ids], ]
    setkey(motif.scores, motif)
    motif.scores <- motif.tbl[motif.scores]
    setkey(motif.scores, snpid)
    setkey(snp.tbl, snpid)
    motif.scores <- snp.tbl[motif.scores]
    motif.scores[ref_strand == "+", ref_match_seq := substr(ref_seq, ref_start, ref_end)]
    motif.scores[ref_strand == "-", ref_match_seq := substr(ref_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
    motif.scores[snp_strand == "+", snp_match_seq := substr(snp_seq, snp_start, snp_end)]
    motif.scores[snp_strand == "-", snp_match_seq := substr(snp_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]
    motif.scores[ref_strand == "+", snp_seq_ref_match := substr(snp_seq, ref_start, ref_end)]
    motif.scores[ref_strand == "-", snp_seq_ref_match := substr(snp_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
    motif.scores[snp_strand == "+", ref_seq_snp_match := substr(ref_seq, snp_start, snp_end)]
    motif.scores[snp_strand == "-", ref_seq_snp_match := substr(ref_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]

    motif.scores[, list(snpid,
                        motif,
                        ref_seq,
                        snp_seq,
                        motif_len,
                        ref_start,
                        ref_end,
                        ref_strand,
                        snp_start,
                        snp_end,
                        snp_strand,
                        log_lik_ref,
                        log_lik_snp,
                        log_lik_ratio,
                        log_enhance_odds,
                        log_reduce_odds,
			IUPAC,
                        ref_match_seq,
                        snp_match_seq,
                        ref_seq_snp_match,
                        snp_seq_ref_match)]
  }

  motif_score_tbl <- motif_score_par[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      motif_score_tbl <- rbind(motif_score_tbl,
                               motif_score_par[[i]])
    }
  }
  
  return(motif_score_tbl)
}

#' @name ComputePValues
#' @title Compute p values.
#' @description TODO.
#' @param motif.lib A list object with the output format of function 'LoadMotifLibrary'.
#' @param snp.info A list object with the output format of function 'LoadSNPData'.
#' @param motif.scores A data.table object containing at least the following columns:
#' \tabular{ll}{
#' motif \tab The name of the motif.\cr
#' log_lik_ref \tab The log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab The log-likelihood score for the SNP allele.\cr
#' }
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @details TODO.
#' @return A data.table extending 'motif.scores' by the following additional columns:
#' \tabular{ll}{
#' pval_ref \tab P values for scores on the reference allele.\cr
#' pval_snp \tab P values for scores on the SNP allele.\cr
#' pval_diff \tab P values for the difference in scores between the reference and the SNP alleles.\cr
#' }
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' ComputePValues(motif_library, snpInfo, motif_scores$motif.scores, ncores = 4) 
#' @import doMC Rcpp data.table
#' @useDynLib atSNP
#' @export
ComputePValues <- function(motif.lib, snp.info, motif.scores, ncores = 1, getPlot = FALSE) {
  ncores <- min(c(ncores, length(motif.lib)))
    registerDoMC(ncores)
    results <- as.list(seq_along(motif.lib))
    nsets <- as.integer(length(motif.lib) / ncores)
    if(FALSE) {
        for(i in seq(nsets)) {
            message(i)
            if(i < nsets) {
                ids <- seq(ncores) + (i - 1) * ncores
            } else {
                ids <- ((nsets - 1) * ncores + 1) : length(motif.lib)
            }
        }
    }
    results <- foreach(motifid = seq_along(motif.lib)) %dopar% {
      rowids <- which(motif.scores$motif == names(motif.lib)[motifid])
      scores <- cbind(motif.scores$log_lik_ref[rowids],
                    motif.scores$log_lik_snp[rowids])
    pwm <- motif.lib[[motifid]]
    pwm[pwm < 1e-10] <- 1e-10
    wei.mat <- pwm
    for(i in seq(nrow(wei.mat))) {
      for(j in seq(ncol(wei.mat))) {
        wei.mat[i, j] <- exp(mean(log(pwm[i, j] / pwm[i, -j])))
      }
    }
    set.seed(motifid)
    
    p <- 5 / nrow(scores)
    
    m <- 20
    b <- (1 / p) ^ ( 1 / sum(seq(m)))
    allp <- rep(1, m + 1)
    step <- b
    for(k in rev(seq(m))) {
      allp[k] <- allp[k + 1] / step
      step <- step * b
    }
    allp <- allp[-(m + 1)]
    
    score.p <- quantile(c(scores), 1 - allp)
    
    pval_a <- pval_cond <- matrix(1, nrow = nrow(scores), ncol = 4)
    for(l in seq(length(allp) + 1)) {
      if(l == 1) {
        score.upp = max(scores) + 1
      } else if(l <= length(allp)) {
          score.upp = score.p[l - 1]
      } else {
          score.upp = quantile(c(scores), 0.2)
      }
      if(l >= length(allp)) {
        score.low = min(scores) - 1
      } else {
        score.low = score.p[l + 1]
      }
      compute.id <- which(scores < score.upp & scores >= score.low)
      if(length(compute.id) == 0) {
        next
      }
      if(l < length(allp) + 1) {
          theta <- .Call("test_find_theta", pwm, snp.info$prior, snp.info$transition, score.p[l], package = "atSNP")
      } else {
          theta <- 0
      }
      pval_a.new <- .Call("test_p_value", pwm, snp.info$prior, snp.info$transition, scores[compute.id], theta, package = "atSNP")
      pval_cond.new <- .structure(pval_a.new[, 4 + seq(4)])
      pval_a.new <- .structure(pval_a.new[, seq(4)])
      if(FALSE) {
        update.id <- which(scores < score.p[l] & scores >= score.low)
        pval_a[update.id] <- pval_a.new[compute.id %in% update.id, 1]
        if(l > 1) {
          update.id <- which(scores < score.upp & scores >= score.p[l])
          pval_a[update.id] <- (pval_a.new[compute.id %in% update.id, 1] + pval_a[update.id]) / 2
        }
      }
      if(TRUE) {
        update.id <- which(pval_a.new[, 2] < pval_a[, 3:4][compute.id])
        pval_a[compute.id[update.id]] <- pval_a.new[update.id, 1]
        pval_a[compute.id[update.id] + 2 * nrow(pval_a)] <- pval_a.new[update.id, 2]
        update.id <- which(pval_cond.new[, 2] < pval_cond[, 3:4][compute.id])
        pval_cond[compute.id[update.id]] <- pval_cond.new[update.id, 1]
        pval_cond[compute.id[update.id] + 2 * nrow(pval_cond)] <- pval_cond.new[update.id, 2]
      }
      if(FALSE) {
        if(length(update.id) > 0) {
          pval_a[update.id] <- pval_a.new[update.id]
          update.id <- update.id + nrow(pval_a) * 2
          pval_a[update.id] <- pval_a.new[update.id]
        }
      }
    }
    
    score_diff <- apply(scores, 1, function(x) abs(diff(x)))
    score.p <- round(quantile(score_diff, c((seq(8) + 1) / 10, 0.9 + seq(9) / 100)))
    score.p <- c(score.p,
                 seq(round(quantile(score_diff, 0.1) + 1),
                     round(quantile(score_diff, 0.9)),
                     by = 2)
                 )
    score.p <- rev(sort(unique(score.p)))
    
    pval_diff <- matrix(1, nrow = length(score_diff), ncol = 2)
    for(l in seq_along(score.p)) {
      if(l == 1) {
        score.upp <- max(score_diff) + 1
      } else {
        score.upp <- score.p[l - 1]
      }
      if(l == length(score.p)) {
        score.low <- min(score_diff) - 1
      } else {
        score.low <- score.p[l + 1]
      }
      compute.id <- which(score_diff < score.upp & score_diff >= score.low)
      if(length(compute.id) == 0) {
        next
      }
      pval_diff.new <- .Call("test_p_value_change", pwm,
                             wei.mat, pwm + 0.25, snp.info$prior,
                             snp.info$transition, score_diff[compute.id],
                             score.p[l], package = "atSNP")
      pval_diff.new <- .structure_diff(pval_diff.new)
      update.id <- which(pval_diff.new[, 2] < pval_diff[compute.id, 2])
      pval_diff[compute.id[update.id], ] <- pval_diff.new[update.id, ]
      ## print(summary(pval_diff.new[,1]))
    }
    
    ## Force the p-values to be increasing
    pval_a[, 1] <- sort(pval_a[, 1])[rank(-scores[,1])]
    pval_a[, 2] <- sort(pval_a[, 2])[rank(-scores[,2])]
    pval_cond[, 1] <- sort(pval_cond[, 1])[rank(-scores[,1])]
    pval_cond[, 2] <- sort(pval_cond[, 2])[rank(-scores[,2])]
    pval_diff[, 1] <- sort(pval_diff[,1])[rank(-score_diff)]
    pval_a[pval_a[, 1] > 1, 1] <- 1
    pval_a[pval_a[, 2] > 1, 2] <- 1
    pval_cond[pval_cond[, 1] > 1, 1] <- 1
    pval_cond[pval_cond[, 2] > 1, 2] <- 1
    pval_diff[pval_diff[, 1] > 1, 1] <- 1
    message("Finished testing the ", motifid, "th motif")
    ##    save(list = ls(), file = paste("/p/keles/ENCODE-CHARGE/volume2/SNP/test/motif", motifid, ".Rda", sep= ""))
    if(getPlot) {
      plotdat <- data.frame(
                            score = c(scores),
                            p.value = c(pval_a[, seq(2)]),
                            var = c(pval_a[, 3:4]),
                            Allele = rep(c("ref", "snp"), each = nrow(scores))
                            )
      plotdat.diff <- data.frame(
                                 score = score_diff,
                                 p.value = pval_diff[,1],
                                 var = pval_diff[,2]
                                 )
      plotdat <- unique(plotdat)
      plotdat.diff <- unique(plotdat.diff)
      localenv <- environment()
      options(warn = -1)
      pdf(paste("/p/keles/ENCODE-CHARGE/volume2/SNP/test/motif", motifid, ".pdf", sep = ""), width = 10, height = 10)
      id <- which(rank(plotdat$p.value[plotdat$Allele == "ref"]) <= 500)
      print(ggplot(aes(x = score, y = p.value), data = plotdat[plotdat$Allele == "ref", ], environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(motif.lib)[motifid], "ref"))) 
      id <- which(rank(plotdat$p.value[plotdat$Allele == "snp"]) <= 500)
      print(ggplot(aes(x = score, y = p.value), data = plotdat[plotdat$Allele == "snp", ], environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(motif.lib)[motifid], "SNP")))
      print(ggplot(aes(x = score, y = p.value), data = plotdat.diff, environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(motif.lib)[motifid], " Change")))
      dev.off()
    }
    
    list(rowids = rowids,
         pval_a = pval_a,
         pval_cond = pval_cond,
         pval_diff = pval_diff)
  }
  
  for(i in seq(length(results))) {
    motif.scores[results[[i]]$rowids, pval_ref := results[[i]]$pval_a[, 1]]
    motif.scores[results[[i]]$rowids, pval_snp := results[[i]]$pval_a[, 2]]
    motif.scores[results[[i]]$rowids, pval_cond_ref := results[[i]]$pval_cond[, 1]]
    motif.scores[results[[i]]$rowids, pval_cond_snp := results[[i]]$pval_cond[, 2]]
    motif.scores[results[[i]]$rowids, pval_diff := results[[i]]$pval_diff[, 1]]
  }
  return(motif.scores)
}

#' @name GetIUPACSequence
#' @title Get the IUPAC sequence of a motif.
#' @description Convert the posotion weight matrix of a motif to the IUPAC sequence.
#' @param pwm The position weight matrix, with the columns representing A, C, G, T.
#' @param prob The probability threshold. Default: 0.25.
#' @return A character string.
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' GetIUPACSequence(motif_library[[1]], prob = 0.2)
#' @export
GetIUPACSequence <- function(pwm, prob = 0.25) {
  iupac.table <-
    c(".", "A", "C", "M", "G", "R", "S", "V", "T", "W", "Y", "H", "K", "D", "B", "N")
  iupac.value <- (pwm >= prob) %*% c(1, 2, 4, 8) + 1
  return(paste(iupac.table[iupac.value], collapse = ""))
}
