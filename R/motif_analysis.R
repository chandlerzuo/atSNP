#' @name LoadMotifLibrary
#' @title Load position weight matrices.
#' @description Load the file for position weight matrices for motifs.
#' @param filename A file containing MEME format: http://meme.nbcr.net/meme/doc/meme-format.html.
#' @details This function reads the MEME formatted file containing motif information and convert them into a list of position weight matrices.
#' @return A list object of two components:
#' \tabular{ll}{
#' matrix \tab A list of position weight matrices.\cr
#' prior \tab A vector of the prior distribution for ACGT.\cr}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' pwms <- LoadMotifLibrary("/p/keles/ENCODE-CHARGE/volume1/ENCODE-Motifs/encode_motifs_for_fimo.txt")
#' }
#' @useDynLib MotifAnalysis
#' @export
LoadMotifLibrary <- function(filename) {
    lines <- readLines(filename)
    bkngLineNum <- grep("Background letter frequencies", lines) + 1
    priorProb <- as.numeric(strsplit(lines[bkngLineNum], " ")[[1]][seq(4) * 2])
    motifLineNums <- grep("MOTIF", lines)
    motifnames <-
        sapply(strsplit(lines[motifLineNums], " "), function(x) x[2])
    allmats <- as.list(seq_along(motifnames))

    for(matrixId in seq_along(motifLineNums)) {
        motifLineNum <- motifLineNums[matrixId]
        nrows <- as.integer(strsplit(lines[motifLineNum + 1], " ")[[1]][6])
        pwm <-
          t(matrix(as.numeric(unlist(strsplit(lines[seq(nrows) + motifLineNum + 1], " "))), nrow = 4))
        ## pwm <- read.table(filename, skip = motifLineNum + 1, nrows = nrows)
        allmats[[matrixId]] <- pwm
    }
    names(allmats) <- motifnames
    return(
        list(matrix = allmats,
             prior = priorProb))
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
#' @useDynLib MotifAnalysis
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
  transition <- .Call("transition_matrix", sequences, package = "MotifAnalysis")
  prior <- apply(transition, 1, sum)
  prior <- prior / sum(prior)
  transition <- transition / apply(transition, 1, sum)
  names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  return(list(
              sequence_matrix= sequences,
              a1 = a1,
              a2 = a2,
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
#' @return A list object with the following components:
#' \tabular{ll}{
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' on_odds \tab Enhanced log-odds for motif binding by SNP.\cr
#' off_odds \tab Reduced log-odds for motif binding by SNP.\cr
#' match_a1 \tab The starting positions of the best match for a1 sequences. Negative numbers correspond to the reverse strand.\cr
#' match_a2 \tab The starting positions of the best match for a2 sequences. Negative numbers correspond to the reverse strand.\cr
#' }
#' Each component is a matrix, with columns representing the motifs and the rows representing the SNPs.
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples \dontrun{}
#' @useDynLib MotifAnalysis
#' @export
ComputeMotifScore <- function(motif.lib, snp.info, ncores = 1) {
  ## check arguments
  if(sum(!c("prior", "matrix") %in% names(motif.lib)) > 0) {
    stop("Error: 'motif.lib' must contain both 'matrix' and'prior'.")
  }
  if(sum(!unlist(sapply(motif.lib$matrix, is.matrix))) > 0 | sum(unlist(sapply(motif.lib$matrix, ncol)) != 4) > 0) {
    stop("Error: 'motif.lib$matrix' must be a list of numeric matrices each with 4 columns.")
  }
  if(sum(!c("sequence_matrix", "a1", "a2") %in% names(snp.info)) > 0) {
    stop("Error: 'snp.info' must contain three components: 'a1', 'a2', 'sequence_matrix'.")
  }
  if(ncol(snp.info$sequence_matrix) != length(snp.info$a1) | length(snp.info$a1) != length(snp.info$a2)) {
    stop("Error: the number of columns of 'snp.info$sequence_matrix', the length of 'snp.info$a1' and the length of 'snp.info$a2' must be the same.")
  }
  if(sum(sort(unique(c(c(snp.info$sequence_matrix), snp.info$a1, snp.info$a2))) != seq(4)) > 0) {
    stop("Error: 'snp.info$sequence_matrix', 'snp.info$a1', 'snp.info$a2' can only contain entries in 1, 2, 3, 4.")
  }
  
  library(doMC)
  registerDoMC(ncores)

  motif_score_par <- foreach(i = seq(ncores)) %dopar% {
    k <- as.integer(length(snp.info$a1) / ncores)
    if(i < ncores) {
      ids <- seq(k) + k * (i - 1)
    } else {
      ids <- (k * (ncores - 1) + 1):length(snp.info$a1) 
    }
    this.snp.info <- list(sequence_matrix = snp.info$sequence_matrix[, ids],
                          a1 = snp.info$a1[ids], a2 = snp.info$a2[ids])
    .Call("motif_score", motif.lib, this.snp.info, package = "MotifAnalysis")
  }

  motif_scores <- motif_score_par[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      for(j in seq_along(motif_scores)) {
        motif_scores[[j]] <- rbind(motif_scores[[j]], motif_score_par[[i]][[j]])
      }
    }
  }
  for(i in seq_along(motif_scores)) {
    rownames(motif_scores[[i]]) <- colnames(snp.info$sequence_matrix)
    colnames(motif_scores[[i]]) <- names(motif.lib$matrix)
  }
  return(motif_scores)
  
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}
