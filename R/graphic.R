#' @name dtMotifMatch
#' @title Compute the augmented matching subsequence on SNP and reference allele
#' s.
#' @description Calculate the best matching augmented subsequences on both SNP 
#' and reference alleles for motifs. Obtain extra unmatching position on the 
#' best matching augmented subsequence of the reference and SNP alleles.
#' @param motif.lib A list of named position weight matrices.
#' @param snp.tbl A data.frame with the following information:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleobase sequence.\cr
#' snp_seq \tab SNP allele nucleobase sequence.\cr
#' ref_seq_rev \tab Reference allele nucleobase sequence on the reverse 
#' strand.\cr
#' snp_seq_rev \tab SNP allele nucleobase sequence on the reverse strand.\cr}
#' @param motif.scores A data.frame with the following information:
#' \tabular{cc}{
#' motif \tab Name of the motif.\cr
#' motif_len \tab Length of the motif.\cr
#' ref_start, ref_end, ref_strand \tab Location of the best matching subsequence
#'  on the reference allele.\cr
#' snp_start, snp_end, snp_strand \tab Location of the best matching subsequence
#'  on the SNP allele.\cr
#' log_lik_ref \tab Log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab Log-likelihood score for the SNP allele.\cr
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' log_enhance_odds \tab Difference in log-likelihood ratio between SNP allele 
#' and reference allele based on the best matching subsequence on the reference 
#' allele.\cr
#' log_reduce_odds \tab Difference in log-likelihood ratio between reference 
#' allele and SNP allele based on the best matching subsequence on the SNP 
#' allele.\cr
#' }
#' @param snpids A subset of snpids to compute the subsequences. Default: NULL, 
#' when all snps are computed.
#' @param motifs A subset of motifs to compute the subsequences. Default: NULL, 
#' when all motifs are computed.
#' @param ncores The number of cores used for parallel computing. Default: 10
#' @return A data.frame containing all columns from the function, 
#' \code{\link{MatchSubsequence}}. In addition, the following columns are added:
#' \tabular{ll}{
#' snp_ref_start, snp_ref_end, snp_ref_length \tab Location and Length of the 
#' best matching augmented subsequence on both the reference and SNP allele.\cr
#' ref_aug_match_seq_forward \tab Best matching augmented subsequence or its 
#' corresponding sequence to the forward strand on the reference allele.\cr 
#' snp_aug_match_seq_forward \tab Best matching augmented subsequence or its 
#' corresponding sequence to the forward strand on the SNP allele.\cr 
#' ref_aug_match_seq_reverse \tab Best matching augmented subsequence or its 
#' corresponding sequence to the reverse strand on the reference allele.\cr 
#' snp_aug_match_seq_reverse \tab Best matching augmented subsequence or its 
#' corresponding sequence to the reverse strand on the SNP allele.\cr 
#' ref_location \tab SNP location of the best matching augmented subsequence on 
#' the reference allele. Starting from zero. \cr
#' snp_location \tab SNP location of the best matching augmented subsequence on 
#' the SNP allele. Starting from zero. \cr
#' ref_extra_pwm_left \tab Left extra unmatching position on the best matching 
#' augmented subsequence of the reference allele. \cr
#' ref_extra_pwm_right \tab Right extra unmatching position on the best matching
#'  augmented subsequence of the reference allele. \cr
#' snp_extra_pwm_left \tab Left extra unmatching position on the best matching 
#' augmented subsequence of the SNP allele. \cr
#' snp_extra_pwm_right \tab Right extra unmatching position on the best matching
#'  augmented subsequence of the SNP allele. \cr
#' }
#' @author Sunyoung Shin\email{sunyoung.shin@@utdallas.edu}
#' @examples
#' data(example)
#' dtMotifMatch(motif_scores$snp.tbl, motif_scores$motif.scores, 
#' snpids=motif_scores$snp.tbl$snpid, motifs=motif_scores$motif.scores$motif[1],
#'  motif.lib = motif_library)
#' @import data.table
#' @export
dtMotifMatch<-function(snp.tbl, motif.scores, snpids=NULL, motifs=NULL, 
                       motif.lib, ncores=2)
  {
  if(checkSNPids(snpids))
    {
      stop("snpids must be a vector of class character or NULL.")
    } else if (checkMotifs(motifs)) {
      stop("motifs must be a vector of class character or NULL.")
    }
  #warning for ncores, motif.lib etc.
  snp.tbl<-as.data.table(snp.tbl)
  ncores.v1 <- min(ncores, length(snpids) * length(motifs))
  ncores.v2<-ifelse(ncores.v1==0, ncores, ncores.v1)
  sequence.half.window.size <- (nchar(snp.tbl[1, ref_seq]) - 1) / 2
  motif.match <- MatchSubsequence(snp.tbl = snp.tbl, motif.scores = motif.scores, snpids = snpids, motifs = motifs, ncores = ncores.v2, motif.lib = motif.lib)
  motif.match.dt<-as.data.table(motif.match)
  ##Augmentation of SNP and reference sequences###
  motif.match.dt[, len_seq := nchar(ref_seq)]
  motif.match.dt[,snp_ref_start := apply(cbind(ref_start, snp_start), 1, min)]
  motif.match.dt[,snp_ref_end := apply(cbind(ref_end, snp_end), 1, max)]
  motif.match.dt[,snp_ref_length := snp_ref_end - snp_ref_start + 1]
  
  motif.match.dt[, ref_aug_match_seq_forward := substr(ref_seq, snp_ref_start, snp_ref_end)]
  motif.match.dt[, ref_aug_match_seq_reverse:= apply(as.matrix(ref_aug_match_seq_forward), 1, .find_reverse)]
  motif.match.dt[, snp_aug_match_seq_forward := substr(snp_seq, snp_ref_start, snp_ref_end)]
  motif.match.dt[, snp_aug_match_seq_reverse:= apply(as.matrix(snp_aug_match_seq_forward), 1, .find_reverse)]
  
  ##The starting position of the motif in the augmented sequences
  motif.match.dt[ref_strand == "+", ref_location := (len_seq-1)/2 + 1 - snp_ref_start]
  motif.match.dt[ref_strand == "-", ref_location := snp_ref_end - (len_seq - 1) / 2 - 1]
  motif.match.dt[snp_strand=="+", snp_location:=(len_seq-1)/2+1-snp_ref_start]
  motif.match.dt[snp_strand=="-", snp_location:=snp_ref_end-(len_seq-1)/2-1]
  motif.match.dt[, len_seq := NULL]
  
  ##PWM Location Adjustment Value for reference and SNP
  motif.match.dt[ref_strand == "+", ref_extra_pwm_left := ref_start-snp_ref_start]
  motif.match.dt[ref_strand == "-", ref_extra_pwm_left := snp_ref_end-ref_end]
  motif.match.dt[ref_strand == "+", ref_extra_pwm_right := snp_ref_end-ref_end]
  motif.match.dt[ref_strand == "-", ref_extra_pwm_right := ref_start-snp_ref_start]
  motif.match.dt[snp_strand == "+", snp_extra_pwm_left := snp_start-snp_ref_start]
  motif.match.dt[snp_strand == "-", snp_extra_pwm_left := snp_ref_end-snp_end]
  motif.match.dt[snp_strand == "+", snp_extra_pwm_right := snp_ref_end-snp_end]
  motif.match.dt[snp_strand == "-", snp_extra_pwm_right := snp_start-snp_ref_start]
  setkey(motif.match.dt, snpid)
  return(motif.match.dt)
}

#' @name plotMotifMatch
#' @title Plot sequence logos of the position weight matrix of the motif and 
#' sequences of its corresponding best matching augmented subsequence on the 
#' reference and SNP allele.
#' @description Plot the best matching augmented subsequences on the reference 
#' and SNP alleles. Plot sequence logos of the position weight matrix of the 
#' motif to the corresponding positions of the best matching subsequences on the
#'  references and SNP alleles.
#' @param snp.tbl A data.frame with the following information:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleobase sequence.\cr
#' snp_seq \tab SNP allele nucleobase sequence.\cr
#' ref_seq_rev \tab Reference allele nucleobase sequence on the reverse 
#' strand.\cr
#' snp_seq_rev \tab SNP allele nucleobase sequence on the reverse strand.\cr}
#' @param motif.scores A data.frame with the following information:
#' \tabular{cc}{
#' motif \tab Name of the motif.\cr
#' motif_len \tab Length of the motif.\cr
#' ref_start, ref_end, ref_strand \tab Location of the best matching subsequence
#'  on the reference allele.\cr
#' snp_start, snp_end, snp_strand \tab Location of the best matching subsequence
#'  on the SNP allele.\cr
#' log_lik_ref \tab Log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab Log-likelihood score for the SNP allele.\cr
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' log_enhance_odds \tab Difference in log-likelihood ratio between SNP allele 
#' and reference allele based on the best matching subsequence on the reference 
#' allele.\cr
#' log_reduce_odds \tab Difference in log-likelihood ratio between reference 
#' allele and SNP allele based on the best matching subsequence on the SNP 
#' allele.\cr
#' }
#' @param snpid A snpid to plot the sequences on the reference and SNP alleles
#' @param snp snp nucleotide of the snpid
#' @param motif A motif to match the sequences with its position weight matrix
#' @param motif.lib A list of position weight matrices
#' @param cex.main The size of the main title.
#' @param ... Other parameters passed to plotMotifLogo.
#' @return Sequence logo stacks: Reference subsequences, sequence logo of 
#' reference allele matching potision weight matrix, SNP subsequences, sequence 
#' logo of SNP allele matching potision weight matrix
#' @author Sunyoung Shin\email{sunyoung.shin@@utdallas.edu}
#' @examples
#' data(example)
#' plotMotifMatch(motif_scores$snp.tbl, motif_scores$motif.scores, 
#' snpid=motif_scores$snp.tbl$snpid[1], 
#' motif=motif_scores$motif.scores$motif[1], motif.lib = motif_library)
#' @import data.table 
#' @importFrom motifStack plotMotifLogo pcm2pfm 
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics arrows mtext par segments title
#' @importFrom stats quantile var
#' @importFrom utils data read.table write.table
#' @export
plotMotifMatch<-function(snp.tbl, motif.scores, snpid, snp=NULL, motif, motif.lib, cex.main = 2, ...) {
  if (!is(snpid, "character") | length(snpid)!=1) {
    stop("snpid must be a character")
  }
  if (!is(motif, "character") | length(motif)!=1) {
    stop("motif must be a character")
  }
  if(sum(! motif %in% names(motif.lib)) > 0) {
    stop("The motif is not included in 'motif.lib'.")
  }

  motif.match.dt <- dtMotifMatch(snp.tbl, motif.scores, snpids=snpid, motifs=motif, ncores = 1, motif.lib = motif.lib)  
if(nrow(motif.match.dt)>1) {
      if(is.null(snp)) {
      stop(paste("Multiple nucleobases on the SNP genome. Add snp= ", paste(motif.match.dt$snpbase, collapse=" or "), sep=""))
  } else {
  motif.match.dt<-motif.match.dt[snpbase==snp]
  }
  }
  
  ##Convert ACGT to 1234
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  ref_aug_match_seq_forward_code <- codes[strsplit(motif.match.dt[,ref_aug_match_seq_forward], "")[[1]]]
  ref_aug_match_seq_reverse_code <- codes[strsplit(motif.match.dt[,ref_aug_match_seq_reverse], "")[[1]]]
  snp_aug_match_seq_forward_code <- codes[strsplit(motif.match.dt[,snp_aug_match_seq_forward], "")[[1]]]
  snp_aug_match_seq_reverse_code <- codes[strsplit(motif.match.dt[,snp_aug_match_seq_reverse], "")[[1]]]
  
  ##Convert 1234 to (1000)(0100)(0010)(0001)
  codes.vec <- diag(4)
  rownames(codes.vec) <- c("A", "C", "G", "T")
  ref_aug_match_pwm_forward<- mapply(function(i) codes.vec[,i], as.list(ref_aug_match_seq_forward_code))
  ref_aug_match_pwm_reverse<- mapply(function(i) codes.vec[,i], as.list(ref_aug_match_seq_reverse_code))
  snp_aug_match_pwm_forward<- mapply(function(i) codes.vec[,i], as.list(snp_aug_match_seq_forward_code))
  snp_aug_match_pwm_reverse<- mapply(function(i) codes.vec[,i], as.list(snp_aug_match_seq_reverse_code))
  
  ##(3,2) to Augmented PWM: ___PWM__
  ref_aug_pwm <- cbind(matrix(0, 4, motif.match.dt[, ref_extra_pwm_left]), t(get(motif.match.dt[, motif], motif.lib)), matrix(0, 4, motif.match.dt[, ref_extra_pwm_right]))
  rownames(ref_aug_pwm) <- c("A", "C", "G", "T")
  snp_aug_pwm <- cbind(matrix(0, 4, motif.match.dt[, snp_extra_pwm_left]), t(get(motif.match.dt[, motif], motif.lib)), matrix(0, 4, motif.match.dt[, snp_extra_pwm_right]))
  rownames(snp_aug_pwm) <- c("A", "C", "G", "T")

  snp_loc <- motif.match.dt$ref_location
  revert.columns <- function(mat) {
    mat[, rev(seq(ncol(mat)))]
  }

  ref_aug_match_pwm<-ref_aug_match_pwm_forward
  snp_aug_match_pwm<-snp_aug_match_pwm_forward

  if(motif.match.dt$ref_strand == "-") {
    ref_aug_pwm <- revert.columns(ref_aug_pwm)
    snp_loc <- ncol(ref_aug_match_pwm_forward) - 1 - snp_loc
    ref_aug_match_pwm<-ref_aug_match_pwm_reverse
  }
  if(motif.match.dt$snp_strand == "-") {
    snp_aug_pwm <- revert.columns(snp_aug_pwm)
    snp_aug_match_pwm<-snp_aug_match_pwm_reverse
  }
  
  par(mfrow=c(4,1), oma=c(1,1,4,1))
  par(mar=c(1.5, 3, 4, 2))
  plotMotifLogo(pcm2pfm(ref_aug_pwm), "Best match to the reference genome", yaxis=FALSE, xaxis=FALSE, xlab="", ylab="PWM", ...)
if(motif.match.dt$ref_strand=='+') {
arrows((min(which(colSums(ref_aug_pwm)!=0))-1)/ncol(ref_aug_pwm), -0.17, max(which(colSums(ref_aug_pwm)!=0))/ncol(ref_aug_pwm), -0.17, length = 0.1, angle = 15, code = 2, col = "blue", lwd = 1.5, xpd=NA)
  mtext("5'", 1, adj=(min(which(colSums(ref_aug_pwm)!=0))-1)/ncol(ref_aug_pwm), padj=1, col="blue", cex=1) 
  mtext("3'", 1, adj=max(which(colSums(ref_aug_pwm)!=0))/ncol(ref_aug_pwm), padj=1, col="blue", cex=1)
} else {
arrows(max(which(colSums(ref_aug_pwm)!=0))/ncol(ref_aug_pwm), -0.17, (min(which(colSums(ref_aug_pwm)!=0))-1)/ncol(ref_aug_pwm), -0.17, length = 0.1, angle = 15, code = 2, col = "blue", lwd = 1.5, xpd=NA)
  mtext("5'", 1, adj=max(which(colSums(ref_aug_pwm)!=0))/ncol(ref_aug_pwm), padj=1, col="blue", cex=1) 
  mtext("3'", 1, adj=(min(which(colSums(ref_aug_pwm)!=0))-1)/ncol(ref_aug_pwm), padj=1, col="blue", cex=1) 
}
  par(mar = c(4, 3, 1.5, 2))
plotMotifLogo(pcm2pfm(ref_aug_match_pwm), font="mono,Courier", yaxis=FALSE, xlab="", ylab=paste("(", motif.match.dt$ref_strand, ")", sep=""), ...)
segments(snp_loc/motif.match.dt[,snp_ref_length], 0, snp_loc/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments(snp_loc/motif.match.dt[,snp_ref_length], 1, (snp_loc+1)/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments((snp_loc+1)/motif.match.dt[,snp_ref_length], 0, (snp_loc+1)/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments(snp_loc/motif.match.dt[,snp_ref_length], 0, (snp_loc+1)/motif.match.dt[,snp_ref_length], 0, col="blue", lty=3, lwd=2)
  if(motif.match.dt$ref_strand=="+")   {
  mtext("5'", 1,  adj=0, padj=1, col="blue", cex=1) 
  mtext("3'", 1,  adj=1, padj=1, col="blue", cex=1)
} else {
  mtext("3'", 1, adj=0, padj=1, col="blue", cex=1) 
  mtext("5'", 1, adj=1, padj=1, col="blue", cex=1) 
  }
par(mar=c(1.5, 3, 4, 2))      
plotMotifLogo(pcm2pfm(snp_aug_match_pwm), "Best match to the SNP genome", font="mono,Courier", yaxis=FALSE, xlab="", ylab=paste("(", motif.match.dt$snp_strand, ")", sep=""), ...)
segments(snp_loc/motif.match.dt[,snp_ref_length], 0, snp_loc/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments(snp_loc/motif.match.dt[,snp_ref_length], 1, (snp_loc+1)/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments((snp_loc+1)/motif.match.dt[,snp_ref_length], 0, (snp_loc+1)/motif.match.dt[,snp_ref_length], 1, col="blue", lty=3, lwd=2)
segments(snp_loc/motif.match.dt[,snp_ref_length], 0, (snp_loc+1)/motif.match.dt[,snp_ref_length], 0, col="blue", lty=3, lwd=2)
  if(motif.match.dt$snp_strand=="+")   {
  mtext("5'", 1,  adj=0, padj=1, col="blue", cex=1) 
  mtext("3'", 1,  adj=1, padj=1, col="blue", cex=1)
} else {
  mtext("3'", 1, adj=0, padj=1, col="blue", cex=1) 
  mtext("5'", 1, adj=1, padj=1, col="blue", cex=1) 
  }
par(mar=c(4, 3, 1.5, 2))
plotMotifLogo(pcm2pfm(snp_aug_pwm), yaxis=FALSE, xaxis=FALSE, xlab="", ylab="PWM", ...)
if(motif.match.dt$snp_strand=='+') {
arrows((min(which(colSums(snp_aug_pwm)!=0))-1)/ncol(snp_aug_pwm), -0.17, max(which(colSums(snp_aug_pwm)!=0))/ncol(snp_aug_pwm), -0.17, length = 0.1, angle = 15, code = 2, col = "blue", lwd = 1.5, xpd=NA)
  mtext("5'", 1, adj=(min(which(colSums(snp_aug_pwm)!=0))-1)/ncol(snp_aug_pwm), padj=1, col="blue", cex=1) 
  mtext("3'", 1, adj=max(which(colSums(snp_aug_pwm)!=0))/ncol(snp_aug_pwm), padj=1, col="blue", cex=1)
} else {
arrows(max(which(colSums(snp_aug_pwm)!=0))/ncol(snp_aug_pwm), -0.17, (min(which(colSums(snp_aug_pwm)!=0))-1)/ncol(snp_aug_pwm), -0.17, length = 0.1, angle = 15, code = 2, col = "blue", lwd = 1.5, xpd=NA)
  mtext("5'", 1, adj=max(which(colSums(snp_aug_pwm)!=0))/ncol(snp_aug_pwm), padj=1, col="blue", cex=1) 
  mtext("3'", 1, adj=(min(which(colSums(snp_aug_pwm)!=0))-1)/ncol(snp_aug_pwm), padj=1, col="blue", cex=1)
}
title(main=paste(motif.match.dt[,motif], " Motif Scan for ", motif.match.dt[,snpid], sep=""), outer=TRUE, cex.main=cex.main)
}

.find_reverse <- function(sequence) {
  if(length(sequence) > 0) {
    codes <- seq(4)
    names(codes) <- c("A", "C", "G", "T")
    return(paste(names(codes)[5 - codes[strsplit(sequence, split = "")[[1]]]], collapse = ""))
  }
}
