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
#' motif.lib = motif_library)
#' @import data.table
#' @export
dtMotifMatch<-function(snp.tbl, motif.scores, snpids=NULL, motifs=NULL, 
                       motif.lib, ncores=2) {
  if(checkSNPids(snpids))
    {
      stop("snpids must be a vector of class character or NULL.")
    } else if (checkMotifs(motifs)) {
      stop("motifs must be a vector of class character or NULL.")
    }
  
  if(length(setdiff(snpids, motif.scores$snpid))!=0)
  {
    stop("snpids are not found in motif.scores.")
  } else if (length(setdiff(motifs, motif.scores$motif))!=0) {
    stop("motifs are not found in motif.scores.")
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
#' @param motif.match a single row ofdtMotifMatch output in data.frame format
#' @param motif.lib A list of position weight matrices
#' @param cex.main The size of the main title.
#' @param ... Other parameters passed to plotMotifLogo.
#' @return Sequence logo stacks: Reference subsequences, sequence logo of 
#' reference allele matching potision weight matrix, SNP subsequences, sequence 
#' logo of SNP allele matching potision weight matrix
#' @author Sunyoung Shin\email{sunyoung.shin@@utdallas.edu}
#' @examples
#' data(example)
#' plotMotifMatch(motif_match, motif.lib = motif_library)
#' @import grid
#' @importFrom motifStack plotMotifLogo pcm2pfm 
#' @importFrom grDevices graphics.off pdf
#' @importFrom stats quantile var
#' @importFrom utils data read.table write.table
#' @export
plotMotifMatch<-function(motif.match, motif.lib, cex.main = 2, ...) {
  if (!is(motif.match$snpid, "character") | length(motif.match$snpid)!=1) {
    stop("snpid must be a character")
  }
  if (!is(motif.match$motif, "character") | length(motif.match$motif)!=1) {
    stop("motif must be a character")
  }
  if(sum(! motif.match$motif %in% names(motif.lib)) > 0) {
    stop("The motif is not included in 'motif.lib'.")
  }

  if(nrow(motif.match)>1) {
      stop(paste("Pick a single row of dtMotifMatch output."))
  }
  
  motif.pwm <- t(get(motif.match$motif, motif.lib))
  ##Convert ACGT to 1234
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  ref_aug_match_seq_forward_code <- codes[strsplit(motif.match$ref_aug_match_seq_forward, "")[[1]]]
  ref_aug_match_seq_reverse_code <- codes[strsplit(motif.match$ref_aug_match_seq_reverse, "")[[1]]]
  snp_aug_match_seq_forward_code <- codes[strsplit(motif.match$snp_aug_match_seq_forward, "")[[1]]]
  snp_aug_match_seq_reverse_code <- codes[strsplit(motif.match$snp_aug_match_seq_reverse, "")[[1]]]
  
  ##Convert 1234 to (1000)(0100)(0010)(0001)
  codes.vec <- diag(4)
  rownames(codes.vec) <- c("A", "C", "G", "T")
  ref_aug_match_pwm_forward<- mapply(function(i) codes.vec[,i], as.list(ref_aug_match_seq_forward_code))
  ref_aug_match_pwm_reverse<- mapply(function(i) codes.vec[,i], as.list(ref_aug_match_seq_reverse_code))
  snp_aug_match_pwm_forward<- mapply(function(i) codes.vec[,i], as.list(snp_aug_match_seq_forward_code))
  snp_aug_match_pwm_reverse<- mapply(function(i) codes.vec[,i], as.list(snp_aug_match_seq_reverse_code))
  
  ##(3,2) to Augmented PWM: ___PWM__
  ref_aug_pwm <- cbind(matrix(0, 4, motif.match$ref_extra_pwm_left), motif.pwm, matrix(0, 4, motif.match$ref_extra_pwm_right))
  rownames(ref_aug_pwm) <- c("A", "C", "G", "T")
  snp_aug_pwm <- cbind(matrix(0, 4, motif.match$snp_extra_pwm_left), motif.pwm, matrix(0, 4, motif.match$snp_extra_pwm_right))
  rownames(snp_aug_pwm) <- c("A", "C", "G", "T")

  snp_loc <- motif.match$ref_location
  revert.columns <- function(mat) {
    mat[, rev(seq(ncol(mat)))]
  }
  
  ref_aug_match_pwm<-ref_aug_match_pwm_forward
  snp_aug_match_pwm<-snp_aug_match_pwm_forward

  if(motif.match$ref_strand == "-") {
    ref_aug_pwm <- revert.columns(ref_aug_pwm)
    snp_loc <- ncol(ref_aug_match_pwm_forward) - 1 - snp_loc
    ref_aug_match_pwm<-ref_aug_match_pwm_reverse
  }
  if(motif.match$snp_strand == "-") {
    snp_aug_pwm <- revert.columns(snp_aug_pwm)
    snp_aug_match_pwm<-snp_aug_match_pwm_reverse
  }
  
  pushViewport(viewport(y =unit(.5, "npc") - unit(2, "lines"), height = unit(1, "npc") - unit(3, "lines")))
  pushViewport(viewport(y=.875, height=.25))
  plotMotifLogo(pcm2pfm(ref_aug_pwm), "Best match to the reference genome", yaxis=FALSE,
                xaxis=FALSE, xlab="", ylab="PWM", newpage=FALSE, margins = c(1.5, 3, 2, 2))
  if(motif.match$ref_strand=='+') {
    grid.lines(x=c(convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE), 1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)),
               y=unit(1, "lines"),
               gp=gpar(col = "blue", lwd = 1.5, xpd=NA), arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "last"))
    grid.text("3'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
    grid.text("5'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
  } else {
    grid.lines(x=c(convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE), 1 - convertUnit(unit(2, "lines"), "npc", valueOnly = TRUE)),
               y=unit(1, "lines"), gp=gpar(col = "blue", lwd = 1.5, xpd=NA),
               arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "first"))
    grid.text("5'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
    grid.text("3'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
  }
  popViewport()
  pushViewport(viewport(y=.625, height=.25))
  #par(mar = c(4, 3, 1.5, 2))
  plotMotifLogo(pcm2pfm(ref_aug_match_pwm), font="mono,Courier", yaxis=FALSE, xlab="",
                ylab=paste("(", motif.match$ref_strand, ")", sep=""),  newpage=FALSE, margins = c(2, 3, 1.5, 2))
  pushViewport(plotViewport(margins = c(2, 3, 1.5, 2)))
  grid.rect(x=(snp_loc+.5)/motif.match$snp_ref_length, width = 1/motif.match$snp_ref_length, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
  popViewport()
  if(motif.match$ref_strand=="+")   {
    grid.text("3'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(2.5, "lines"))
    grid.text("5'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(2.5, "lines"))
  } else {
    grid.text("5'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(2.5, "lines"))
    grid.text("3'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(2.5, "lines"))
  }
  popViewport()
  pushViewport(viewport(y=.375, height=.25))
  #par(mar=c(1.5, 3, 4, 2))     
  plotMotifLogo(pcm2pfm(snp_aug_match_pwm), "Best match to the SNP genome", font="mono,Courier",
                yaxis=FALSE, xlab="", ylab=paste("(", motif.match$snp_strand, ")", sep=""), newpage=FALSE, margins = c(1.5, 3, 2, 2))
  pushViewport(plotViewport(margins = c(1.5, 3, 2, 2)))
  grid.rect(x=(snp_loc+.5)/motif.match$snp_ref_length, width = 1/motif.match$snp_ref_length, gp=gpar(col="blue", lty=3, lwd=2, fill=NA))
  popViewport()
  if(motif.match$snp_strand=="+")   {
    grid.text("3'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
    grid.text("5'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
  } else {
    grid.text("5'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
    grid.text("3'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(.5, "lines"))
  }
  popViewport()
  pushViewport(viewport(y=.125, height=.25))
  #par(mar=c(4, 3, 1.5, 2))
  plotMotifLogo(pcm2pfm(snp_aug_pwm), yaxis=FALSE, xaxis=FALSE, xlab="", ylab="PWM", newpage=FALSE, margins = c(2, 3, 1.5, 2))
  if(motif.match$snp_strand=='+') {
    grid.lines(x=c(convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE), 1 - convertUnit(unit(1, "lines"), "npc", valueOnly = TRUE)),
               y=unit(1.5, "lines"), gp=gpar(col = "blue", lwd = 1.5, xpd=NA),
               arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "last"))
    grid.text("3'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(1, "lines"))
    grid.text("5'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(1, "lines"))
  } else {
    grid.lines(x=c(convertUnit(unit(3, "lines"), "npc", valueOnly = TRUE), 1 - convertUnit(unit(1, "lines"), "npc", valueOnly = TRUE)),
               y=unit(1.5, "lines"), gp=gpar(col = "blue", lwd = 1.5, xpd=NA),
               arrow=arrow(length = unit(0.1, "inches"), angle = 15, ends = "first"))
    grid.text("5'", x=unit(1, "npc") - unit(1, "lines"), gp=gpar(col="blue", cex=1), y=unit(1, "lines"))
    grid.text("3'", x=unit(2, "lines"), gp=gpar(col="blue", cex=1), y=unit(1, "lines"))
  }
  popViewport()
  popViewport()
  grid.text(label = paste(motif.match$motif, " Motif Scan for ", motif.match$snpid, sep=""),
            y=unit(1, "npc") - unit(1.5, "lines"),
            gp=gpar(cex.main=cex.main, fontface="bold"))
  graphics.off()
  
}

.find_reverse <- function(sequence) {
  if(length(sequence) > 0) {
    codes <- seq(4)
    names(codes) <- c("A", "C", "G", "T")
    return(paste(names(codes)[5 - codes[strsplit(sequence, split = "")[[1]]]], collapse = ""))
  }
}
