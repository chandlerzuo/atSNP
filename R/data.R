#' @name motif_library
#' @title A sample motif library.
#' @description A list of the position weight matrices corresponding to motifs, loaded by 'data(example)'.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
#' @import testthat
NULL

#' @name snpInfo
#' @title A data set for SNP information.
#' @description This list object loaded by 'data(example)' contains three fields:
#' \tabular{ll}{
#' sequence_matrix \tab A sequence matrix, coded by 1-A, 2-C, 3-G, 4-T, with each column corresponding to a subsequence of 61 bp around one SNP.\cr
#' transition \tab The transition matrix used in Markov model.\cr
#' prior \tab The stationary distribution used in the Markov model.\cr
#' }
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name motif_scores
#' @title Scores for the sample snp data computed based on the motif data.
#' @description This list object loaded by 'data(example)' contains two fields:
#' \tabular{ll}{
#' snp.tbl \tab A data.table containing the sequence of nucleobases around each SNP.\cr
#' motif.scores \tab A data.table containing the likelihood scores computed for each SNP and each motif.\cr
#'}
#' @docType data
#' @format A data.table object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name motif_encode
#' @title A motif library provided by the ENCODE consortium.
#' @description This motif library is obtained from http://compbio.mit.edu/encode-motifs/motifs-toscan.txt.gz. Loaded by 'data(encode_motif)'.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name prior
#' @title Default stationary distribution for nucleotide sequences in the reference genome.
#' @description This parameter is fitted using 61bp windowns around the SNPs in the NHGRI catalog. Loaded by 'data(default_par)'.
#' @docType data
#' @format A numeric vector.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name transition
#' @title Default transition probability matrix for nucleotide sequences in the reference genome.
#' @description This parameter is fitted using 61bp windowns around the SNPs in the NHGRI catalog. Loaded by 'data(default_par)'.
#' @docType data
#' @format A 4 by 4 numeric matrix.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL
