#' @name motif_library
#' @title A sample motif library.
#' @description A list of the position weight matrices corresponding to motifs.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
#' @import testthat
NULL

#' @name snpInfo
#' @title A data set for SNP information.
#' @description This data contains three fields:
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
#' @description This list object contains two fields:
#' \tabular{ll}{
#' snp.tbl \tab A data.table containing the sequence of nucleobases around each SNP.\cr
#' motif.scores \tab A data.table containing the likelihood scores computed for each SNP and each motif.\cr
#'}
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL
