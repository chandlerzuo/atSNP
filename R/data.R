#' @name motif_library
#' @title A sample motif library.
#' @description A list of the position weight matrices corresponding to motifs, loaded by 'data(example)'.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
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

#' @name snp_tbl
#' @title A data frame for SNP information.
#' @description This data frame is loaded by 'data(example)'. It is a table including the following columns:
#' \tabular{ll}{
#' chr \tab The chromosome.\cr
#' snp \tab The SNP location coordinate.\cr
#' snpid \tab The SNP label.\cr
#' a1,a2 \tab The nucleotide on the reference and SNP allele.\cr
#' }
#' @docType data
#' @format A data.frame object.
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

#' @name encode_motif
#' @title A motif library containing 2065 motifs downloaded from http://compbio.mit.edu/encode-motifs/motifs.txt.
#' @description This motif library can be loaded by 'data(encode_library)'.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name encode_motifinfo
#' @title The information for the motif library downloaded from http://compbio.mit.edu/encode-motifs/motifs.txt.
#' @description This is a character vector that be loaded by 'data(encode_library)'. The names of this vector are the same as the names for 'encode_motif'. The entries of this vector are the corresponding motif information parsed from the raw file.
#' @docType data
#' @format A character vector.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name jaspar_motif
#' @title A motif library containing 593 motifs downloaded from http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt.
#' @description This motif library can be loaded by 'data(jaspar_library)'.
#' @docType data
#' @format A list object.
#' @author Chandler Zuo zuo@@stat.wisc.edu
NULL

#' @name jaspar_motifinfo
#' @title The information for the motif library downloaded from http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt.
#' @description This is a character vector that be loaded by 'data(jaspar_library)'. The names of this vector are the same as the names for 'jaspar_motif'. The entries of this vector are the corresponding motif information parsed from the raw file.
#' @docType data
#' @format A character vector.
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
