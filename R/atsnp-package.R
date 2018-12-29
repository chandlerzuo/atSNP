#' atSNP: affinity tests for regulatory SNP detection
#'
#' @description atSNP implements  the affinity test for large sets of SNP-motif interactions using the importance sampling algorithm.
#' Users may identify SNPs that potentially may affect binding affinity of transcription factors. 
#' Given a set of SNPs and a library of motif position weight matrices (PWMs), 
#' atSNP provides two main functions for analyzing SNP effects:
#' (i) the binding affinity score for each allele and each PWM and 
#' the p-values for allele-specific binding affinity scores 
#' (ii) the p-values for affinity score changes between the two alleles for each SNP.
#' Compared to other bioinformatics tools that provide similar functionalities, 
#' atSNP is highly scalable.
#' 
#' The atSNP main functions are:
#' \enumerate{
#' \item \code{\link{LoadMotifLibrary}} - Load position weight matrices
#' \item \code{\link{LoadSNPData}} - Load the SNP information and code the genome sequences around the SNP locations
#' \item \code{\link{LoadFastaData}} - Load the SNP data from fasta files
#' \item \code{\link{ComputeMotifScore}} - Compute the scores for SNP effects on motifs
#' \item \code{\link{ComputePValues}} - Compute p-values for affinity scores
#' }
#'
#' Some helper functions are:
#' \enumerate{
#' \item \code{\link{MatchSubsequence}} - Compute the matching subsequence
#' \item \code{\link{GetIUPACSequence}} - Get the IUPAC sequence of a motif
#' \item \code{\link{dtMotifMatch}} - Compute the augmented matching subsequence on SNP and reference alleles
#' }
#'
#' The composite logo plotting function is:
#' \enumerate{
#' \item \code{\link{plotMotifMatch}} - Plot sequence logos of the position weight matrix of the motif and sequences of its corresponding best matching augmented subsequence on the reference and SNP allele
#' }
#'
#' @references
#' Zuo, Chandler, Shin, Sunyoung, and Keles, Sunduz. (2015). atSNP: Transcription factor binding affinity testing for regulatory SNP detection. Bioinformatics 31 (20): 3353-5.
#' 
#' @name atSNP-package
#' @aliases atSNP-package
#' @docType package
#' @author Chandler Zuo Sunyoung Shin \email{sunyoung.shin@@utdallas.edu} 
#' @keywords GenomeAnnotation MotifAnnotation LogoPlot
#' @importFrom BiocParallel bpmapply MulticoreParam
#' @importFrom motifStack plotMotifLogo pcm2pfm 
#' @import Rcpp data.table BSgenome
#' @seealso atSNP vignette for more information
NULL

