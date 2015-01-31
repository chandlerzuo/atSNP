## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
  BiocStyle::latex()

## ----eval=FALSE, echo=TRUE,results="hide"-----------------------------------------------
#  library(devtools)
#  install_github("chandlerzuo/atSNP")

## ----include=TRUE,echo=TRUE,eval=FALSE,results="markup"---------------------------------
#    source("http://bioconductor.org/biocLite.R")
#    biocLite("BSgenome.Hsapiens.UCSC.hg19")

## ----include=FALSE,eval=TRUE, echo=FALSE, results="hide"--------------------------------
source("http://bioconductor.org/biocLite.R")
if (!require("BSgenome.Hsapiens.UCSC.hg19",character.only = TRUE))
{
  biocLite("BSgenome.Hsapiens.UCSC.hg19",suppressAutoUpdate=TRUE)
  if(!require("BSgenome.Hsapiens.UCSC.hg19",character.only = TRUE)) stop("Package not found")
}
tidy.opt = list(width.cutoff = 60)

## ----eval=TRUE,echo=FALSE,results="hide"------------------------------------------------
  library(IRanges)
  library(BSgenome)

## ----eval=TRUE, echo=TRUE, results = "markup"-------------------------------------------
library(atSNP)

## ----eval=TRUE, echo=TRUE, results = "markup"-------------------------------------------
data(encode_motif)
length(motif_encode)
motif_encode[seq(3)]

## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------
motif_encode[[1]]
GetIUPACSequence(motif_encode[[1]])

## ----eval=TRUE, echo=TRUE, results = "markup",tidy=TRUE---------------------------------
length(motif_info)
head(motif_info)

## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------
motif_info[names(motif_encode[1])]

## ----eval=FALSE, echo=TRUE, results="hide"----------------------------------------------
#  
#  pwms <- LoadMotifLibrary(
#   "http://meme.nbcr.net/meme/examples/sample-dna-motif.meme-io")
#  pwms <- LoadMotifLibrary(
#   "http://compbio.mit.edu/encode-motifs/motifs.txt",
#   tag = ">", transpose = FALSE, field = 1,
#   sep = c("\t", " ", ">"), skipcols = 1,
#   skiprows = 1, pseudocount = 0)
#  pwms <- LoadMotifLibrary(
#   "http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt",
#   tag = "/NAME",skiprows = 1, skipcols = 0, transpose = FALSE,
#   field = 2)
#  pwms <- LoadMotifLibrary(
#   "http://jaspar.genereg.net/html/DOWNLOAD/ARCHIVE/JASPAR2010/all_data/matrix_only/matrix.txt",
#   tag = ">", skiprows = 1, skipcols = 1, transpose = TRUE,
#   field = 1, sep = c("\t", " ", "\\[", "\\]", ">"),
#   pseudocount = 1)
#  pwms <- LoadMotifLibrary(
#   "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt",
#   tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1,
#   sep = c(">", "\t", " "), pseudocount = 1)
#  
#  ## pwms <- LoadMotifLibrary(
#  ##  "http://gibbs.biomed.ucf.edu/PreDREM/download/nonredundantmotif.transfac",
#  ##  tag = "DE", skiprows = 1, skipcols = 1,
#  ##  transpose = FALSE, field = 2, sep = "\t")
#  

## ----eval=TRUE, echo=TRUE, results="markup",tidy=FALSE----------------------------------

  data(example)
  write.table(snp_tbl, file = "test_snp_file.txt",
            row.names = FALSE, quote = FALSE)
  snp_info <- LoadSNPData("test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                          half.window.size = 30, default.par = TRUE,mutation = FALSE)
  ncol(snp_info$sequence) == nrow(snp_tbl)


## ----eval=TRUE, echo=TRUE, results="markup",tidy=FALSE----------------------------------

  mutation_info <- LoadSNPData("test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                               half.window.size = 30, default.par = TRUE, mutation = TRUE)
  ncol(mutation_info$sequence) == nrow(snp_tbl)
  file.remove("test_snp_file.txt")


## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------

  data(example)
  names(motif_library)
  str(snpInfo)
## to look at the motif information
  data(encode_motif)
  motif_info[names(motif_library)]


## ----eval=TRUE, echo=TRUE, results="markup"---------------------------------------------

  motif_score <- ComputeMotifScore(motif_library, snpInfo, ncores = 2)
  motif_score$snp.tbl
  motif_score$motif.scores[, list(snpid, motif, log_lik_ref,
                                log_lik_snp, log_lik_ratio)]

## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------
  
  motif.scores <- ComputePValues(motif.lib = motif_library, snp.info = snpInfo,
                                 motif.scores = motif_scores$motif.scores,
				 ncores = 2)
  motif.scores[, list(snpid, motif, pval_ref, pval_snp, pval_rank, pval_diff)]


## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------
  
  match_result <- MatchSubsequence(snp.tbl = motif_scores$snp.tbl,
                                 motif.scores = motif.scores,
                                 motif.lib = motif_library,
                                 snpids = c("rs10910078", "rs4486391"),
                                 motifs = names(motif_library)[1:2],
                                 ncores = 2)
  match_result[, list(snpid, motif, IUPAC, ref_match_seq, snp_match_seq)]


## ----include=TRUE,eval=TRUE, echo=TRUE,fig.align="center",dpi=600,fig.width=6,fig.height=6----

  plotMotifMatch(snp.tbl = motif_scores$snp.tbl,
               motif.scores = motif_scores$motif.scores,
               snpid = motif_scores$snp.tbl$snpid[1],
               motif.lib = motif_library,
               motif = motif_scores$motif.scores$motif[1])


## ----eval=TRUE,echo=FALSE,results="markup",cache=FALSE----------------------------------
print(sessionInfo())

