## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"---------------------------------
  BiocStyle::latex()

## ----eval=FALSE, echo=TRUE,results="hide"-----------------------------------------------
#  library(devtools)
#  install_github("chandlerzuo/atSNP")

## ----include=TRUE,echo=TRUE,eval=FALSE,results="markup"---------------------------------
#    source("http://bioconductor.org/biocLite.R")
#    biocLite("BSgenome.Hsapiens.UCSC.hg19")
#    biocLite("SNPlocs.Hsapiens.dbSNP.20120608")

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
data(encode_library)
length(encode_motif)
encode_motif[1]

## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------
encode_motif[[1]]
GetIUPACSequence(encode_motif[[1]])

## ----eval=TRUE, echo=TRUE, results = "markup",tidy=TRUE---------------------------------
length(encode_motifinfo)
head(encode_motifinfo)

## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------
encode_motifinfo[names(encode_motif[1])]

## ----eval=TRUE, echo = TRUE, results = "markup", tidy = TRUE----------------------------
data(jaspar_library)
jaspar_motif[[1]]
jaspar_motifinfo[names(jaspar_motif[1])]

## ----eval=FALSE, echo=TRUE, results="hide"----------------------------------------------
#  ## Source: http://meme.nbcr.net/meme/doc/examples/sample-dna-motif.meme-io
#  pwms <- LoadMotifLibrary(
#   "http://pages.stat.wisc.edu/~keles/atSNP-Data/sample-dna-motif.meme-io.txt")
#  
#  ## Source: http://compbio.mit.edu/encode-motifs/motifs.txt
#  pwms <- LoadMotifLibrary(
#   "http://pages.stat.wisc.edu/~keles/atSNP-Data/motifs.txt",
#   tag = ">", transpose = FALSE, field = 1,
#   sep = c("\t", " ", ">"), skipcols = 1,
#   skiprows = 1, pseudocount = 0)
#  
#  ## Source: http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt
#  pwms <- LoadMotifLibrary(
#   "http://pages.stat.wisc.edu/~keles/atSNP-Data/JASPAR_motifs_H_sapiens.txt",
#   tag = "/NAME",skiprows = 1, skipcols = 0, transpose = FALSE,
#   field = 2)
#  
#  ## Source: http://jaspar.genereg.net/html/DOWNLOAD/ARCHIVE/JASPAR2010/all_data/matrix_only/matrix.txt
#  pwms <- LoadMotifLibrary(
#   "http://pages.stat.wisc.edu/~keles/atSNP-Data/matrix.txt",
#   tag = ">", skiprows = 1, skipcols = 1, transpose = TRUE,
#   field = 1, sep = c("\t", " ", "\\[", "\\]", ">"),
#   pseudocount = 1)
#  
#  ## Source: http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt
#  pwms <- LoadMotifLibrary(
#   "http://pages.stat.wisc.edu/~keles/atSNP-Data/pfm_vertebrates.txt",
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
                        half.window.size = 30, default.par = TRUE, mutation = FALSE)
ncol(snp_info$sequence) == nrow(snp_tbl)
snp_info$rsid.rm


## ----eval=TRUE, echo=TRUE, results="markup",tidy=FALSE----------------------------------

  mutation_info <- LoadSNPData("test_snp_file.txt", genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                               half.window.size = 30, default.par = TRUE, mutation = TRUE)
  ncol(mutation_info$sequence) == nrow(snp_tbl)
  file.remove("test_snp_file.txt")


## ----eval=TRUE, echo=TRUE, results="markup", tidy=FALSE---------------------------------

snp_info1 <- LoadSNPData(snpids = c("rs5050", "rs616488", "rs11249433",
                           "rs182799", "rs12565013", "rs11208590"),
                         genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
  			 snp.lib = "SNPlocs.Hsapiens.dbSNP.20120608",
 			 half.window.size = 30,
			 default.par = TRUE,
			 mutation = FALSE)

## ----eval=TRUE, echo=TRUE, results="markup", tidy=TRUE----------------------------------
snp_info1$rsid.missing
snp_info1$rsid.duplicate
snp_info1$rsid.rm

## ----eval=TRUE, echo = TRUE, results="markup",tidy=FALSE--------------------------------
snp_info2 <- LoadFastaData("http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_1.fasta",
                           "http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_2.fasta",
                           default.par = TRUE)

## ----eval=TRUE, echo=TRUE, results="markup",tidy=TRUE-----------------------------------

  data(example)
  names(motif_library)
  str(snpInfo)
## to look at the motif information
  data(encode_library)
  encode_motifinfo[names(motif_library)]


## ----eval=TRUE, echo=TRUE, results="markup"---------------------------------------------
  atsnp.scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 2)
  head(atsnp.scores$snp.tbl)
  head(atsnp.scores$motif.scores[, list(snpid, motif, log_lik_ref,
                                log_lik_snp, log_lik_ratio)])

## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------
  
  atsnp.result <- ComputePValues(motif.lib = motif_library, snp.info = snpInfo,
                                 motif.scores = atsnp.scores$motif.scores,
				 ncores = 2)
  head(atsnp.result[, list(snpid, motif, pval_ref, pval_snp, pval_rank, pval_diff)])


## ----eval=TRUE, echo = TRUE, results="markup"-------------------------------------------
head(atsnp.result[order(pval_rank), list(snpid, motif, pval_ref, pval_snp, pval_rank)])

## ----eval=TRUE, echo = TRUE, results = "markup"-----------------------------------------
atsnp.result[pval_rank <= 0.1, list(snpid, motif, pval_ref, pval_snp, pval_rank)]

## ----eval=TRUE, echo = FALSE, results="hide"--------------------------------------------
atsnp.result[, pval_rank_bh := p.adjust(pval_rank, method = "BH")]

## ----eval=TRUE, echo = FALSE, results="markup"------------------------------------------
atsnp.result[, list(snpid, motif, pval_rank, pval_rank_bh)]

## ----eval=FALSE, echo =TRUE,results="markup"--------------------------------------------
#  library(qvalue)
#  atsnp.result[, qval_rank := qvalue(pval_rank)$qvalues]

## ----eval=FALSE, echo =TRUE,results="markup"--------------------------------------------
#  atsnp.result[, pval_rank_bh := p.adjust(pval_rank, method = "BH"), by = motif]
#  atsnp.result[, qval_rank := qvalue(pval_rank, pi0=1)$qvalues, by = motif]

## ----eval=FALSE, echo =TRUE,results="markup"--------------------------------------------
#  atsnp.result[, pval_rank_bh := p.adjust(pval_rank, method = "BH"), by = snpid]
#  atsnp.result[, qval_rank := qvalue(pval_rank, pi0=1)$qvalues, by = snpid]

## ----eval=TRUE,echo=TRUE,results="markup"-----------------------------------------------
  
  match_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl,
                                 motif.scores = atsnp.result,
                                 motif.lib = motif_library,
                                 snpids = c("rs10910078", "rs4486391"),
                                 motifs = names(motif_library)[1:2],
                                 ncores = 2)
  match_result[, list(snpid, motif, IUPAC, ref_match_seq, snp_match_seq)]


## ----include=TRUE,eval=TRUE, echo=TRUE,fig.align="center",dpi=600,fig.width=6,fig.height=6----

  plotMotifMatch(snp.tbl = atsnp.scores$snp.tbl,
               motif.scores = atsnp.scores$motif.scores,
               snpid = atsnp.scores$snp.tbl$snpid[1],
                motif = atsnp.scores$motif.scores$motif[1],
		motif.lib = motif_library)

## ----eval=TRUE,echo=FALSE,results="markup",cache=FALSE----------------------------------
print(sessionInfo())

