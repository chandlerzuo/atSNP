#' @name LoadMotifLibrary
#' @title Load position weight matrices.
#' @description Load the file for position weight matrices for motifs.
#' @param filename a MEME format file name.
#' @param urlname URL containing a MEME format file.
#' @param tag A string that marks the description line of the position weight 
#' matrix.
#' @param skiprows Number of description lines before each position weight 
#' matrix.
#' @param skipcols Number of columns to be skipped in the position weight 
#' matrix.
#' @param transpose If TRUE (default), then the position weight matrix should 
#' have 4 columns. Otherwise, it should have 4 rows.
#' @param field The index of the field in the description line, seperated by 
#' space, that indicates the motif name.
#' @param sep A vector of chars for the string separators to parse each lines of
#'  the matrix. Default: c(" ", "\\t").
#' @param pseudocount An integer for the pseudocount added to each of the 
#' original matrices. Default: 0. Recommended to be 1 if the original matrices 
#' are position frequency matrices.
#' @details This function reads the formatted file containing motif information 
#' and convert them into a list of position weight matrices. The list of 
#' arguments should provide enough flexibility of importing a varying number of 
#' formats. Some examples are the following:
#' For MEME format, the suggested arguments are: tag = 'Motif', skiprows = 2, 
#' skipcols = 0, transpose = FALSE, field = 2, sep = ' ';
#' For motif files from JOHNSON lab (i.e. 
#' http://johnsonlab.ucsf.edu/mochi_files/JASPAR_motifs_H_sapiens.txt), 
#' the suggested arguments are: tag = '/NAME', skiprows = 1, skipcols = 0, 
#' transpose = FALSE, field = 2, sep = "\\t";
#' For JASPAR pfm matrices (i.e. http://jaspar.genereg.net/download/CORE/JASPAR
#' 2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt), the suggested arguments
#' are: tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1, 
#' sep = "\\t"; For the TRANSFAC library provided by UCF bioinformatics groups 
#' (i.e. http://gibbs.biomed.ucf.edu/PreDREM/download/nonredundantmotif.transfac
#' ), the suggested arguments are: tag = "DE", skiprows = 1, skipcols = 1, 
#' transpose = FALSE, field = 2, sep = "\\t".
#' @return A list object of position weight matrices.
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo 
#' \email{chandler.c.zuo@@gmail.com}
#' @examples
#' pwms <- LoadMotifLibrary(
#' urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/pfm_vertebrates.txt", 
#' tag = ">", transpose = FALSE, field = 1, sep = c("\t", " ", ">"), 
#' skipcols = 1, skiprows = 1, pseudocount = 1)
#' @useDynLib atSNP
#' @import BiocFileCache
#' @import rappdirs
#' @export
LoadMotifLibrary <- function(filename = NULL, urlname = NULL, tag = "MOTIF", transpose = FALSE, field = 2, sep = c("\t", " "), skipcols = 0, skiprows = 2, pseudocount = 0) {
  if ( is.null(filename) & is.null(urlname) )  {
    stop("one argument among 'filename' and 'urlname' should be provided.")
  } else {
    if(!is.null(filename) & !is.null(urlname))
      stop("only one argument among 'filename' and 'urlname' should be provided.")
  }
  if(is.null(filename)) {
    bfc <- BiocFileCache(cache = user_cache_dir(appname = "BiocFileCache"), ask = FALSE)
    rid <- bfcrid(bfcquery(bfc, query=basename(urlname), exact=TRUE, field="rname"))
    if (!length(rid))
      rid <- names(bfcadd(bfc, rname=basename(urlname), urlname))
    lines<-readLines(bfcrpath(rids=rid))
  } else {
    lines<-readLines(filename)
  }

    motifLineNums <- grep(tag, lines)
  if(length(myStrSplit(lines[motifLineNums[1]], sep)[[1]]) >= field) {
    motifnames <-
      sapply(myStrSplit(lines[motifLineNums], sep), function(x) x[field])
  } else {
    motifnames <-
      sapply(myStrSplit(lines[motifLineNums], sep), function(x) x[field])
  }
  allmats <- as.list(seq_along(motifnames))

  for(matrixId in seq_along(motifLineNums)) {
    motifLineNum <- motifLineNums[matrixId] + skiprows
    if(!transpose) {
      pwm <- NULL
      nrows <- 0
      tmp <- myStrSplit(lines[nrows + motifLineNum], split = sep)[[1]]
      tmp <- tmp[nchar(tmp) > 0]
      while(length(tmp) >= 4 + skipcols) {
        tmp <- as.numeric(tmp[skipcols + seq(4)])
        if(sum(is.na(tmp)) == 0) {
          pwm <- rbind(pwm, tmp)
        }
        nrows <- nrows + 1
        tmp <- myStrSplit(lines[nrows + motifLineNum], split = sep)[[1]]
        tmp <- tmp[nchar(tmp) > 0]
      }
    } else {
      nrows <- 4
      if(skipcols == 0) {
        pwm <-
          matrix(as.numeric(unlist(myStrSplit(lines[seq(nrows) + motifLineNum - 1], split = sep))), ncol = 4)
      } else {
        pwm <-
          matrix(as.numeric(unlist(sapply(myStrSplit(lines[seq(nrows) + motifLineNum - 1], split = sep), function(x) x[-seq(skipcols)]))), ncol = 4)
      }
    }
    pwm <- pwm + pseudocount
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
#' @title Load the SNP information and code the genome sequences around the SNP 
#' locations.
#' @description Load the SNP data.
#' @param filename A table containing the SNP information. Must contain at least
#'  five columns with exactly the following names:
#' \tabular{ll}{
#' chr \tab chromosome.\cr
#' snp \tab The nucleotide position of the SNP.\cr
#' snpid \tab The names of the SNPs.\cr
#' a1 \tab The deoxyribose for one allele.\cr
#' a2 \tab The deoxyribose for the other allele.\cr
#' }
#' If this file exists already, it is used to extract the SNP information. 
#' Otherwise, SNP information extracted using argument 'snpids' is outputted to 
#' this file.
#' @param snpids A vector of rs ids for the SNPs. This argument is overidden 
#' if the file with name \code{filename} exists.
#' @param snp.lib A string of the library name to obtain the SNP information 
#' based on rs ids. Default: "SNPlocs.Hsapiens.dbSNP144.GRCh38".
#' @param genome.lib A string of the library name for the genome version. 
#' Default: "BSgenome.Hsapiens.UCSC.hg38".
#' @param half.window.size An integer for the half window size around the SNP 
#' within which the motifs are matched. Default: 30.
#' @param default.par A boolean for whether using the default Markov parameters.
#'  Default: FALSE.
#' @param mutation A boolean for whether this is mutation data. See details for 
#' more information. Default: FALSE.
#' @param ... Other parameters passed to \code{\link[utils]{read.table}}.
#' @details This function extracts the nucleotide sequence within a window 
#' around each SNP and code them using 1-A, 2-C, 3-G, 4-T.\cr
#' There are two ways of obtaining the nucleotide sequences. If \code{filename} 
#' is not NULL and the file exists, it should contain the positions and alleles 
#' for each SNP. Based on such information, the sequences around SNP positions 
#' are extracted using the Bioconductor annotation package specified by 
#' \code{genome.lib}. Users should make sure that this annotation package 
#' corresponds to the correct species and genome version of the actual data. 
#' Alternatively, users can also provide a vector of rs ids via the argument 
#' \code{snpids}. The SNP locations and allele information is then obtained via 
#' the Bioconductor annotation package specified by \code{snp.lib}, and passed 
#' on to the package specified by \code{genome.lib} to further obtain the 
#' nucleotide sequences.\cr
#' If \code{mutation=FALSE} (default), this function assumes that the data is 
#' for SNP analysis, and the reference genome should be consistent with either 
#' the a1 or a2 nucleotide. When extracting the genome sequence around each SNP 
#' position, this function compares the nucleotide at the SNP location on the 
#' reference genome with both a1 and a2 to distinguish between the reference 
#' allele and the SNP allele. If the nucleotide extracted from the reference 
#' genome does not match either a1 or a2, the SNP is discarded. The discarded 
#' SNPs are in the 'rsid.rm' field in the output.\cr
#' Alternatively, if \code{mutation=TRUE}, this function assumes that the data 
#' is for general single nucleotide mutation analysis. After extracting the 
#' genome sequence around each SNP position, it replaces the nucleotide at the 
#' SNP location by the a1 nucleotide as the 'reference' allele sequence, and by 
#' the a2 nucleotide as the 'snp' allele sequence. It does NOT discard the 
#' sequence even if neither a1 or a2 matches the reference genome. When this 
#' data set is used in other functions, such as \code{\link{ComputeMotifScore}},
#'  \code{\link{ComputePValues}}, all the results (i.e. affinity scores and 
#'  their p-values) for the reference allele are indeed for the a1 allele, and 
#'  results for the SNP allele are indeed for the a2 allele.\cr
#' If the input is a list of rsid's, the SNP information extracted from 
#' \code{snp.lib} may contain more than two alleles for a single location. For 
#' such cases, \code{\link{LoadSNPData}} first extracts all pairs of alleles 
#' associated with those locations. If 'mutation=TRUE', all those pairs are 
#' considered as pairs of reference and SNP alleles, and their information is 
#' contained in 'sequence_matrix', 'a1', 'a2' and 'snpid'. If 'mutation=FALSE', 
#' \code{\link{LoadSNPData}} further filters these pairs based on whether one 
#' allele matches to the reference genome nucleotide extracted from 
#' \code{genome.lib}. Only those pairs with one allele matching the reference 
#' genome nucleotide is considered as pairs of reference and SNP alleles, with 
#' their information contained in 'sequence_matrix', 'a1', 'a2' and 'snpid'.\cr
#' @return A list object containing the following components:
#' \tabular{ll}{
#' sequence_matrix \tab A list of integer vectors representing the deroxyribose 
#' sequence around each SNP.\cr
#' a1 \tab An integer vector for the deroxyribose at the SNP location on the 
#' reference genome.\cr
#' a2 \tab An integer vector for the deroxyribose at the SNP location on the 
#' SNP genome.\cr
#' snpid \tab A string vector for the SNP rsids.\cr
#' rsid.missing \tab If the data source is a list of rsids, this field records 
#' rsids for SNPs that are discarded because they are not in the SNPlocs package.\cr
#' rsid.duplicate \tab If the data source is a list of rsids, this field records
#'  rsids for SNPs that based on the SNPlocs package, this locus has more than 
#'  2 alleles. \cr
#' rsid.na \tab This field records rsids for SNPs that are discarded because the
#'  nucleotide sequences contain none ACGT characters.\cr
#' rsid.rm \tab If the data source is a table and \code{mutation=FALSE}, this 
#' field records rsids for SNPs that are discarded because the nucleotide on the
#'  reference genome matches neither 'a1' or 'a2' in the data source.\cr
#' }
#' The results are coded as: "A"-1, "C"-2, "G"-3, "T"-4.
#' @author Chandler Zuo \email{chandler.c.zuo@@gmail.com}
#' @examples
#' \dontrun{LoadSNPData(snpids = c("rs53576", "rs7412"), 
#' genome.lib ="BSgenome.Hsapiens.UCSC.hg38", snp.lib = 
#' "SNPlocs.Hsapiens.dbSNP144.GRCh38", half.window.size = 30, default.par = TRUE
#' , mutation = FALSE)}
#' @import BSgenome
#' @useDynLib atSNP
#' @export
LoadSNPData <- function(filename = NULL, genome.lib = "BSgenome.Hsapiens.UCSC.hg38",
                        snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38",
                        snpids = NULL, half.window.size = 30, default.par = FALSE,
                        mutation = FALSE, ...) {
  IUPAC_CODE_MAP = NULL
  useFile <- FALSE
  rsid.rm <- rsid.missing <- rsid.duplicate <- rsid.na <- NULL
  if(!is.null(filename)) {
    if(file.exists(filename)) {
      useFile <- TRUE
    }
  }
  if(useFile) {
    if(!is.null(snpids)) {
      message("Warning: load SNP information from 'filename' only. The argument 'snpids' is overridden.")
    }
    tbl <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, ...)
    tbl<-tbl[c("snpid", "chr", "snp", "a1", "a2")]
    ## check if the input file has the required information
    if(sum(!c("snp", "chr", "a1", "a2", "snpid") %in% names(tbl)) > 0) {
      stop("'filename' must be a table containing 'snp' and 'chr' columns.")
    }
    snpid.index <- seq(nrow(tbl))
  } else {
    if(is.null(snpids)) {
      stop("either 'snpids' should be a vector, or 'filename' should be the file name that contains the SNP information.")
    }
    snpid.index <- seq_along(snpids)
    ## load the corresponding snp library
    library(package = snp.lib, character.only = TRUE)
    rsid.missing.all <- NULL
    snps <- get(snp.lib)
    snp.loc <- tryCatch({snpsById(snps, snpids, ifnotfound="error")}, error = function(e) return(e$message))
    ## remove rsids not included in the database
    if(is(snp.loc, "character")) {
      rsid.missing <- myStrSplit(snp.loc, split = c(": ", "\n"))[[1]][-1]
      rsid.missing.all <- myStrSplit(rsid.missing, split = c(",", " "))[[1]]
      snpids <- snpids[!snpids %in% rsid.missing]
      snp.loc <- snpsById(snps, snpids, ifnotfound="drop")
    } 
    
    if(!is.null(rsid.missing.all)) {
      message("Warning: the following rsids are not included in the database and discarded: ")
      message(paste(rsid.missing.all, collapse = ", "))
      rsid.missing <- rsid.missing.all
    }
    
    snp.alleles <- IUPAC_CODE_MAP[snp.loc@elementMetadata@listData$alleles_as_ambig]
    snp.strands<-as.character(snp.loc@strand)
    if(sum(nchar(snp.alleles) > 2) > 0) {
      message("Warning: the following SNPs have more than 2 alleles. All pairs of nucleotides are considered as pairs of the SNP and the reference allele:")
      rsid.duplicate <- snpids[nchar(snp.alleles) > 2]
      message(paste(rsid.duplicate, collapse = ", "))
    }
    ## retain only SNPs with >= 2 alleles
    tbl <- NULL
    for(nalleles in 2:4) {
      ids <- which(sapply(snp.alleles, nchar) == nalleles)
      if(length(ids) == 0) {
        next
      }
      snp.loc.n <- snp.loc[ids]
      snp.alleles.n <- snp.alleles[ids]
      snp.ids.n <- snpids[ids]
      snp.alleles.n <- strsplit(snp.alleles.n, "")
      snp.strands.n <- snp.strands[ids]
      ## get all pairs of alleles
      for(i_allele1 in seq(nalleles - 1)) {
        for(i_allele2 in (i_allele1 + 1):nalleles) {
          a1 <- sapply(snp.alleles.n, function(x) x[i_allele1])
          a2 <- sapply(snp.alleles.n, function(x) x[i_allele2])
          
          ## revert the alleles on the reverse strand
          id.rev <- which(snp.strands.n == "-")
          if(length(id.rev) > 0) {
            rev.codes <- c("A", "C", "G", "T")
            names(rev.codes) <- rev(rev.codes)
            a1[id.rev] <- rev.codes[a1[id.rev]]
            a2[id.rev] <- rev.codes[a2[id.rev]]
          }
          tbl <- rbind(tbl,
                       data.frame(snp = as.numeric(snp.loc.n@ranges), 
                                  chr = as.character(paste0("chr", gsub("ch", "", snp.loc.n@seqnames))), 			  
                                  a1 = as.character(a1),
                                  a2 = as.character(a2),
                                  snpid = as.character(snp.ids.n),
                                  index = snpid.index[ids])
          )
        }
      }
    }
    
    tbl <- tbl[order(tbl$index), ]
    if(!is.null(filename)) {
      write.table(tbl[, c('snp', 'chr', 'a1', 'a2', 'snpid')], file = filename, row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    snpid.index <- tbl$index
  }
  
  ## load the corresponding genome version
  library(package = genome.lib, character.only = TRUE)
  species<-get(strsplit(genome.lib, "[.]")[[1]][2])
  seqvec <- getSeq(species,
                   as.character(tbl$chr),
                   start=tbl$snp - half.window.size,
                   end=tbl$snp + half.window.size,
                   as.character=TRUE)
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  sequences <- sapply(seqvec, function(x) codes[strsplit(x, "")[[1]]])
  sequences <- matrix(sequences, ncol = length(seqvec))
  snpid.output <- as.character(tbl$snpid)
  rownames(sequences) <- NULL
  a1 <- codes[as.character(tbl$a1)]
  a2 <- codes[as.character(tbl$a2)]
  names(a1) <- names(a2) <- NULL
  keep.id <- which(colSums(is.na(sequences)) == 0)
  if(length(keep.id) < nrow(tbl)) {
    message("Warning: the following rows are discarded because the reference genome sequences contain non ACGT characters:")
    rsid.na <- tbl[-keep.id, ]$snpid
    print(tbl[-keep.id, ])
  }
  ## remove sequences containing non ACGT characters
  sequences <- sequences[, keep.id, drop = FALSE]
  a1 <- a1[keep.id]
  a2 <- a2[keep.id]
  snpid.output <- snpid.output[keep.id]
  snpid.index <- snpid.index[keep.id]
  ## whether use the default parameters
  if(!default.par) {
    transition <- .Call("transition_matrix", sequences, package = "atSNP")
    prior <- rowSums(transition)
    prior <- prior / sum(prior)
    transition <- transition / rowSums(transition)
    names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  } else {
    data(default_par, envir = environment())
  }
  if(!mutation) {
    ## real SNP data
    a1.ref.base.id <- which(a1 == sequences[half.window.size + 1, ])
    a2.ref.base.id <- which(a2 == sequences[half.window.size + 1, ])
    ## store SNPs that have the same base in SNP and REF alleles only once
    a2.ref.base.id <- a2.ref.base.id[!a2.ref.base.id %in% a1.ref.base.id]
    discard.id <- setdiff(seq_along(a1), c(a1.ref.base.id, a2.ref.base.id))
    if(length(discard.id) > 0) {
      message("Warning: the following sequences are discarded because the reference nucleotide matches to neither a1 nor a2:")
      rsid.rm <- unique(as.character(tbl[keep.id[discard.id], ]$snpid))
      message("snpid\tchr\tsnp\ta1\ta2")
      message(paste(apply(tbl[keep.id[discard.id], c("snpid", "chr", "snp", "a1", "a2")], 1, function(x) paste(x, collapse = "\t")), collapse = "\n"))
    }
  } else {
    ## single nucleotide mutation data
    a1.ref.base.id <- seq_along(a1)
    a2.ref.base.id <- numeric(0)
  }
  sequences <- sequences[, c(a1.ref.base.id, a2.ref.base.id), drop = FALSE]
  snpid.output <- snpid.output[c(a1.ref.base.id, a2.ref.base.id)]
  ref.base <- c(a1[a1.ref.base.id], a2[a2.ref.base.id])
  snp.base <- c(a2[a1.ref.base.id], a1[a2.ref.base.id])
  snpid.index <- snpid.index[c(a1.ref.base.id, a2.ref.base.id)]
  ## Keep the order of SNPs as in the input file
  if(useFile) {
    output.index = seq(ncol(sequences))
  } else {
    output.index = order(snpid.index)
  }
  return(list(
    sequence_matrix= matrix(sequences[, output.index], nrow=2*half.window.size+1),
    ref_base = ref.base[output.index],
    snp_base = snp.base[output.index],
    snpids = snpid.output[output.index],
    transition = transition,
    prior = prior,
    rsid.na = rsid.na,
    rsid.rm = rsid.rm,
    rsid.duplicate = rsid.duplicate,
    rsid.missing = rsid.missing
  ))
}

#' @name LoadFastaData
#' @title Load the SNP data from fasta files.
#' @description Load SNP data.
#' @param ref.filename a fastq file name for the reference allele sequences.
#' @param snp.filename a fastq file name for the SNP allele sequences.
#'@param ref.urlname URL of a fastq file for the reference allele sequences.
#' @param snp.urlname URL of a fastq file for the SNP allele sequences.
#' @param snpids SNP IDs			  
#' @param default.par A boolean for whether using the default Markov parameters.
#'  Default: FALSE.
#' @return A list object containing the following components:
#' \tabular{ll}{
#' sequence_matrix \tab A list of integer vectors representing the deroxyribose 
#' sequence around each SNP.\cr
#' a1 \tab An integer vector for the deroxyribose at the SNP location on the 
#' reference genome.\cr
#' a2 \tab An integer vector for the deroxyribose at the SNP location on the SNP
#'  genome.\cr
#' }
#' The results are coded as: "A"-1, "C"-2, "G"-3, "T"-4.
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo 
#' \email{chandler.c.zuo@@gmail.com}
#' @examples LoadFastaData(
#' ref.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_1.fasta",
#' snp.urlname="http://pages.stat.wisc.edu/~keles/atSNP-Data/sample_2.fasta")
#' @useDynLib atSNP
#' @import BiocFileCache
#' @import rappdirs
#' @export
LoadFastaData <- function(ref.filename = NULL, snp.filename = NULL, 
                          ref.urlname = NULL, snp.urlname = NULL,
                          snpids = NULL, default.par = FALSE) {
  if ( is.null(ref.filename) & is.null(ref.urlname) )  {
    stop("one argument among 'ref.filename' and 'ref.urlname' should be provided.")
  } else {
    if(!is.null(ref.filename) & !is.null(ref.urlname))
      stop("only one argument among 'ref.filename' and 'ref.urlname' should be provided.")
  }
  if(is.null(ref.filename)) {
    bfc <- BiocFileCache(cache = user_cache_dir(appname = "BiocFileCache"), ask = FALSE)
    ref.rid <- bfcrid(bfcquery(bfc, query=basename(ref.urlname), exact=TRUE, field="rname"))
    if (!length(ref.rid))
      ref.rid <- names(bfcadd(bfc, rname=basename(ref.urlname), ref.urlname))
    refdat <- read.table(bfcrpath(rids=ref.rid))
  } else {
    refdat <- read.table(ref.filename)
  }
  
  if ( is.null(snp.filename) & is.null(snp.urlname) )  {
    stop("one argument among 'snp.filename' and 'snp.urlname' should be provided.")
  } else {
    if(!is.null(snp.filename) & !is.null(snp.urlname))
      stop("only one argument among 'snp.filename' and 'snp.urlname' should be provided.")
  }
  if(is.null(snp.filename)) {
    bfc <- BiocFileCache(cache = user_cache_dir(appname = "BiocFileCache"), ask = FALSE)
    snp.rid <- bfcrid(bfcquery(bfc, query=basename(snp.urlname), exact=TRUE, field="rname"))
    if (!length(snp.rid))
      snp.rid <- names(bfcadd(bfc, rname=basename(snp.urlname), snp.urlname))
    snpdat <- read.table(bfcrpath(rids=snp.rid))
  } else {
    snpdat <- read.table(snp.filename)
  }

  if(nrow(refdat) != nrow(snpdat)) {
    stop("'ref.data' and 'snp.data' should have the same number of rows.")
  }
  n <- nrow(refdat)
  if(n%%2==1) {
    stop("'ref.data' and 'snp.data' should have an even number of rows.")
    }

  ids <- 2 * seq(n / 2)
  refseqs <- as.character(refdat[ids, 1])
  snpseqs <- as.character(snpdat[ids, 1])

  if(is.null(snpids)) {
   snpids<-gsub(">", "", as.character(refdat[ids-1,1]))
 }

  if(var(sapply(refseqs, nchar)) != 0) {
    stop("sequences in 'ref.data' have different lengths.")
  }
  if(var(sapply(snpseqs, nchar)) != 0) {
    stop("sequences in 'snp.data' have different lengths.")
  }
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  refmat <- sapply(refseqs,
                   function(x)
                   codes[strsplit(x, "")[[1]]])
  snpmat <- sapply(snpseqs,
                   function(x)
                   codes[strsplit(x, "")[[1]]])
  colnames(refmat) <- colnames(snpmat) <- rownames(refmat) <- rownames(snpmat) <- NULL
  m <- nrow(refmat)
  if(nrow(refmat) != nrow(snpmat)) {
    stop("the sequences for the SNP alleles and the reference alleles have different lengths.")
  }
  if(m %% 2 == 0) {
    stop("each sequence must have an odd number of length.")
  }

  id.na1 <- which(apply(refmat, 2, function(x) sum(is.na(x))) > 0)
  id.na2 <- which(is.na(snpmat[(m + 1) / 2, ]))
  id.na <- union(id.na1, id.na2)
  if(length(id.na) > 0) {
    message("Warning: the following sequences are discarded, because they include bases other than A, C, G, T: ")
    message(paste(id.na, collapse = ", "))
  }

  id.wrong <- which(apply((refmat != snpmat)[-(m + 1) / 2, ], 2, sum) > 0)
  if(length(id.wrong) > 0) {
    message("Warning: the following sequences are discarded, because they have unidentical nucleotides between the SNP and the reference allele at positions other than the central location: ")
    message(paste(id.wrong, collapse = ", "))
  }

  ids <- union(id.wrong, id.na)
  if(length(ids) > 0) {
    sequences <- refmat[, -ids, drop = FALSE]
    ref.base <- refmat[(m + 1) / 2, -ids]
    snp.base <- snpmat[(m + 1) / 2, -ids]
  } else {
    sequences <- refmat
    ref.base <- refmat[(m + 1) / 2, ]
    snp.base <- snpmat[(m + 1) / 2, ]
  }

  if(!default.par) {
    transition <- .Call("transition_matrix", sequences, package = "atSNP")
    prior <- rowSums(transition)
    prior <- prior / sum(prior)
    transition <- transition / rowSums(transition)
    names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  } else {
    data(default_par, envir = environment())
  }
  colnames(sequences)<-snpids

  return(list(
              sequence_matrix= sequences,
              ref_base = ref.base,
              snp_base = snp.base,
              transition = transition,
              prior = prior, snpids = snpids
              ))
}

#' @name ComputeMotifScore
#' @title Compute the scores for SNP effects on motifs.
#' @description Compute the log-likelihood scores for motifs.
#' @param motif.lib A list object with the output format of function 
#' \code{\link{LoadMotifLibrary}}.
#' @param snp.info A list object with the output format of function 
#' \code{\link{LoadSNPData}}.
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @details This function computes the binding affinity scores for both alleles 
#' at each SNP window. For each pair of SNP and motif, it finds the subsequence 
#' from both strand that maximizes the affinity binding score. It returns both 
#' the matching positions and the maximized affinity scores.
#' @return A list of two data.frame's. Field \code{snp.tbl} contains:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleotide sequence.\cr
#' snp_seq \tab SNP allele nucleotide sequence.\cr
#' ref_seq_rev \tab Reference allele nucleotide sequence on the reverse 
#' strand.\cr
#' snp_seq_rev \tab SNP allele nucleotide sequence on the reverse strand.\cr}
#' Field \code{motif.score} contains:
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
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo
#' \email{chandler.c.zuo@@gmail.com}
#' @examples
#' data(example)
#' ComputeMotifScore(motif_library, snpInfo, ncores = 2)
#' @useDynLib atSNP
#' @import data.table
#' @importFrom BiocParallel bpmapply MulticoreParam SnowParam
#' @export
ComputeMotifScore <- function(motif.lib, snp.info, ncores = 1) {
  ## check arguments
  if(sum(!unlist(sapply(motif.lib, is.matrix))) > 0 | sum(unlist(sapply(motif.lib, ncol)) != 4) > 0) {
    stop("'motif.lib' must be a list of numeric matrices each with 4 columns.")
  }
  if(sum(!c("sequence_matrix", "snp_base", "ref_base") %in% names(snp.info)) > 0) {
    stop("'snp.info' must contain three components: 'ref_base', 'snp_base', 'sequence_matrix'.")
  }
  if(ncol(snp.info$sequence_matrix) != length(snp.info$ref_base) | length(snp.info$ref_base) != length(snp.info$snp_base)) {
    stop("the number of columns of 'snp.info$sequence_matrix', the length of 'snp.info$ref_base' and the length of 'snp.info$snp_base' must be the same.")
  }
  if( ! all( sort(unique(c(c(snp.info$sequence_matrix), snp.info$ref_base, snp.info$snp_base))) %in% seq(4) )) {
    stop("'snp.info$sequence_matrix', 'snp.info$ref_base', 'snp.info$snp_base' can only contain entries in 1, 2, 3, 4.")
  }
  if(nrow(snp.info$sequence_matrix) / 2 == as.integer(nrow(snp.info$sequence_matrix) / 2)) {
    stop("'snp.info$sequence_matrix' must have an odd number of rows so that the central row refers to the SNP nucleotide.")
  }

  motifs <- names(motif.lib)
  snpids <- snp.info$snpids
  snpbases<-ifelse(snp.info$snp_base==1, "A", ifelse(snp.info$snp_base==2, "C", ifelse(snp.info$snp_base==3, "G", "T")))	
# nsnps <- ncol(snp.info$sequence_matrix)
  nmotifs <- length(motif.lib)
  len_seq <- nrow(snp.info$sequence_matrix)
  snp.info$sequence_matrix[(len_seq+1)/2,]<-snp.info$ref_base
  ncores <- min(c(ncores, length(snp.info$ref_base)))

#  startParallel(ncores)
  k <- as.integer(length(snp.info$ref_base) / ncores)
  if(ncores > 1) {
    if(Sys.info()[["sysname"]] == "Windows"){
      snow <- SnowParam(workers = ncores, type = "SOCK")
      motif_score_par_list <- bpmapply(function(x) motif_score_par(i=x, par.k=k, par.ncores=ncores, par.motifs=motifs, par.nmotifs=nmotifs, par.snpids=snpids,  par.snpbases=snpbases, par.len_seq=len_seq, par.motif.lib=motif.lib, par.snp.info=snp.info), seq(ncores), BPPARAM = snow,SIMPLIFY = FALSE)
    } else {
      motif_score_par_list <- bpmapply(function(x) motif_score_par(i=x, par.k=k, par.ncores=ncores, par.motifs=motifs, par.nmotifs=nmotifs, par.snpids=snpids,  par.snpbases=snpbases, par.len_seq=len_seq, par.motif.lib=motif.lib, par.snp.info=snp.info), seq(ncores), BPPARAM = MulticoreParam(workers = ncores),
                                SIMPLIFY = FALSE)
    } 
  }
    else {
    motif_score_par_list <- list(motif_score_par(i=1, par.k=k, par.ncores=ncores, par.motifs=motifs, par.nmotifs=nmotifs, 
                                                 par.snpids=snpids,  par.snpbases=snpbases, par.len_seq=len_seq, 
                                                 par.motif.lib=motif.lib, par.snp.info=snp.info))
   }
#  endParallel()

  motif.scores_dt <- motif_score_par_list[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      motif.scores_dt <- rbind(motif.scores_dt, motif_score_par_list[[i]])
    }
  }
#  setkey(motif.scores_dt, motif, snpid, snpbase)
  motif.scores<-as.data.frame(motif.scores_dt, stringAsFactors=FALSE)
  motif.scores<-motif.scores[order(motif.scores$motif, motif.scores$snpid, motif.scores$snpbase), ]
  # motif.scores[with(motif.scores, order(motif, snpid, snpbase)), ]
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

#  snp_tbl <- data.table(snpid = snpids,
#                        ref_seq = ref_seqs,
#                        snp_seq = snp_seqs,
#                        ref_seq_rev = ref_seqs_rev,
#                        snp_seq_rev = snp_seqs_rev,
#                        snpbase=snpbases)
#  setkey(snp_tbl, snpid, snpbase)
  snp_tbl <- data.frame(snpid = snpids, ref_seq = ref_seqs, snp_seq = snp_seqs, ref_seq_rev = ref_seqs_rev, snp_seq_rev = snp_seqs_rev, snpbase = snpbases, stringsAsFactors=FALSE)
  snp_tbl[with(snp_tbl, order(snpid, snpbase)),]
  return(list(snp.tbl = snp_tbl, motif.scores = motif.scores))
}


#' @name MatchSubsequence
#' @title Compute the matching subsequence.
#' @description This function combines the SNP set, the motif library and the 
#' affinity score table and produce the matching subsequence found at each SNP 
#' location for each motif.
#' @param snp.tbl A data.frame with the following information:
#' \tabular{cc}{
#' snpid \tab SNP id.\cr
#' ref_seq \tab Reference allele nucleotide sequence.\cr
#' snp_seq \tab SNP allele nucleotide sequence.\cr
#' ref_seq_rev \tab Reference allele nucleotide sequence on the reverse 
#' strand.\cr
#' snp_seq_rev \tab SNP allele nucleotide sequence on the reverse strand.\cr}
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
#' @param motif.lib A list of the position weight matrices for the motifs.
#' @param snpids A subset of snpids to compute the subsequences. Default: NULL, 
#' when all snps are computed.
#' @param motifs A subset of motifs to compute the subsequences. Default: NULL, 
#' when all motifs are computed.
#' @param ncores The number of cores used for parallel computing.
#' @return A data.frame containing all columns in both \code{snp.tbl} and 
#' \code{motif.scores}. In addition, the following columns are added:
#' \tabular{ll}{
#' ref_match_seq \tab Best matching subsequence on the reference allele.\cr
#' snp_match_seq \tab Best matching subsequence on the SNP allele.\cr
#' ref_seq_snp_match \tab Subsequence on the reference allele corresponding to 
#' the best matching location on the SNP allele.\cr
#' snp_seq_ref_match \tab Subsequence on the SNP allele corresponding to the 
#' best matching location on the reference allele.\cr
#' }
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo 
#' \email{chandler.c.zuo@@gmail.com}
#' @examples
#' data(example)
#' MatchSubsequence(motif_scores$snp.tbl, motif_scores$motif.scores, 
#' motif_library, ncores=2)
#' @useDynLib atSNP
#' @import data.table
#' @importFrom BiocParallel bpmapply MulticoreParam SnowParam
#' @export
MatchSubsequence <- function(snp.tbl, motif.scores, motif.lib, snpids = NULL, motifs = NULL, ncores = 1) {
  motif = NULL
  snpid = NULL
  snpbase = NULL
  len_seq = NULL
  ref_seq = NULL
  motif = NULL
  if(is.null(snpids)) {
    snpids <- unique(snp.tbl$snpid)
  }
  if(is.null(motifs)) {
    motifs <- unique(motif.scores$motif)
  }
  if(sum(! motifs %in% names(motif.lib)) > 0) {
    motif.discard <- setdiff(motifs, unique(names(motif.lib)))
    message("Warning: the following motifs are not included in 'motif.lib' and are discarded: ")
    message(paste(motif.discard, collapse = ", "))
    motifs <- setdiff(motifs, motif.discard)
  }
  if(sum(! snpids %in% motif.scores$snpid) > 0) {
    snp.discard <- setdiff(snpids, unique(motif.scores$snpid))
    message("Warning: the following snpids are not included in 'motif.scores' and are discarded: ")
    message(paste(snp.discard, collapse = ", "))
    snpids <- setdiff(snpids, snp.discard)
  }
  snpids <- unique(snpids)
  motifs <- unique(motifs)
  motif.scores_dt<-as.data.table(motif.scores)
  setkey(motif.scores_dt, motif, snpid, snpbase)
  snp.tbl <- as.data.table(snp.tbl)
  setkey(snp.tbl, snpid, snpbase)
  motif.scores <- motif.scores_dt[snpid %in% snpids & motif %in% motifs, ]
  snp.tbl <- snp.tbl[snpid %in% snpids, ]
  snp.tbl[, len_seq := nchar(ref_seq)]

  ## get the IUPAC subsequence for the motifs
  motif.tbl <- data.table(
    motif = motifs,
    IUPAC = sapply(motif.lib[motifs],
      function(x) GetIUPACSequence(x, prob = 0.25))
  )
  setkey(motif.tbl, motif)
  setkey(snp.tbl, snpid, snpbase)
  
  ncores <- min(c(ncores, length(snpids)))

#  startParallel(ncores)
    k <- as.integer(length(snpids) / ncores)

      if(ncores > 1) {
    if(Sys.info()[["sysname"]] == "Windows"){
      snow <- SnowParam(workers = ncores, type = "SOCK")
      motif_score_par_list<-bpmapply(function(x) match_subseq_par(i=x, par.k=k, par.ncores=ncores, par.snp.tbl=snp.tbl, par.snpids=snpids, par.motif.scores=motif.scores, par.motif=motif, par.motif.tbl=motif.tbl), seq(ncores), BPPARAM = snow,SIMPLIFY = FALSE)
    }else{
      motif_score_par_list<-bpmapply(function(x) match_subseq_par(i=x, par.k=k, par.ncores=ncores, par.snp.tbl=snp.tbl, par.snpids=snpids, par.motif.scores=motif.scores, par.motif=motif, par.motif.tbl=motif.tbl), seq(ncores), BPPARAM = MulticoreParam(workers = ncores),
                                SIMPLIFY = FALSE)
    }
  } else 
  {
    motif_score_par_list<-list(match_subseq_par(1, par.k=k, par.ncores=ncores, par.snp.tbl=snp.tbl, 
                                                par.snpids=snpids, par.motif.scores=motif.scores,
                                                par.motif=motif, par.motif.tbl=motif.tbl))
  }
  
#  endParallel()

  motif_score_tbl <- motif_score_par_list[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      motif_score_tbl <- rbind(motif_score_tbl,
                               motif_score_par_list[[i]])
    }
  }
 motif_score_tbl<-as.data.frame(motif_score_tbl, stringAsFacotrs=FALSE)
  return(motif_score_tbl)
}

#' @name ComputePValues
#' @title Compute p-values for affinity scores.
#' @description This function computes the p-values for allele-specific affinity
#'  scores and between-allele affinity score changes using the importance 
#'  sampling technique.
#' @param motif.lib A list object with the output format of function 
#' \code{\link{LoadMotifLibrary}}.
#' @param snp.info A list object with the output format of function 
#' \code{\link{LoadSNPData}}.
#' @param motif.scores A data.frame object containing at least the following 
#' columns:
#' \tabular{ll}{
#' motif \tab The name of the motif.\cr
#' log_lik_ref \tab The log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab The log-likelihood score for the SNP allele.\cr
#' }
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @param  testing.mc Monte Carlo sample size of 200 is considered. Do not
#'  change the default unless conducting a quick test. Default: FALSE
#' @param figdir A string for the path to print p-value plots for monitoring 
#' results. Default: NULL (no figure).
#' @return A data.frame extending \code{motif.scores} by the following 
#' additional columns:
#' \tabular{ll}{
#' pval_ref \tab P-values for scores on the reference allele.\cr
#' pval_snp \tab P-values for scores on the SNP allele.\cr
#' pval_cond_ref \tab Conditional p-values for scores on the reference 
#' allele.\cr
#' pval_cond_snp \tab Conditional p-values for scores on the SNP allele.\cr
#' pval_diff \tab P-values for the difference in scores between the reference 
#' and the SNP alleles.\cr
#' pval_rank \tab P-values for the log rank ratio between the reference and the 
#' SNP alleles.\cr
#' }
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo 
#' \email{chandler.c.zuo@@gmail.com}
#' @examples
#' data(example)
#' ComputePValues(motif_library, snpInfo, motif_scores$motif.scores, ncores = 2, testing.mc=TRUE)
#' @import Rcpp
#' @import data.table
#' @importFrom BiocParallel bpmapply MulticoreParam SnowParam
#' @useDynLib atSNP
#' @export
ComputePValues <- function(motif.lib, snp.info, motif.scores, ncores = 1, testing.mc=FALSE, figdir = NULL) {
  ncores <- min(c(ncores, length(motif.lib)))
  motif = NULL
  snpid = NULL
  snpbase = NULL
  pval_ref = NULL
  pval_snp = NULL
  pval_cond_ref = NULL
  pval_cond_snp = NULL
  pval_diff = NULL
  pval_rank = NULL
  
#  startParallel(ncores)
  
  results <- as.list(seq_along(motif.lib))
  nsets <- as.integer(length(motif.lib) / ncores)
  motif.scores <- as.data.table(motif.scores)
  prior <- snp.info$prior
  transition <- snp.info$transition
  
  if(Sys.info()[["sysname"]] == "Windows"){
    snow <- SnowParam(workers = ncores, type = "SOCK")
    results<-bpmapply(function(x) results_motif_par(i=x, par.prior=prior, par.transition=transition, par.motif.lib=motif.lib, par.motif.scores=motif.scores, par.testing.mc=testing.mc, par.figdir=figdir), seq_along(motif.lib), BPPARAM = snow,SIMPLIFY = FALSE)
  }else{
    results<-bpmapply(function(x) results_motif_par(i=x, par.prior=prior, par.transition=transition, par.motif.lib=motif.lib, par.motif.scores=motif.scores, par.testing.mc=testing.mc, par.figdir=figdir), seq_along(motif.lib), BPPARAM = MulticoreParam(workers = ncores),
                      SIMPLIFY = FALSE)
  }
#  endParallel()

  motif.scores_dt<-as.data.table(motif.scores)
  setkey(motif.scores_dt, motif, snpid, snpbase)
  
  for(i in seq_along(results)) {
    motif.scores_dt[results[[i]]$rowids, pval_ref := results[[i]]$pval_a[, 1]]
    motif.scores_dt[results[[i]]$rowids, pval_snp := results[[i]]$pval_a[, 2]]
    motif.scores_dt[results[[i]]$rowids, pval_cond_ref := results[[i]]$pval_cond[, 1]]
    motif.scores_dt[results[[i]]$rowids, pval_cond_snp := results[[i]]$pval_cond[, 2]]
    motif.scores_dt[results[[i]]$rowids, pval_diff := results[[i]]$pval_diff[, 1]]
    motif.scores_dt[results[[i]]$rowids, pval_rank := results[[i]]$pval_rank[, 1]]
  }
  motif.scores.pval<-as.data.frame(motif.scores_dt, stringsAsFactors=FALSE)
  return(motif.scores.pval)
}

#' @name GetIUPACSequence
#' @title Get the IUPAC sequence of a motif.
#' @description Convert the posotion weight matrix of a motif to the IUPAC 
#' sequence.
#' @param pwm The position weight matrix, with the columns representing 
#' A, C, G, T.
#' @param prob The probability threshold. Default: 0.25.
#' @return A character string.
#' @author Sunyoung Shin \email{sunyoung.shin@@utdallas.edu}, Chandler Zuo 
#' \email{chandler.c.zuo@@gmail.com}
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
