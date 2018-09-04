.structure <- function(pval_mat) {
    if(is.matrix(pval_mat)) {
      if(nrow(pval_mat) > 1) {
        id <- apply(pval_mat[, c(2, 4)], 1, which.min)
        return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
                     pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
      }
    }
    
    pval_mat <- as.vector(pval_mat)
    id <- which.min(pval_mat[c(2, 4)])
    return(cbind(pval_mat[c(1, 3)][id],
                 pval_mat[c(2, 4)][id]))
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}

myStrSplit <- function(x, split) {
  ret <- list(seq_along(x))
  for(i in seq_along(x)) {
    ret[[i]] <- x[i]
    for(sp in split) {
      ret[[i]] <- unlist(strsplit(ret[[i]], sp))
      ret[[i]] <- ret[[i]][nchar(ret[[i]]) > 0]
      if(length(ret[[i]]) == 0)
        break
    }
  }
  return(ret)
}

#' @import doParallel
startParallel <- function(ncores) {
  if(.Platform$OS.type == "unix") {
    registerDoParallel(ncores)
  } else {
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    return("cl")
  }
}

#' @import doParallel
endParallel <- function() {
  if(.Platform$OS.type != "unix") {
    stopCluster(cl)
  }
}
