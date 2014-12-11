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
