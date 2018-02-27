show_nonzero <- function(X){
  X <- as.matrix(X)
  return(X[which(rowSums(abs(X)) > 0),] )
}
