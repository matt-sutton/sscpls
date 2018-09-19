show_nonzero <- function(X){
  X <- as.matrix(X)
  rownames(X) <- 1:nrow(X)
  return(X[which(rowSums(abs(X)) > 0),] )
}


#-- Internal funciton for vector normalisation --#
normalise <- function(x, xnorm){
  xnorm <- drop(xnorm + 1*(xnorm == 0))
  return(x/(xnorm))
}

#-- Internal funciton to get projection given basis X --#
get_proj <- function(X){
  X <- as.matrix(X)
  if(all(X == 0)) return(tcrossprod(X))

  L <- t(suppressWarnings(chol(crossprod(X), pivot = TRUE)))
  r <- attr(L, "rank")
  piv <- attr(L, "pivot")
  Qt <- forwardsolve(L, t(X[, piv]), r)

  ## P = QQ'
  H <- crossprod(Qt)
  return(H)
}

# Computes x = A^+ * b using an SVD
# (slow but stable)

svdsolve <- function(A,b,rtol) {
  s = svd(A)
  di = s$d
  ii = di>rtol
  di[ii] = 1/di[ii]
  di[!ii] = 0
  return(list(x=s$v%*%(di*(t(s$u)%*%b)),q=sum(ii)))
}

