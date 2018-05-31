#-- Some simple internal funcitons --#
prox_soft <- function(x, lambda){
  return(sign(x)*pmax(0, abs(x) - lambda))
}

prox_l2 <- function(u){
  nrm <- drop(sqrt(crossprod(u)))
  if(nrm  < 1 ){
    return(u)
  }  else {
    return(u/nrm)
  }
}

#-- Projection operator --#

proj_compositional <- function(u, proj){
  if(all(u == 0)) {return(u)}
  simplex_root<- function(eta, v){
    res <- proj%*%(v - eta) - 1e-8
    return(sum(res))
  }
  eta <- try(uniroot(simplex_root, interval = c(-max(abs(u)), max(abs(u))), v=u)$root)
  return(proj%*%(u - eta))
}

proj_orthogonal <- function(u, proj){
    return(proj%*%u)
}
