Cal_MSPE <- function(x, y, oldX, Best){
  # Take independent dataset, and estimated Beta
  
  mux <- colMeans(oldX)
  newx <- x
  newx <- scale(newx, mux, scale = F)
  
  pred <- newx %*% Best + matrix(1, nrow = nrow(newx), 1) %*% colMeans(y[,,drop=F])
  
  return(mean(sd(y - newx%*%Best)))
}

Cal_ME <- function(B, B_True, V){
  # Take independent dataset, and estimated Beta
  return( t(B - B_True)%*%V%*%(B - B_True) )
}

feature_selection <- function(B_est, B_True, tol = 10^-6){
  
  B <- as.factor((abs(c(0,B_est, B_True)) == 0)+0.0);
  B <- B[-1]
  p <- length(B_est)
  Bhat <- B[1:p]; Bref <- B[1:p + p]
  
  sens <- caret::sensitivity(Bhat, Bref)
  spec <- caret::specificity(Bhat, Bref)
  mccr <- mccr::mccr(Bhat, Bref) 
  
  return(c(sens, spec, mccr))
}

Cal_ME <- function(B, B_True, V){
  # Take independent dataset, and estimated Beta
  return( t(B - B_True)%*%V%*%(B - B_True) )
}
