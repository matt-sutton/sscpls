cv_naive_spls <- function (x, y, fold = 10, K, lambda, deflate = "simpls") {
  
  #-- Initalise --#
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  foldi <- split(sample(1:n), rep(1:fold, length = n))
  sdmat <- mspemat <- matrix(0, length(lambda), length(K))
  
  #-- loop the lambda --#
  for (i in 1:length(lambda)) {
    mspemati <- matrix(0, fold, length(K))
    for (j in 1:fold) {
      omit <- foldi[[j]]
      object <- naive_spls(x[-omit, , drop = FALSE], y[-omit, , drop = FALSE], ncomp = max(K),
                       lambda = lambda[i], scale = F,deflate = deflate)
      mux <- colMeans(x[-omit, , drop = FALSE])
      newx <- x[omit, , drop = FALSE]
      newx <- scale(newx, mux, scale = F)
      betamat <- object$Beta
      
      for (k in K) {
        pred <- pred <- newx %*% betamat[[k]] + matrix(1, nrow = nrow(newx), 1) %*% colMeans(y[-omit,,drop=F])
        mspemati[j, (k - min(K) + 1)] <- mean(apply((y[omit,] - pred)^2, 2, mean))
      }
    }
    mspemat[i, ] <- apply(mspemati, 2, mean)
    sdmat[i, ] <- sqrt(apply(mspemati, 2, var))/sqrt(fold)
  }
  
  minpmse <- min(mspemat)
  se <- max(sdmat[which.min(mspemat)])
  rownames(mspemat) <- lambda
  colnames(mspemat) <- K
  mspecol <- apply(mspemat, 2, min)
  msperow <- apply(mspemat, 1, min)
  K.opt <- min(K[mspecol == minpmse])
  lambda.opt <- max(lambda[msperow == minpmse])
  lambda.adj <- max(lambda[msperow <= 1.1*minpmse])
  lambda.1se <- max(lambda[msperow <= (1 + se)*minpmse])
  
  rownames(mspemat) <- paste("lambda=", lambda)
  colnames(mspemat) <- paste("K =", K)
  cv <- list(mspemat = mspemat, sdmat = sdmat,
             lambda.opt = lambda.opt, lambda.adj = lambda.adj, lambda.1se = lambda.1se,K = K.opt)
  invisible(cv)
}
