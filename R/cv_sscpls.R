cv_sscpls_par <- function (x, y, fold = 10, K, lambda, abstol = 1e-03, reltol= 1e-03, max_itter = 500, compositional = FALSE, foldi =NULL, scale = F) {

  #-- Initalise --#
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  if(is.null(foldi)){
    foldi <- split(sample(1:n), rep(1:fold, length = n))
  }
  sdmat <- mspemat <- matrix(0, length(lambda), length(K))

  #-- loop the lambda --#
  for (i in 1:length(lambda)) {

    #-- parallel compute of cross validation --#
    mspemati <- foreach(j = 1:fold, .packages = c("sscpls"), .combine = rbind) %dopar% {
      omit <- foldi[[j]]
      mspes <- rep(0, length(K))
      object <- sscpls(x[-omit, , drop = FALSE], y[-omit, , drop = FALSE], ncomp = max(K),compositional = compositional,
                       lambda = lambda[i], scale = scale, abstol = abstol, reltol= reltol, max_itter = max_itter)

      newx <- x[omit, , drop = FALSE]
      #newx <- scale(newx, object$xmeans, scale = object$x_scale)
      betamat <- object$Beta

      for (k in K) {
        pred <- pred <- newx %*% betamat[[k]]# + matrix(1, nrow = nrow(newx), 1) %*% colMeans(y[-omit,,drop=F])
        mspes[ k - min(K) + 1 ] <- mean(apply((y[omit,] - pred)^2, 2, mean))
      }
      mspes
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

cv_sscpls <- function (x, y, fold = 10, K, lambda, abstol = 1e-03, reltol= 1e-03, max_itter = 100, compositional = FALSE, foldi=NULL, scale = F) {

  #-- Initalise --#
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  if(is.null(foldi)){
    foldi <- split(sample(1:n), rep(1:fold, length = n))
  }
  sdmat <- mspemat <- matrix(0, length(lambda), length(K))

  #-- loop the lambda --#
  for (i in 1:length(lambda)) {
    mspemati <- matrix(0, fold, length(K))
    for (j in 1:fold) {
      omit <- foldi[[j]]
      object <- sscpls(x[-omit, , drop = FALSE], y[-omit, , drop = FALSE], ncomp = max(K), compositional = compositional,
                       lambda = lambda[i], scale = F, abstol = abstol, reltol= reltol, max_itter = max_itter)
      newx <- x[omit, , drop = FALSE]
      newx <- scale(newx, object$xmeans, scale = object$x_scale)
      betamat <- object$Beta

      for (k in K) {
        pred <- newx %*% betamat[[k]] + matrix(1, nrow = nrow(newx), 1) %*% object$mu
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
             lambda.opt = lambda.opt, lambda.adj = lambda.adj, lambda.1se = lambda.1se,K = K.opt, foldi=foldi)
  invisible(cv)
}
