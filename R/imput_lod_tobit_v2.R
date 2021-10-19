#' Title
#'
#' @param time
#' @param delta
#' @param covs
#'
#' @return
#' @export
#'
#' @examples
imput_lod_tobit_v2 <- function(time, delta, covs, r) {#, tfull, s, a, r) {
  if (!all.equal(dim(time), dim(delta))) stop("'time' and 'delta' must be of the same dimension")

  p <- ncol(time)
  q <- ncol(covs)
  res <- array(dim = dim(time))
  scale <- rep(NA, p)
  beta <- array(dim = c(q+1, p))
  # tres <- array(dim = dim(time))
  for (i in 1:p) {
    # Fit linear model for censored data
    fit <- survival::survreg(survival::Surv(time[,i], delta[,i]) ~ covs,
                             dist = 'gaussian')
    res[,i] <- resid(fit)
    scale[i] <- fit$scale
    beta[,i] <- coef(fit)
    # beta[,i] <- a[,i]
    # res[,i] <- c(time[,i] - c(beta[,i] %*% t(cbind(1, covs))))
    # tres[,i] <- c(tfull[,i] - c(beta[,i] %*% t(cbind(1, covs))))
  }

  # scaled residuals covariance matrix
  S <- var(res %*% diag(scale/apply(res, 2, sd)))
  # S <- s

  # everyone that has at least one censored outcome
  idx <- which(rowSums(1 - delta) > 0)
  l = length(idx)
  mi.res <- array(dim = c(l*r, p+2))

  for (i in 1:l) {
    j0 <- which(delta[idx[i],] == 0)
    xbi <- c(1, covs[idx[i],]) %*% beta
    new.res <- tmvtnorm::rtmvnorm(r, lower = res[idx[i],],
                                  sigma = diag(diag(as.matrix(S))),
                                  algorithm = 'gibbs')
    # new.res <- mvtnorm::rmvnorm(r, sigma = diag(diag(as.matrix(S))))
    # new.res <- t(t(new.res) + res[idx[i],])

    mi.res[1:r+r*(i-1),1] <- 1:r
    mi.res[1:r+r*(i-1),2] <- rep(idx[i], r)
    mi.res[1:r+r*(i-1),j0+2] <- t(xbi[j0] + t(new.res[,j0]))
  }


  return(mi.res)
  # return(list(tres, nres, new.time))
  # return(new.time)
}
