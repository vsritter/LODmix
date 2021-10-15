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
imput_lod_tobit_v2 <- function(time, delta, covs, tfull, s, a, r) {
  if (!all.equal(dim(time), dim(delta))) stop("'time' and 'delta' must be of the same dimension")

  p <- ncol(time)
  q <- ncol(covs)
  res <- array(dim = dim(time))
  scale <- rep(NA, p)
  beta <- array(dim = c(q+1, p))
  tres <- array(dim = dim(time))
  for (i in 1:p) {
    # Fit linear model for censored data
    fit <- survival::survreg(survival::Surv(time[,i], delta[,i]) ~ covs,
                             dist = 'gaussian')
    # res[,i] <- resid(fit)
    # scale[i] <- fit$scale
    # beta[,i] <- coef(fit)
    beta[,i] <- a[,i]
    res[,i] <- c(time[,i] - c(beta[,i] %*% t(cbind(1, covs))))
    tres[,i] <- c(tfull[,i] - c(beta[,i] %*% t(cbind(1, covs))))
  }

  # scaled residuals covariance matrix
  # S <- var(res %*% diag(scale/apply(res, 2, sd)))
  S <- s

  # everyone that has at least one censored outcome
  idx <- which(rowSums(1 - delta) > 0)
  l = length(idx)
  mi.res <- array(dim = c(l*r, p+2))

  for (i in 1:l) {
    new.res <- tmvtnorm::rtmvnorm(r, lower = res[idx[i],],
                                  sigma = diag(diag(as.matrix(S))),
                                  algorithm = 'gibbs')

    mi.res[1:r+l*(i-1),1] <- rep(i, r)
    mi.res[1:r+l*(i-1),2] <- idx[i]
    mi.res[1:r+l*(i-1),3:(p+2)] <- new.res
  }

  return(mi.res)
  # return(list(tres, nres, new.time))
  # return(new.time)
}
