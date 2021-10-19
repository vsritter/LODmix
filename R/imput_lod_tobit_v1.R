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
imput_lod_tobit_v1 <- function(time, delta, covs) {#, tfull, s, a) {
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
  new.time <- time
  # nres <- tres

  for (i in idx) {
    j0 <- which(delta[i,] == 0)
    j1 <- which(delta[i,] == 1)

    # if (length(j1) > 0) {
    #   s00 <- as.matrix(S[j0,j0])
    #   s11 <- as.matrix(S[j1,j1])
    #   s01 <- as.matrix(S[j0,j1])
    #
    #   if (length(j0) == 1)
    #     s01 <- t(s01)
    #
    #   cond.S <- s00 - s01 %*% solve(s11) %*% t(s01)
    #   new.res <- tmvtnorm::rtmvnorm(1, lower = res[i,j0],
    #                                 sigma = diag(length(j0))*diag(as.matrix(cond.S)),
    #                                 algorithm = 'gibbs')
    #
    #   new.time[i,j0] <- c(1, covs[i,]) %*% beta[,j0] + new.res
    # } else {
    #   new.res <- tmvtnorm::rtmvnorm(1, lower = res[i,],
    #                                 sigma = diag(p)*diag(as.matrix(S)),
    #                                 algorithm = 'gibbs')
    #
    #   new.time[i,] <- c(1, covs[i,]) %*% beta + new.res
    # }

    # sd = diag(as.matrix(S))^.5
    # a = res[i,]/sd
    # fa = dnorm(a)
    # z = 1 - pnorm(a, lower.tail = T)
    # v = sd^2*(1+(a*fa-0)/z-((fa-0)/z)^2)

    # new.res <- fa/z*sd

    new.res <- tmvtnorm::rtmvnorm(1, lower = res[i,],
                                  sigma = diag(diag(as.matrix(S))),
                                  algorithm = 'gibbs')

    # nres[i,j0] <- new.res[j0]

    new.time[i,j0] <- c(1, covs[i,]) %*% beta[,j0] + new.res[j0]
  }

  # return(list(tres, nres, new.time))
  return(new.time)
}
