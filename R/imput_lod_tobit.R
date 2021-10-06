#' Title
#'
#' @param time
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
imput_lod_tobit <- function(time, delta, covs) {
  if (!all.equal(dim(time), dim(delta))) stop("'time' and 'delta' must be of the same dimension")

  p <- ncol(time)
  res <- array(dim = dim(time))
  scale <- rep(NA, p)
  for (i in 1:p) {
    # Fit linear model for censored data
    fit <- survival::survreg(survival::Surv(time[,i], delta[,i]) ~ covs, dist = 'gaussian')
    res[,i] <- resid(fit)
    scale[i] <- fit$scale
  }

  # scaled residuals covariance matrix
  S <- var(res %*% diag(scale/apply(res, 2, sd)))

  # Get censored residuals
  idx0 <- which(delta[,i] == 0)
  res0 <- resid(fit)[idx0]

  return(S)
}
