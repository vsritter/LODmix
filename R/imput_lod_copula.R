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
imput_lod_copula <- function(time, delta, covs) {
  if (!all.equal(dim(time), dim(delta))) stop("'time' and 'delta' must be of the same dimension")

  p <- ncol(time)
  q <- ncol(covs)
  beta <- array(dim = c(q, p))
  res <- array(dim = dim(time))
  km <- vector('list', length = p)
  for (i in 1:p) {
    # Rank AFT coefficients
    beta[,i] <- coef_rankAFT(y = time[,i], x = covs, delta = delta[,i],
                             intercept = F)

    # Residuals
    res[,i] <- c(time[,i] - c(beta[,i] %*% t(covs)))

    # Residual survival
    aux.km <- survf(res[,i], delta[,i])
    aux.km <- cbind(aux.km$time[aux.km$prob > 0],
                    aux.km$S[aux.km$prob > 0])

    km[[i]] <- rbind(c(min(res[,i]) - 0.0001, 1), aux.km)
    # km[[i]] <- aux.km

    # If KM does not reach 0, impute with Exp(lambda)
    if(min(aux.km[,2]) > 0) {
      lambda <- -log(min(aux.km[,2]))/max(aux.km[,1])
      aux.s <- seq(min(aux.km[,2]) + 0.0001, 0.0001, length.out = 50)
      aux.t <- qexp(aux.s, rate = lambda, lower.tail = F)

      km[[i]] <- rbind(km[[i]], cbind(aux.t, aux.s))
    }
  }

  # everyone that has at least one censored outcome
  idx <- which(rowSums(1 - delta) > 0)
  new.time <- time

  u <- rCopula(nrow(res), empCopula(pobs(res)))

  for (i in idx) {
    j0 <- which(delta[i,] == 0)
    j1 <- which(delta[i,] == 1)

    new.res <- sapply(j0, function(x) {
      s0 = approx(km[[x]][,1], km[[x]][,2], xout = res[i,x])$y

      if (is.na(s0)) {
        print('NA')
        return(NA)
      }

      cond.s = km[[x]][,2]/s0
      approx(cond.s, km[[x]][,1], xout = u[i,x])$y
    })

    new.time[i,j0] <- covs[i,] %*% beta[,j0] + new.res
  }

  return(new.time)
}
