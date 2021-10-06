#' Title
#'
#' @param x
#' @param y
#' @param delta
#' @param intercept
#'
#' @return
#' @export
#'
#' @examples
coef_rankAFT <- function(x, y, delta, intercept = T) {
  ynew <- 1000 * (length(y))^2
  data1 <- data.frame(y, x)
  options(contrasts = c("contr.treatment", "contr.poly"))
  tempfit <- lm(y ~ ., x = TRUE, y = TRUE, data = data1)
  x <- as.matrix(tempfit$x[, -1])
  xn <- dimnames(x)[[2]]
  y <- tempfit$y
  dimnum <- dim(x)
  n1 <- dimnum[1]
  n2 <- dimnum[2]
  yy0 <- rep(y, rep(n1, n1))
  delta1 <- rep(delta, rep(n1, n1))
  yy1 <- rep(y, n1)
  yy <- delta1 * (yy0 - yy1)
  xx0 <- matrix(rep(as.vector(x), rep(n1, n1 * n2)), nrow = n1 * n1)
  xx1 <- t(matrix(rep(as.vector(t(x)), n1), nrow = n2))
  xx <- xx0 - xx1
  xxdif <- xx * delta1
  xnew <- apply(xxdif, 2, sum)
  xnew <- rbind(xxdif, -xnew)
  yynew <- c(yy, ynew)
  dimnames(xnew) <- list(NULL, xn)
  beta <- quantreg::rq(yynew ~ xnew - 1, method = "fn")$coef
  #beta <- fit$coef

  if (intercept) {
    res <- y-as.matrix(x)%*%beta
    # K-M estimate of the survival function of the error term (at the estimated beta)
    KM.fit <- survf(res,delta)
    # KM.fit$time is ordered
    time <- unique(sort(res))
    #surv <- KM.fit$sur # estimated survival function for the residual error
    jump <- KM.fit$prob #-surv # drop the last one if it is not a failure (in that case, surv[q-1]=surv[q])
    #jump <- c(1,surv[-q])-c(surv[-q],0) #force last one to be failure if it is not a failure
    alpha <- time%*%jump # intercept estimator of the AFT model

    #  predmatrix <- x - t(matrix(rep(apply(x, 2, mean), n1), ncol = n1))
    #  residuals <- y - predmatrix %*% as.matrix(betag)

    return(c(alpha, beta))
  } else {
    return(beta)
  }
}
