#' Title
#'
#' @param n
#' @param p
#' @param q
#' @param error.dist
#' @param alpha
#' @param return.T
#' @param mu.X
#' @param keep.full.data
#' @param LOD
#' @param LOD.p
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
gen_env_mix <- function(n, p, q, error.dist = 'mvnorm', alpha, sigma,
                        return.T = TRUE, mu.X = NULL,
                        return.full.data = FALSE,
                        LOD = NULL, LOD.p = NULL) {

  if (is.null(mu.X))
    mu.X <- rep(0, q)

  if (error.dist == 'mvnorm')
    X <- mvtnorm::rmvnorm(n, mean = mu.X, sigma = diag(q))

  xi <- mvtnorm::rmvnorm(n, sigma = as.matrix(sigma))
  Time <- cbind(1, X) %*% alpha + xi

  if (return.full.data)
    Tfull <- Time

  if (is.null(LOD)) {
    if (is.null(LOD.p)) stop("when 'LOD = NULL', LOD.p must be provided")

    limit <- vector('double', length = p)
    for (i in 1:p)
      limit[i] <- quantile(Time[,i], probs = 1 - LOD.p[i])

    return(limit)
  }

  Delta <- array(0L, dim = dim(Time))
  for (i in 1:p) {
    Time[,i] <- ifelse(Time[,i] >= LOD[i], LOD[i], Time[,i])
    Delta[,i] <- ifelse(Time[,i] >= LOD[i], 0L, 1L)
  }

  if (return.full.data)
    return(list('X' = X,
                'full.data' = Tfull,
                'lod.data' = Time,
                'delta' = Delta))

  return(list('X' = X,
              'lod.data' = Time,
              'delta' = Delta))
}
