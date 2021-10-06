#' Title
#'
#' @param time
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
survf <- function(time, delta){
  ndim <- length(time)
  uniquetime1 <- unique(sort(time))
  npts1 <- length(uniquetime1)

  tmp1 <- matrix(rep(time,each = npts1), npts1, ndim)
  tmp3 <- matrix(rep(delta,each = npts1), npts1, ndim)

  lambda1 <- apply(ifelse(tmp1 == uniquetime1 & tmp3 == 1,1,0),1,sum)/apply(ifelse(tmp1 >= uniquetime1,1,0),1,sum)
  S <- cumprod(1 - lambda1)
  jump <- -diff(c(1,S))
  return(list(time = uniquetime1, S = S, prob = jump))
}
