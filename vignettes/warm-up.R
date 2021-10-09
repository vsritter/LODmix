# rm(list = ls()); gc()




# SETUP -------------------------------------------------------------------

set.seed(1234)

p = 5
q = 2

sig <- Matrix::bdiag(
  matrix(.25, nrow = 2, ncol = 2),
  matrix(.75, nrow = 3, ncol = 3))
Matrix::diag(sig) <- 1

sig <- 1/2*sig


alph <- -matrix(runif(p*(q+1), .1, .5), ncol = p, nrow = q+1) %>% round(2)

prop = sample(1 - seq(.5, .9, length.out = p))


# SIMULATE DATA -----------------------------------------------------------


pop_LODs <- gen_env_mix(n = 10^5, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        LOD = NULL, LOD.p = prop)

lod_data <- gen_env_mix(n = 400, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        return.full.data = TRUE,
                        LOD = pop_LODs)

colMeans(1 - lod_data$delta)

# GAUSSIAN IMPUTATION -----------------------------------------------------

imp_time <- imput_lod_tobit(lod_data$lod.data, lod_data$delta, lod_data$X)


library(survival)
s0 <- survfit(Surv(lod_data$full.data[,2], rep(1, nrow(lod_data$full.data))) ~ 1)
s1 <- survfit(Surv(lod_data$lod.data[,2], lod_data$delta[,2]) ~ 1)
s2 <- survfit(Surv(imp_time[,2], rep(1, nrow(imp_time))) ~ 1)

colMeans(1 - lod_data$delta)

plot(s0, conf.int = F)
lines(s1, conf.int = F, col = 2)
lines(s2, conf.int = F, col = 3)







time = lod_data$lod.data
delta = lod_data$delta
covs = lod_data$X

# -------------------------------------------------------------------------

p <- ncol(time)
q <- ncol(covs)
res <- array(dim = dim(time))
scale <- rep(NA, p)
beta <- array(dim = c(q+1, p))
for (i in 1:p) {
  # Fit linear model for censored data
  fit <- survival::survreg(survival::Surv(time[,i], delta[,i]) ~ covs,
                           dist = 'gaussian')
  res[,i] <- resid(fit)
  scale[i] <- fit$scale
  beta[,i] <- coef(fit)
}

# scaled residuals covariance matrix
S <- var(res %*% diag(scale/apply(res, 2, sd)))

# everyone that has at least one censored outcome
idx <- which(rowSums(1 - delta) > 0)
new.time <- time

for (i in idx) {
  j0 <- which(delta[i,] == 0)
  j1 <- which(delta[i,] == 1)

  if (length(j1) > 0) {
    s00 <- as.matrix(S[j0,j0])
    s11 <- as.matrix(S[j1,j1])
    s01 <- as.matrix(S[j0,j1])

    if (length(j0) == 1)
      s01 <- t(s01)

    cond.S <- s00 - s01 %*% solve(s11) %*% t(s01)
    new.res <- tmvtnorm::rtmvnorm(1, lower = res[i,j0], sigma = cond.S,
                                  algorithm = 'gibbs')

    new.time[i,j0] <- c(1, covs[i,]) %*% beta[,j0] + new.res
  } else {
    new.res <- tmvtnorm::rtmvnorm(1, lower = res[i,], sigma = S,
                                  algorithm = 'gibbs')

    new.time[i,] <- c(1, covs[i,]) %*% beta[,j0] + new.res
  }
}

new.time %>% tail()
time %>% tail()



mvtnorm::pmvnorm(upper = c(1.96, 1.96), sigma = diag(2), keepAttr = F)













