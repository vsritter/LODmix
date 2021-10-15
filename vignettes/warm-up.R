rm(list = ls()); gc()

library(tidyverse)
library(survival)
library(LODmix)


# SETUP -------------------------------------------------------------------

set.seed(1234)

p = 5
q = 2

sig <- Matrix::bdiag(
  matrix(.1, nrow = 2, ncol = 2),
  matrix(.75, nrow = 3, ncol = 3))
Matrix::diag(sig) <- 1

sig <- diag(p)
sig <- 1/4^2*sig

alph <- -matrix(runif(p*(q+1), .1, .5), ncol = p, nrow = q+1) %>% round(2)

prop = sample(1 - seq(.5, .9, length.out = p))


# SIMULATE DATA -----------------------------------------------------------


pop_LODs <- gen_env_mix(n = 10^6, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        LOD = NULL, LOD.p = prop)

lod_data <- gen_env_mix(n = 200, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        return.full.data = TRUE,
                        LOD = pop_LODs)


# time = lod_data$lod.data
# delta = lod_data$delta
# covs = lod_data$X


mi_time <- imput_lod_tobit_v2(lod_data$lod.data, lod_data$delta, lod_data$X,
                              lod_data$full.data, sig, alph, r = 500)

coord <- which(lod_data$delta == 0, arr.ind = T
               )
lapply(1:r, function(x) {
  mi <- mi_time[,1]
  miT <- mi_time[mi == x,-c(1, 2)]


  theta[3,] <- lm(Y ~ cbind(X, Z))$coef
})

lm(Y ~ cbind(X, Z))$coef



# GAUSSIAN IMPUTATION -----------------------------------------------------

# imp1_time <- imput_lod_tobit_v1(lod_data$lod.data, lod_data$delta, lod_data$X,
#                                 lod_data$full.data, sig, alph)
# imp2_time <- imput_lod_copula(lod_data$lod.data, lod_data$delta, lod_data$X)

data.frame(
  mix = which(lod_data$delta >= 0, arr.ind = T)[,2],
  tres = imp1_time[[1]][lod_data$delta >= 0],
  nres = imp1_time[[2]][lod_data$delta >= 0]) %>%
  ggplot() +
  geom_density(aes(x = tres)) +
  geom_density(aes(x = nres), color = 2) +
  facet_wrap(vars(mix), ncol = 2, scales = 'free')

data.frame(
  mix = which(lod_data$delta == 0, arr.ind = T)[,2],
  tres = imp1_time[[1]][lod_data$delta == 0],
  nres = imp1_time[[2]][lod_data$delta == 0]) %>%
  ggplot() +
  geom_density(aes(x = tres)) +
  geom_density(aes(x = nres), color = 2) +
  facet_wrap(vars(mix), ncol = 2, scales = 'free')

data.frame(
  mix = which(lod_data$delta == 0, arr.ind = T)[,2],
  tres = lod_data$full.data[lod_data$delta == 0],
  nres = imp1_time[[3]][lod_data$delta == 0]) %>%
  ggplot() +
  geom_density(aes(x = tres)) +
  geom_density(aes(x = nres), color = 2) +
  facet_wrap(vars(mix), ncol = 2, scales = 'free')

data.frame(
  mix = which(lod_data$delta == 0, arr.ind = T)[,2],
  tres = lod_data$full.data[lod_data$delta == 0],
  nres = imp1_time[[3]][lod_data$delta == 0]) %>%
  ggplot() +
  geom_point(aes(x = nres, y = tres)) +
  facet_wrap(vars(mix), ncol = 2, scales = 'free')


mix <- 5

s0 <- survfit(Surv(lod_data$full.data[,mix], rep(1, nrow(lod_data$full.data))) ~ 1)
s1 <- survfit(Surv(lod_data$lod.data[,mix], lod_data$delta[,mix]) ~ 1)
s2 <- survfit(Surv(imp1_time[[3]][,mix], rep(1, nrow(imp1_time[[3]]))) ~ 1)
s3 <- survfit(Surv(imp2_time[,mix], rep(1, nrow(imp2_time))) ~ 1)

colMeans(1 - lod_data$delta)

plot(s0, conf.int = F)
lines(s1, conf.int = F, col = 2)
lines(s2, conf.int = F, col = 3)
abline(v = pop_LODs[mix])
lines(s3, conf.int = F, col = 4)



theta <- matrix(NA, ncol = 1+p+q, nrow = 4)

X <- lod_data$X
Z <- -lod_data$full.data
Y <- cbind(1, X, Z) %>% rowSums() + rnorm(nrow(X))

# full beta
theta[1,] <- lm(Y ~ cbind(X, Z))$coef

# LOD beta
Z <- -lod_data$lod.data
theta[2,] <- lm(Y ~ cbind(X, Z))$coef

# mvn imputation
Z <- -imp1_time[[3]]
theta[3,] <- lm(Y ~ cbind(X, Z))$coef

# copula imputation
Z <- -imp2_time
theta[4,] <- lm(Y ~ cbind(X, Z))$coef

colnames(theta) <- c('Intercept', paste0('X', 1:q), paste0('Z', 1:p))
round(theta-1, 4)


bind_cols(
  data.frame(lod_data$full.data) %>% rename_with(~gsub('X', 'Z1', .)),
  data.frame(imp1_time[[3]]) %>% rename_with(~gsub('X', 'Z2', .)),
  data.frame(imp2_time) %>% rename_with(~gsub('X', 'Z3', .))) %>%
  pivot_longer(everything()) %>%
  mutate(mix = substr(name, 3, 3),
         name = substr(name, 1, 2)) %>%
  ggplot() +
  geom_density(aes(x = value, color = name), size = 1) +
  geom_vline(aes(xintercept = pop_LODs),
             data = data.frame(pop_LODs, mix = 1:5)) +
  facet_wrap(vars(mix), ncol = 2, scales = 'free') +
  theme_bw()



# -------------------------------------------------------------------------

rep_sim_lod <- function(n, p, q, sig) {
  lod_data <- gen_env_mix(n = n, p = p, q = q, error.dist = 'mvnorm',
                          alpha = alph, sigma = sig,
                          return.full.data = TRUE,
                          LOD = pop_LODs)

  theta <- matrix(NA, ncol = 1+p+q, nrow = 4)

  X <- lod_data$X
  Z <- -lod_data$full.data
  Y <- cbind(1, X, Z) %>% rowSums() + rnorm(nrow(X))

  # full beta
  theta[1,] <- lm(Y ~ cbind(X, Z))$coef

  # LOD beta
  Z <- -lod_data$lod.data
  theta[2,] <- lm(Y ~ cbind(X, Z))$coef

  # mvn imputation
  mvn_time <- imput_lod_tobit(lod_data$lod.data, lod_data$delta, lod_data$X)
  Z <- -mvn_time
  theta[3,] <- lm(Y ~ cbind(X, Z))$coef

  # copula imputation
  cpl_time <- imput_lod_copula(lod_data$lod.data, lod_data$delta, lod_data$X)
  Z <- -cpl_time
  theta[4,] <- lm(Y ~ cbind(X, Z))$coef

  colnames(theta) <- paste0('theta', 1:ncol(theta))
  return(data.frame('method' = 1:4, theta))
}


res1 <- lapply(1:20, function(x) {
  print(x)
  rep_sim_lod(n = 1000, p = 5, q = 2, sig = sig)
}) %>% bind_rows()


res2 <- lapply(1:50, function(x) {
  print(x)
  rep_sim_lod(n = 1000, p = 5, q = 2, sig = diag(5))
}) %>% bind_rows()


colMeans(1 - lod_data$delta)

res1 %>%
  group_by(method) %>%
  summarise_all(list(b = ~mean(.)-1)) %>%
  round(4)

res2 %>%
  group_by(method) %>%
  summarise_all(list(b = ~mean(.)-1)) %>%
  round(4)


# saveRDS(res1, 'res1.rds')
# saveRDS(res2, 'res2.rds')



# -------------------------------------------------------------------------



