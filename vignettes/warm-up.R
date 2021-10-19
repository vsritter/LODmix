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

# sig <- diag(p)
sig <- 1/4^2*sig

alph <- -matrix(runif(p*(q+1), .1, .5), ncol = p, nrow = q+1) %>% round(2)

prop = sample(1 - seq(.5, .9, length.out = p))


# SIMULATE DATA -----------------------------------------------------------


pop_LODs <- gen_env_mix(n = 10^6, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        LOD = NULL, LOD.p = prop)

lod_data <- gen_env_mix(n = 400, p = p, q = q, error.dist = 'mvnorm',
                        alpha = alph, sigma = sig,
                        return.full.data = TRUE,
                        LOD = pop_LODs)


# time = lod_data$lod.data
# delta = lod_data$delta
# covs = lod_data$X

r = 500
mi_time <- imput_lod_tobit_v2(lod_data$lod.data, lod_data$delta, lod_data$X, r = r)


X <- lod_data$X
Z <- -lod_data$full.data
Y <- cbind(1, X, Z) %>% rowSums() + rnorm(nrow(X))

Time <- lod_data$lod.data
theta2 <- lapply(1:r, function(x) {
  Tmi <- mi_time[mi_time[,1] == x,-c(1,2)]
  Time[lod_data$delta == 0] <- Tmi[!is.na(Tmi)]

  Z <- -Time
  lm(Y ~ cbind(X, Z))$coef
}) %>% bind_rows()
colnames(theta2) <- c('int', paste0('demo', 1:q), paste0('mix', 1:p))

theta2 <- theta2 %>%
  summarise_all(mean)

round(theta2-1, 4)



# GAUSSIAN IMPUTATION -----------------------------------------------------

# imp1_time <- imput_lod_tobit_v1(lod_data$lod.data, lod_data$delta, lod_data$X,
#                                 lod_data$full.data, sig, alph)
imp2_time <- imput_lod_copula(lod_data$lod.data, lod_data$delta, lod_data$X)

# data.frame(
#   mix = which(lod_data$delta >= 0, arr.ind = T)[,2],
#   tres = imp1_time[[1]][lod_data$delta >= 0],
#   nres = imp1_time[[2]][lod_data$delta >= 0]) %>%
#   ggplot() +
#   geom_density(aes(x = tres)) +
#   geom_density(aes(x = nres), color = 2) +
#   facet_wrap(vars(mix), ncol = 2, scales = 'free')
#
# data.frame(
#   mix = which(lod_data$delta == 0, arr.ind = T)[,2],
#   tres = imp1_time[[1]][lod_data$delta == 0],
#   nres = imp1_time[[2]][lod_data$delta == 0]) %>%
#   ggplot() +
#   geom_density(aes(x = tres)) +
#   geom_density(aes(x = nres), color = 2) +
#   facet_wrap(vars(mix), ncol = 2, scales = 'free')
#
# data.frame(
#   mix = which(lod_data$delta == 0, arr.ind = T)[,2],
#   tres = lod_data$full.data[lod_data$delta == 0],
#   nres = imp1_time[[3]][lod_data$delta == 0]) %>%
#   ggplot() +
#   geom_density(aes(x = tres)) +
#   geom_density(aes(x = nres), color = 2) +
#   facet_wrap(vars(mix), ncol = 2, scales = 'free')
#
# data.frame(
#   mix = which(lod_data$delta == 0, arr.ind = T)[,2],
#   tres = lod_data$full.data[lod_data$delta == 0],
#   nres = imp1_time[[3]][lod_data$delta == 0]) %>%
#   ggplot() +
#   geom_point(aes(x = nres, y = tres)) +
#   facet_wrap(vars(mix), ncol = 2, scales = 'free')


# mix <- 5
#
# s0 <- survfit(Surv(lod_data$full.data[,mix], rep(1, nrow(lod_data$full.data))) ~ 1)
# s1 <- survfit(Surv(lod_data$lod.data[,mix], lod_data$delta[,mix]) ~ 1)
# s2 <- survfit(Surv(imp1_time[[3]][,mix], rep(1, nrow(imp1_time[[3]]))) ~ 1)
# s3 <- survfit(Surv(imp2_time[,mix], rep(1, nrow(imp2_time))) ~ 1)
#
# colMeans(1 - lod_data$delta)
#
# plot(s0, conf.int = F)
# lines(s1, conf.int = F, col = 2)
# lines(s2, conf.int = F, col = 3)
# abline(v = pop_LODs[mix])
# lines(s3, conf.int = F, col = 4)



theta <- matrix(NA, ncol = 1+p+q, nrow = 4)

X <- lod_data$X
Z <- -lod_data$full.data
# Y <- cbind(1, X, Z) %>% rowSums() + rnorm(nrow(X))

# full beta
theta[1,] <- lm(Y ~ cbind(X, Z))$coef

# LOD beta
Z <- -lod_data$lod.data
theta[2,] <- lm(Y ~ cbind(X, Z))$coef

# mvn imputation
Z <- -imp1_time[[3]]
theta[3,] <- lm(Y ~ cbind(X, Z))$coef

# copula imputation
# Z <- -imp2_time
# theta[4,] <- lm(Y ~ cbind(X, Z))$coef

colnames(theta) <- c('int', paste0('demo', 1:q), paste0('mix', 1:p))

round(theta-1, 4)
round(theta2-1, 4)


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

  theta <- matrix(NA, ncol = 1+p+q, nrow = 5)

  X <- lod_data$X
  Z <- -lod_data$full.data
  Y <- cbind(1, X, Z) %>% rowSums() + rnorm(nrow(X))

  # full beta
  theta[1,] <- lm(Y ~ cbind(X, Z))$coef

  # LOD beta
  Z <- -lod_data$lod.data
  Z[lod_data$delta == 0] <- Z[lod_data$delta == 0]/sqrt(2)
  theta[2,] <- lm(Y ~ cbind(X, Z))$coef

  # mvn imputation
  # mvn_time <- imput_lod_tobit_v1(lod_data$lod.data, lod_data$delta, lod_data$X,
  #                                lod_data$full.data, sig, alph)[[3]]
  mvn_time <- imput_lod_tobit_v1(lod_data$lod.data, lod_data$delta, lod_data$X)
  Z <- -mvn_time
  theta[3,] <- lm(Y ~ cbind(X, Z))$coef

  # copula imputation
  cpl_time <- imput_lod_copula(lod_data$lod.data, lod_data$delta, lod_data$X)
  Z <- -cpl_time
  theta[4,] <- lm(Y ~ cbind(X, Z))$coef

  # MI
  # mi_time <- imput_lod_tobit_v2(lod_data$lod.data, lod_data$delta, lod_data$X,
  #                               lod_data$full.data, sig, alph, r = 500)
  mi_time <- imput_lod_tobit_v2(lod_data$lod.data, lod_data$delta, lod_data$X, r = 500)

  Time <- lod_data$lod.data
  theta2 <- lapply(1:r, function(x) {
    Tmi <- mi_time[mi_time[,1] == x,-c(1,2)]
    Time[lod_data$delta == 0] <- Tmi[!is.na(Tmi)]

    Z <- -Time
    lm(Y ~ cbind(X, Z))$coef
  }) %>% simplify2array()

  theta[5,] <- rowMeans(theta2)

  # colnames(theta) <- c('int', paste0('demo', 1:q), paste0('mix', 1:p))
  return(data.frame('method' = c('Full data', 'LOD.sqrt', 'tmvn', 'Copula', 'MI-like'), theta))
}

library(parallel)
# res1 <- mclapply(1:200, function(x) {
#   rep_sim_lod(n = 400, p = 5, q = 2, sig = sig)
# }, mc.cores = 16) %>% bind_rows()
# colnames(res1)[-1] <- c('int', paste0('demo', 1:q), paste0('mix', 1:p))

# res2 <- mclapply(1:200, function(x) {
#   rep_sim_lod(n = 400, p = 5, q = 2, sig = 1/4^2*diag(5))
# }, mc.cores = 16) %>% bind_rows()

colnames(res1)[-1] <- c('int', paste0('demo', 1:q), paste0('mix', 1:p))
res1 %>%
  group_by(method) %>%
  summarise_all(list(b = ~round(mean(.)-1, 4))) %>%
  rbind(c('Prop. LOD', rep(NA, 1+q), colMeans(1 - lod_data$delta))) %>%
  mutate(i = c(6,2,3,5,4,1)) %>% arrange(i) %>% select(-i)

print.p <- round(colMeans(1 - lod_data$delta)*100, 0)
nm <- c('Intercept',
        paste0('X', 1:q),
        paste0('Z', 1:p, ' (', print.p, '% LOD)'))
colnames(res1)[-1] <- nm


gg <- res1 %>%
  mutate(
    method = gsub(' data|.sqrt', '', method),
    method = factor(method, levels = unique(method)[c(1:3,5,4)])) %>%
  filter(!(method %in% c('Copula', 'MI-like'))) %>%
  pivot_longer(-method) %>%
  mutate(name = factor(name, levels = c(nm[1:5],'',nm[6:8]))) %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = value)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(vars(name), ncol = 3, drop = F) +
  labs(x = NULL, y = 'Bias')

library(gridExtra)
g <- ggplotGrob(gg)
# get the grobs that must be removed
rm_grobs <- g$layout$name %in% c("panel-2-3", "strip-t-3-2")
g$grobs[rm_grobs] <- NULL
g$layout <- g$layout[!rm_grobs, ]
plot(g)


res1 %>%
  mutate(
    method = gsub(' data|.sqrt', '', method),
    method = factor(method, levels = unique(method)[c(1:3,5,4)])) %>%
  pivot_longer(-method) %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = value)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(vars(name), ncol = 4)

res1 %>%
  mutate(
    method = gsub(' data|.sqrt', '', method),
    method = factor(method, levels = unique(method)[c(1:3,5,4)])) %>%
  pivot_longer(-method) %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = value)) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(vars(name), ncol = 4, scales = 'free')




# saveRDS(res1, 'res1.rds')
# saveRDS(res2, 'res2.rds')



# -------------------------------------------------------------------------



