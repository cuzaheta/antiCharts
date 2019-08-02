#---- 
# file: Profile.R
# author: CU
# date:
# description: Implementación de metodología con datos de perfiles no lineales
# input: Datos (suma de scores estandarizados)
# output:
#    - Resultados
#    - 
#---
rm(list = ls())

# packages ----------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(furrr)
library(GGally)
# library(tictoc)

# paths -------------------------------------------------------------------
inPath <- file.path('input')
inFunc <- file.path('src','func')
outPath <- file.path('output')

# functions ---------------------------------------------------------------
pathFunc <- file.path(inFunc, 'funciones.R')
source(pathFunc, encoding = 'UTF8')


# read data ---------------------------------------------------------------
inSumal <- file.path(inPath, 'Datos Simulados.txt')
datosSum <- read.table(inSumal, header = T) %>% 
  as.tibble() %>% 
  set_colnames(paste0('Score',1:ncol(.)))

inAzuc <- file.path(inPath, 'Datos Azucar.txt')
datosAzuc <- read.table(inAzuc, header = T) %>% 
  as.tibble() %>% 
  set_colnames(paste0('Score',1:ncol(.)))


# SumaScores --------------------------------------------------------------
nSum <- nrow(datosSum)

datosSum %>% 
  map(summary)

pairsSum <- ggpairs(datosSum) 

cov(datosSum) 

datosSum %>% 
  map_dbl(var)

# Series 
acf(datosSum) # autocorrelacion no tiene

acf(datosSum[,4])
Box.test(datosSum[,4], 1)

datosSum %>% 
  gather(score, x) %>% 
  ggplot(aes(x)) +
  geom_density(aes(fill = score), alpha = 0.5)

datosSum %>% 
  mutate(obs = seq_along(Score2)) %>% 
  gather(score, x, -obs) %>% 
  ggplot(aes(obs, x, colour=score)) +
  geom_line()

datosSum %>% 
  rowMeans() %>% 
  tibble(obs=seq_along(.), x = .) %>% 
  ggplot(aes(obs, x)) +
  geom_line()

# Extraer la media
vec_media <- colMeans(datosSum)
vec_media

tblScores <- datosSum - rep(1, nrow(datosSum))%*%t(vec_media)
colMeans(tblScores)

# map2_dfc(datosSum, vec_media, ~.x - .y)
# as.tibble(tblScores)


# qiu ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 5
p <- ncol(tblScores)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

g <- distr.g(tblScores, ar, Ss, T)
vec.yQiu <- carta.MACusum(tblScores, ar, Ss, k, F)

plot.ts(vec.yQiu)

# Simulación limites de control
l <- 9.51
tictoc::tic()
future_arl(l, k, g, 1e5, F, ar, Ss)
tictoc::toc()

ls <- seq(9.5,9.52, 0.01)
tbl.ARLqiu <- tibble(ls) %>% 
  mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 50e3, F, ar, Ss)))
tbl.ARLqiu

# Carta de control con limite
tibble(yh = vec.yQiu, obs = seq_along(vec.yQiu)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l) +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Qiu con datos de Suma de Scores',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# Cada cuanto se reinicia? 
set.epsi2 <- apply(tblScores, 1, epsilon_h, ar, Ss, F) %>% t

s1 <- s2 <- numeric(length(g))

vec.y <- numeric(nrow(tblScores))
q.y <- numeric(nrow(tblScores))
for (ii in 1:nrow(tblScores)) {
  e1t <- set.epsi2[ii,]
  qt <- Qt(s1,s2,e1t,g)
  q.y[ii] <- qt
  s12 <- Sh.12(s1,s2,e1t,qt,g,k)
  s1 <- s12$S1
  s2 <- s12$S2
  vec.y[ii] <- max(0,qt - k)
}

sum(q.y <= k)
verti_line <- seq_along(q.y)[q.y <= k]

tibble(yh = vec.yQiu, obs = seq_along(vec.yQiu)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = verti_line, colour = 'red') +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Qiu con datos de Suma de Scores',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# tictoc::tic()
# mean(pbMACusum.RL(8, k, g, 1e4, F, ar, Ss))
# tictoc::toc()

# hwang ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 5
p <- ncol(tblScores)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

g <- distr.g(tblScores, ar, Ss, T)
vec.yHwang <- carta.MACusum(tblScores, ar, Ss, k, T)

plot.ts(vec.yHwang)

# Simulación limites de control
Vmu <- rep(0, ncol(tblScores))
Mcov <- var(tblScores)

l <- 280
tictoc::tic()
future_arl(l, k, g, 1e3, T, ar, Ss)
# mean(pbMACusum.RL(20, k, g, 1e4, T, ar, Ss))
tictoc::toc()

# ls <- seq(279,281, 0.5)
# tbl.ARLqiu <- tibble(ls) %>% 
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 5e3, T, ar, Ss)))
# tbl.ARLqiu

# Carta de control con limite y punto donde toca limite
corteLim <- (seq_along(vec.yHwang)[vec.yHwang >= l])[1]

tibble(yh = vec.yHwang, obs = seq_along(vec.yHwang)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  geom_vline(xintercept = verti_line, linetype="dashed") +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Hwang con datos de Suma de Scores',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# Cada cuanto se reinicia la carta? No se reinicio ni una vez

# Calculo de los limites en base a "distribución-maximo"
set.epsi2 <- apply(tblScores, 1, epsilon_h, ar, Ss, T) %>% t

l <- 432.2
tictoc::tic()
future_arl(l, k, g, 1e4, T, ar, Ss, set.epsi2)
tictoc::toc()

# ls <- seq(432,433, 0.2)
# tbl.ARLqiu <- tibble(ls) %>%
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 1e4, T, ar, Ss, set.epsi2)))
# tbl.ARLqiu

# Carta de control con limite y punto donde toca
corteLim <- (seq_along(vec.yHwang)[vec.yHwang >= l])[1]

tibble(yh = vec.yHwang, obs = seq_along(vec.yHwang)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  geom_vline(xintercept = verti_line, linetype="dashed") +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Hwang con datos de Suma de Scores',
       subtitle = 'Ultimo antirango - limite con "distribución maximo"') +
  theme_classic()


# En el ejemplo del articulo utilizaban limite 9.2 suponiendo normalidad.

# rowSums(set.epsi) %>%
#   plot(seq_along(.),., 'l')

# # Con normalidad - La distribución es la misma si sistema de correlaciones.
# Sigma <- var(tblScores)
# 
# sim1 <- MASS::mvrnorm(4000, rep(0,5), Sigma)
# ind1 <- apply(sim1, 1, epsilon_h, ar, Ss, F) %>% t
# colSums(ind1)/sum(colSums(ind1))
# g

# Azucar ------------------------------------------------------------------
nAzuc <- nrow(datosAzuc)

datosAzuc %>% 
  map(summary)

datosAzuc %>% 
  mutate_all(log) %>% 
  map(summary)

pairsAzuc <- ggpairs(datosAzuc)

cor(datosAzuc)

datosAzuc %>% 
  map_dbl(var)

# serie
acf(datosAzuc) # autocorrelacion dentro y entre

acf(datosAzuc[,7])
Box.test(datosSum[,4], 1)

datosAzuc %>% 
  mutate(obs = seq_along(Score7)) %>% 
  gather(score, x, -obs) %>% 
  ggplot(aes(obs, x, colour = score)) +
  geom_line()

datosAzuc %>% 
  mutate(obs = seq_along(Score7)) %>% 
  gather(score, x, -obs) %>% 
  # mutate(x = -log(x)) %>%
  ggplot(aes(obs, x, colour = score)) +
  geom_line()


# histograma
datosAzuc %>% 
  gather(score, x) %>% 
  mutate(x = sqrt(x)) %>%
  # mutate(x = log(x)) %>%
  ggplot(aes(x)) +
  geom_density(aes(fill = score), alpha = 0.5)

colMeans(datosAzuc)

datosAzuc %>% 
  mutate_all(sqrt) %>% 
  map_dfc(~.x - 1) %>% 
  gather(score, x) %>% 
  ggplot(aes(x)) +
  geom_density(aes(fill = score), alpha = 0.5)

datosAzuc %>% 
  mutate_all(sqrt) %>% 
  map_dfc(~.x - 1) %>% 
  acf()


# VAR
datMAtriAzuc <- datosAzuc %>% 
  mutate_all(sqrt) %>% 
  map_dfc(~.x - 1) %>% 
  as.matrix()

selVAR <- vars::VARselect(datMAtriAzuc, lag.max = 10, type ='trend')
selVAR$selection[1] 

fit00 <- vars::VAR(datMAtriAzuc, p=selVAR$selection[1], type="trend")
summary(fit00)
AIC(fit00)

residAzuc <- residuals(fit00) %>% 
  as.tibble()

acf(residAzuc)

residAzuc %>% 
  gather(score, x) %>% 
  ggplot(aes(x)) +
  geom_density(aes(fill = score), alpha = 0.5)

colMeans(residAzuc)

# Arima 
modelos <- datosAzuc %>%
  mutate_all(sqrt) %>% 
  map_dfc(~.x - 1) %>% 
  map(forecast::auto.arima, allowmean = F)
modelos
residAAri <- modelos %>% 
  map_dfc(residuals)

acf(residAAri)

residAAri %>% 
  gather(score, x) %>% 
  ggplot(aes(x)) +
  geom_density(aes(fill = score), alpha = 0.5)

cor(datosAzuc)
cor(residAzuc)
cor(residAAri)

# residAzuc ---------------------------------------------------------------
# qiu ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 7
p <- ncol(residAzuc)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

g <- distr.g(residAzuc, ar, Ss, T)
vec.yQiu <- carta.MACusum(residAzuc, ar, Ss, k, F)

plot.ts(vec.yQiu)

# Distribución de g con datosAzuc
distr.g(datosAzuc, ar, Ss, T)
distr.g(residAAri, ar, Ss, T)

# Simulación limites de control
l <- 20.4
tictoc::tic()
future_arl(l, k, g, 1e4, F, ar, Ss)
tictoc::toc() # 11.85 sec

ls <- seq(19.8,20.5, 0.1)
tbl.ARLqiu <- tibble(ls) %>% 
  mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 5e3, F, ar, Ss)))
tbl.ARLqiu

# Carta de control con limite con corte
corteLim <- (seq_along(vec.yQiu)[vec.yQiu >= l])[1]

tibble(yh = vec.yQiu, obs = seq_along(vec.yQiu)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Qiu con residuales VAR Azucar',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# Cada cuanto se reinicia? No se reinicio ni una vez

# hwang ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 7
p <- ncol(residAzuc)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

g <- distr.g(residAzuc, ar, Ss, T)
vec.yHwang <- carta.MACusum(residAzuc, ar, Ss, k, T)

plot.ts(vec.yHwang)

# Simulación limites de control
Vmu <- rep(0, ncol(residAzuc))
Mcov <- var(residAzuc)

l <- 7.35
tictoc::tic()
future_arl(l, k, g, 5e3, T, ar, Ss)
tictoc::toc()

# ls <- seq(7.2,7.4, 0.1)
# tbl.ARLqiu <- tibble(ls) %>%
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 1e3, T, ar, Ss)))
# tbl.ARLqiu

# Carta de control con limite y punto donde toca limite
corteLim <- (seq_along(vec.yHwang)[vec.yHwang >= l])[1]

tibble(yh = vec.yHwang, obs = seq_along(vec.yHwang)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  geom_vline(xintercept = verti_line, linetype="dashed") +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Hwang con residuales VAR Azucar',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# Cada cuanto se reinicia la carta? No se reinicio ni una vez


# Calculo de los limites en base a "distribución-maximo"
set.epsi2 <- apply(residAzuc, 1, epsilon_h, ar, Ss, T) %>% t

l <- 9.7
tictoc::tic()
future_arl(l, k, g, 5e3, T, ar, Ss, set.epsi2)
tictoc::toc()

# ls <- seq(9.5,9.9, 0.1)
# tbl.ARLqiu <- tibble(ls) %>%
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 1e4, T, ar, Ss, set.epsi2)))
# tbl.ARLqiu

# Carta de control con limite y punto donde toca
corteLim <- (seq_along(vec.yHwang)[vec.yHwang >= l])[1]

tibble(yh = vec.yHwang, obs = seq_along(vec.yHwang)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  geom_vline(xintercept = verti_line, linetype="dashed") +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Hwang con residuales VAR Azucar',
       subtitle = 'Ultimo antirango - limite en base "distribución"') +
  theme_classic()



# residAAri ---------------------------------------------------------------

# qiu ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 7
p <- ncol(residAAri)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

g <- distr.g(residAAri, ar, Ss, T)
vec.yQiu <- carta.MACusum(residAAri, ar, Ss, k, F)

plot.ts(vec.yQiu)

# Simulación limites de control
l <- 22.3
tictoc::tic()
future_arl(l, k, g, 1e4, F, ar, Ss)
tictoc::toc() # 11.85 sec

# ls <- seq(19.8,20.5, 0.1)
# tbl.ARLqiu <- tibble(ls) %>% 
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 5e3, F, ar, Ss)))
# tbl.ARLqiu

# Carta de control con limite con corte
corteLim <- (seq_along(vec.yQiu)[vec.yQiu >= l])[1]

tibble(yh = vec.yQiu, obs = seq_along(vec.yQiu)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_hline(yintercept = l, colour = 'blue') +
  geom_vline(xintercept = corteLim, colour = 'red') +
  labs(x='Observaciones', y='Estadística monitoreo', 
       title = 'Carta Qiu con residuales ARIMA Azucar',
       subtitle = 'Ultimo antirango') +
  theme_classic()

# Cada cuanto se reinicia? No se reinicio ni una vez










# cor1 <- c(1,.9,.7,
#           .9,1,.4,
#           .7,.4,1) %>% 
#   matrix(3, 3)
# 
# dat1 <- MASS::mvrnorm(n = 400, 
#                       rep(0, 3), 
#                       cor1)  
# cor(dat1)
# acf(dat1)

# library(vars)
# library(astsa)
# x = cbind(cmort, tempr, part)
# plot.ts(x , main = "", xlab = "")
# fitvar1=VAR(x, p=1, type="both")
# summary(fitvar1)