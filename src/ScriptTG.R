#---- 
# file: ScriptTG.R
# author: CU
# date:
# description: Ejemplos desarrollados para el TG
# input: Datos
# output:
#    - Graficas y resultados correspondientes de los ejemplos
#    - 
#---
rm(list = ls())

# packages ----------------------------------------------------------------
library(tidyverse)
library(ggfortify)
library(furrr)
library(tictoc)
# require(pbapply) # for progress bar
# library(listenv)
# library(future)

# paths -------------------------------------------------------------------
inPath <- file.path('input')
inFunc <- file.path('src','func')
outPath <- file.path('output')


# functions ---------------------------------------------------------------
pathFunc <- file.path(inFunc, 'funciones.R')
source(pathFunc, encoding = 'UTF8')


# read data ---------------------------------------------------------------
# global <- listenv::listenv()

# Datos de la tabla 1 del articulo de Hwang
inDatos <- file.path(inPath, "datos.xlsx")
datos1 <- openxlsx::read.xlsx(inDatos, 'tabla1')

# Datos de la tabla 2 del articulo de Hwang
inDatos <- file.path(inPath, "datos.xlsx")
datos2 <- openxlsx::read.xlsx(inDatos, 'tabla2')

# Datos del ejemplo 9.4 del libro de Qiu
# inResid <- file.path(inPath, 'soilResidual.txt')
# tbl.alum <- read.table(inResid) %>% 
#   tbl_df()

# data(Alsmelterdata) # library(NPMVCP)

# Revisión calculo de la estadística --------------------------------------
tbl00 <- datos1 %>% 
  select(-Obs)

# Conjunto S para epsilon_1
p <- 4
q <- 1
Ss1 <- gtools::permutations(p,q,1:p)

# Parametro de la carta
k1 <- 0.5
g1 <- rep(1/4,4)
L_qiu <- 8.029

carta.MACusum(tbl00, c(1), Ss1, k1, hwang = F, g = g1) %>% 
  round(.,2) # Están igual que en el articulo. Carta Qiu

# Carta Hwang
tbl01 <- datos2 %>% 
  select(-Obs)

L_qiu <- 19.23

carta.MACusum(tbl01, c(1), Ss1, k1, hwang = T, g = g1) %>% 
  round(.,2) # Están igual que en el articulo. Carta Hwang

# Por el momento el calculo de la estadística de control se encuentra bien.


# - - - - - - - - - - - Simulaciones - - - - - - - - - - - - - - - - - - - 

# Simulación de tabla 3 hwang ---------------------------------------------

k1 <- 1

# Apartir de una normal multivariada de dimensión 4
mu_0 <- rep(0,4)
sigma_0 <- diag(rep(1,4))

# Se realizaran unos downward shift levels en X_1 de
shift <- c(0,-0.5,-1,-2,-3,-5,-10)

mu_list <- list()
for (sh in shift) {
  mu_aux <- mu_0
  mu_aux[1] <- sh
  mu_list <- append(mu_list, list(mu_aux))
}

# Y se consideran los vectores de antirangos
ar1 <- 1
ar12 <- 1:2
ar123 <- 1:3
ARs <- list(ar1,ar12,ar123)

# Los limites de control para cada carta serían
# Ls_qiu <- c(6.842, 15.6887, 24.845)
# Ls_hwang <- c(24, 36, 50)
Ls_qiu <- c(6.842, 15.6887, 25.31)
Ls_hwang <- c(24, 36, 50.1)

tic1 <- tic()
tbl3_qiu <- tibble(ar = rep(ARs,length(mu_list)),
                l = rep(Ls_qiu, length(mu_list)),
                Vmu = rep(mu_list, each = length(ARs))) %>%
  mutate(arl = pmap(list(Vmu,ar,l), future_OCarl.MAC, sigma_0 ,k1, 1e4))
toc1 <- toc()
tiempo1 <- paste0( round(toc1$toc -toc1$tic, 2), ' segundos')
print(tiempo1) # 1186.22 segundos


tic2 <- tic()
tbl3_hwang <- tibble(ar = rep(ARs,length(mu_list)),
                   l = rep(Ls_hwang, length(mu_list)),
                   Vmu = rep(mu_list, each = length(ARs))) %>%
  mutate(arl = pmap(list(Vmu,ar,l), future_OCarl.MAC, sigma_0 ,k1, 1e4, T))
toc2 <- toc()
tiempo2 <- paste0( round(toc2$toc -toc2$tic, 2), ' segundos')
print(tiempo2) # 846.53 segundos

# Tiempo en total aproximado de 33.88 minutos.

tbl3.0 <- tbl3_qiu %>% 
  mutate(carta = 'qiu') %>% 
  rbind(
    tbl3_hwang %>% 
      mutate(carta = 'hwang')
  )

# outTbl3 <- file.path(outPath, 'tbl3Hwang.rdata')
# save(tbl3.0, file = outTbl3)
# load(outTbl3)

tbl3.1 <- tbl3.0 %>% 
  mutate(arlNames = map(arl, names),
         shift = map_dbl(Vmu, ~round(.x[1], 1)),
         Vmu = map_chr(Vmu, ~paste0('(',paste0(.x, collapse = ','),')')) ,
         ar = map_chr(ar,~paste0('AR',paste0(.x, collapse = '')))
         ) %>% 
  select(ar, arl, arlNames, carta, shift, Vmu)

tbl3.1 %>% 
  unnest() %>% 
  unnest() %>% 
  spread(arlNames, arl) %>% 
  mutate(sdARL = sdRL/1e4 ) 

tbl3.2 <- tbl3.1 %>% 
  unnest() %>% 
  unnest() %>% 
  spread(arlNames, arl) %>% 
  mutate(sdARL = sdRL/1e4,
         ARL_sd = paste0(sprintf("%.2f", ARL),' (',sprintf("%.4f", sdARL),')') 
         ) %>% 
  select(-sdRL, -ARL, -sdARL) %>% 
  spread(carta,ARL_sd) %>% 
  select(ar, shift, Vmu, qiu, hwang) %>% 
  arrange(ar, desc(shift )) %>% 
  select(-shift)

library(knitr)
library(kableExtra)

tbl3.2 %>% 
  kable(., "latex", booktabs = T, align = "c") %>%
  collapse_rows(columns = 1, latex_hline = "major")

# Limites para AR123 ------------------------------------------------------
# Vmu <- rep(0,4)
# Mcov <- diag(rep(1,4))
# 
# # Conjunto S para epsilon_1
# p <- 4
# q <- 3
# Ss123 <- gtools::permutations(p,q,1:p)
# 
# # Parametro de la carta
# k1 <- 1
# g1 <- numeric(perm(p,q)) + 1/perm(p,q)
# L123_qiu <- 24.845
# L123_hwang <- 50
# # AR123 0 172.93(0.0282) 186.43(0.0356)
# 
# # Carta de qiu, buscando el L para AR123 tenga un ARL=200
# list123_qiu <- c(25.31) # 201
# 
# tibble(Lim = list123_qiu) %>% 
#   mutate(arl = map_dbl(Lim, future_arl, k1, g1, 3e4, F, 1:3, Ss123))
# 
# # Carta de hwang, buscando el L para AR123 tenga un ARL=200
# list123_hwang <- c(50.1) # 193
#   
# tibble(Lim = list123_hwang) %>% 
#   mutate(arl = map_dbl(Lim, future_arl, k1, g1, 1e4, T, 1:3, Ss123))

# # A tibble: 3 x 2
# Lim   arl
# <dbl> <dbl>
# 1  51.0   216
# 2  52.0   250
# 3  53.0   294

# Simulación de tabla 4 hwang ---------------------------------------------

k1 <- 1

# Apartir de una normal multivariada de dimensión 4
# mu_0 <- rep(0,4)
sigma_0 <- diag(rep(1,4))

# Se realizaran unos downward shift levels en todas las direcciones
sh1 <- c(-0.5, -4, 4)
sh2 <- c(-2,-4,-4,-2,-3,-4,-4)
sh3 <- c(0,0,-2,-2,-1,-4,-4)
sh4 <- c(rep(0,6),-4)

mu_list <- tibble(v1=rep(sh1,each=7), 
       v2=rep(sh2,3),
       v3=rep(sh3,3),
       v4=rep(sh4,3) ) %>% 
  t() %>% 
  as.tibble() %>% 
  as.list()


# Y se consideran los vectores de antirangos
ar1 <- 1
ar12 <- 1:2
ar123 <- 1:3
ARs <- list(ar1,ar12,ar123)

# Los limites de control para cada carta serían
# Ls_qiu <- c(6.842, 15.6887, 24.845)
# Ls_hwang <- c(24, 36, 50)
Ls_qiu <- c(6.842, 15.6887, 25.31)
Ls_hwang <- c(24, 36, 50.1)



tic14 <- tic()
tbl4_qiu <- tibble(ar = rep(ARs,length(mu_list)),
                   l = rep(Ls_qiu, length(mu_list)),
                   Vmu = rep(mu_list, each = length(ARs))) %>%
  mutate(arl = pmap(list(Vmu,ar,l), future_OCarl.MAC, sigma_0 ,k1, 1e4))
toc14 <- toc()
tiempo14 <- paste0( round(toc14$toc -toc14$tic, 2), ' segundos')
print(tiempo14) # "2746.13 segundos"


tic24 <- tic()
tbl4_hwang <- tibble(ar = rep(ARs,length(mu_list)),
                     l = rep(Ls_hwang, length(mu_list)),
                     Vmu = rep(mu_list, each = length(ARs))) %>%
  mutate(arl = pmap(list(Vmu,ar,l), future_OCarl.MAC, sigma_0 ,k1, 1e4, T))
toc24 <- toc()
tiempo24 <- paste0( round(toc24$toc -toc24$tic, 2), ' segundos')
print(tiempo24) # "294.06 segundos"

# Tiempo en total aproximado de 33.88 minutos.

tbl4.0 <- tbl4_qiu %>% 
  mutate(carta = 'qiu') %>% 
  rbind(
    tbl4_hwang %>% 
      mutate(carta = 'hwang')
  )

# outTbl4 <- file.path(outPath, 'tbl4Hwang.rdata')
# save(tbl4.0, file = outTbl4)
load(outTbl4)

tbl4.1 <- tbl4.0 %>% 
  mutate(arlNames = map(arl, names),
         shift = map_dbl(Vmu, ~round(.x[1], 1)),
         Vmu = map_chr(Vmu, ~paste0('(',paste0(.x, collapse = ','),')')) ,
         ar = map_chr(ar,~paste0('AR',paste0(.x, collapse = '')))
  ) %>% 
  select(ar, arl, arlNames, carta, shift, Vmu)


tbl4.2 <- tbl4.1 %>% 
  unnest() %>% 
  unnest() %>% 
  spread(arlNames, arl) %>% 
  mutate(sdARL = sdRL/1e4,
         ARL_sd = paste0(sprintf("%.2f", ARL),' (',sprintf("%.2f", sdARL),')'),
         nameUp = paste(ar, carta, sep = '_')
  ) %>% 
  select(-sdRL, -ARL, -sdARL, -ar, -carta, -shift) %>% 
  spread(nameUp,ARL_sd) %>%   
  mutate(Vmu = str_replace_all(Vmu,'\\(|\\)','')) %>% 
  separate(Vmu, paste0('V',1:4), ',') %>% 
  mutate_at(vars(starts_with('V')), funs(as.double)) %>% 
  arrange(desc(V1), desc(V4), desc(V3), desc(V2)) %>% 
  mutate(Vmu = paste0('(',paste(V1,V2,V3,V4, sep = ','),')')) %>% 
  select(Vmu, starts_with('AR'))

tbl4.2 %>% 
  kable(.,"latex", booktabs = T, align = "r") %>% 
  add_header_above(c(" ", "AR1" = 2, "AR12" = 2, "AR123" = 2)) 



# tbl3.2 %>% 
#   kable(., "latex", booktabs = T, align = "c") %>%
#   collapse_rows(columns = 1, latex_hline = "major")


# simulation t-multi ------------------------------------------------------
# rmvt <- mvtnorm::rmvt
# rmvnorm <- mvtnorm::rmvnorm
# 
# ## X ~ t_3(0, diag(2))
# x <- rmvt(1000, sigma = diag(2), df = 5) # t_3(0, diag(2)) sample
# plot(x)
# 
# sigma <- matrix(c(4,2,2,3), ncol=2)
# x <- rmvnorm(n=1000, mean=c(0,0), sigma=sigma, method="chol")
# plot(x)


# holmes ------------------------------------------------------------------

inHolmes <- file.path(inPath, 'HolmesData.xlsx')
dataHol <- openxlsx::read.xlsx(inHolmes) %>% 
  mutate_at(vars(-numberSample),~..1/10) %>% 
  select(-numberSample)

cor(dataHol)

# Extraer la media
vec_media <- colMeans(dataHol)
vec_media

dataHolCen <- scale(dataHol,T,F) %>% as.tibble()

dataHolCen <- scale(dataHol) %>% as.tibble()

# dataHolCen %>% 
#   as.tibble() %>% 
#   gather(var, y) %>% 
#   ggplot(aes(y)) +
#   geom_density(aes(fill = var), alpha = 0.5)

# qiu ----------------------------------
datCont <- dataHolCen %>% 
  mutate(M = ifelse(seq_along(M) %in% 11:56, M - 3, M)) 


# Parametros de la carta
k <- 0.5
ar <- 1
p <- ncol(dataHolCen)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

# g <- distr.g(dataHolCen, ar, Ss, T)
g <- rep(0.5, 2)
vec.yQiu <- carta.MACusum(datCont, ar, Ss, k, F, g)

plot.ts(vec.yQiu)

# Simulación limites de control
l <- 3.5
tictoc::tic()
future_arl(l, k, g, 5e3, F, ar, Ss)
tictoc::toc()

# Carta de control con limite
corteLim <- (seq_along(vec.yQiu)[vec.yQiu >= 3.5])[1]

tibble(yh = vec.yQiu, obs = seq_along(vec.yQiu)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = 3.5, linetype="dashed") +
  labs(x='Observación', y='Estadística monitoreo', 
       title = 'Carta MA-CUSUM',
       subtitle = 'Primer antirango') +
  theme_classic()

# hwang ----------------------------------

# Parametros de la carta
k <- 0.5
ar <- 1
p <- ncol(dataHolCen)
q <- length(ar)
Ss <- gtools::permutations(p,q,1:p)

# g <- distr.g(tblScores, ar, Ss, T)
g <- rep(0.5, 2)
vec.yHwang <- carta.MACusum(datCont, ar, Ss, k, T, g)

plot.ts(vec.yHwang)

# Simulación limites de control
Vmu <- rep(0, ncol(datCont))
Mcov <- diag(2)

l <- 9.2
tictoc::tic()
future_arl(l, k, g, 1e3, T, ar, Ss)
# mean(pbMACusum.RL(20, k, g, 1e4, T, ar, Ss))
tictoc::toc()

# ls <- seq(279,281, 0.5)
# tbl.ARLqiu <- tibble(ls) %>% 
#   mutate(ARL = map_dbl(ls,~future_arl(.x, k, g, 5e3, T, ar, Ss)))
# tbl.ARLqiu

# Carta de control con limite y punto donde toca limite
corteLim <- (seq_along(vec.yHwang)[vec.yHwang >= 9.2])[1]

tibble(yh = vec.yHwang, obs = seq_along(vec.yHwang)) %>% 
  ggplot(aes(obs, yh)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = 9.2, linetype="dashed") +
  labs(x='Observación', y='Estadística monitoreo', 
       title = 'Carta HMA-CUSUM',
       subtitle = 'Primer antirango') +
  theme_classic()








