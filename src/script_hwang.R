#---- 
# file: script_hwang.R
# author: CU
# date: 17-06-18
# description: Replicación del paper de HWANG(2016)
# input: 
# output:
#    - 
#    - 
#---
rm(list = ls())

# packages ----------------------------------------------------------------
library(tidyverse)
library(listenv)
# library(furrr)
# library(future)

# paths -------------------------------------------------------------------
inPath <- file.path('input')
outPath <- file.path('output')

# read data ---------------------------------------------------------------
global <- listenv()

inDatos <- file.path(inPath, "datos.xlsx")

datos2 <- openxlsx::read.xlsx(inDatos, 'tabla2')


# functions ---------------------------------------------------------------
# mvrnorm <- MASS::mvrnorm
# perm <- function(n,k){choose(n,k) * factorial(k)}
# 
# Qt <- function(S1_1,S2_1,et,g1=g){
#   A <- diag((S2_1+g1)^(-1))
#   a <- (S1_1-S2_1)+(et-g1)
#   return(t(a)%*%A%*%a %>% as.numeric())
# }
# 
# # Yt <- function(S1,S2){
# #   A <- diag(S2^(-1))
# #   a <- S1-S2
# #   return(t(a)%*%A%*%a %>% as.numeric())
# # }
#   
# Sh.12 <- function(S1_1,S2_1,et,Qt,g1=g,k1=k){
#   if(Qt<=k1){
#     S1 <- S2 <- S1_1*0
#   }else if(Qt>k1){
#     S1 <- ( (S1_1+et)*(Qt-k1) )/(Qt)
#     S2 <- ( (S2_1+g1)*(Qt-k1) )/(Qt)
#   }else{
#     message('Qt revisar')
#   }
#   return(list(S1=S1,S2=S2))
# }
# 
# # Función para el calculo del RL de las dos cartas
# 
# RL2.pb <- function(i, alter =F){
#   setTxtProgressBar(pb, i)
#   RL <- 0
#   flag <- T
#   s1 <- S1.0
#   s2 <- S2.0
#   
#   while(flag){
#     RL <- RL + 1
#     
#     if(!alter){
#       e1t <- rmultinom(1, 1, g)[,1]
#       # obs <- mvrnorm(1, Vmu, Mcov)
#       # e1t <- epsilon_h(obs, vec, Ss2.4)
#     }else{
#       obs <- mvrnorm(1, Vmu, Mcov)
#       e1t <- epsilon_hC(obs, vec, Ss2.4)
#     }
#     qt <- Qt(s1,s2,e1t)
#     # S1_1,S2_1,et,Qt,g1=g,k1=k
#     s12 <- Sh.12(s1,s2,e1t,qt)
#     s1 <- s12$S1
#     s2 <- s12$S2
#     yt <- max(0,qt - k) # 
#     flag <- yt < L
#   }
#   return(RL)
# }
# 
# # Funciones para generar el vector epsilon para la carta 
# # este vector esta en base al conjunto S 
# epsil_h <- function(uno, Ss) apply(Ss, 1, function(x) ifelse(identical(x,uno),1,0))
# epsil_h <- function(uno, Ss) apply(Ss, 1, function(x) ifelse(all(x == uno),1,0))
# 
# # se puede borrar la función epsil_h, meterla en las dos de abajo
# epsilon_h <- function(obs, vec_ep, Ss=NULL){
#   # Original
#   # saca el vector epsilon de S, en base a la observación de una muestra
#   ep <- order(obs)[vec_ep]
#   return(epsil_h(ep, Ss))
# }
#
#epsilon_hC <- function(obs, vec_ep, Ss=NULL){
#  # saca el vector epsilon de S, en base a la propuesta de hwang
#  # ep <- order(obs)[vec_ep]
#  # vep_op <- epsil_h(ep, Ss)
#  vep_op <- epsilon_h(obs,vec_ep, Ss)
#  
#  if(1 %in% vec_ep){
#    return(vep_op*abs(min(obs)))
#  }else if(length(obs) %in% vec_ep){
#    return(vep_op*abs(max(obs)))
#  }else{
#    message('Problemas de dimensión entre la observacion y epsilon')
#  }
#}

# Funciones para el calculo de la estadistica de control 
# Generación aleatoria de la estadistica
yt.pb <- function(i){
  setTxtProgressBar(pb, i)
  s1 <- global$s1
  s2 <- global$s2
  
  e1t <- rmultinom(1, 1, g)[,1] 
  qt <- Qt(s1,s2,e1t)
  s12 <- Sh.12(s1,s2,e1t,qt)
  global$s1 <- s12$S1
  global$s2 <- s12$S2
  yt <- qt - k
  
  return(yt)
}

yt2.pb <- function(i, alter = F){
  setTxtProgressBar(pb, i)
  s1 <- global$s1
  s2 <- global$s2
  
  if(!alter){
    # e1t <- rmultinom(1, 1, g)[,1]
    obs <- mvrnorm(1, Vmu, Mcov)
    e1t <- epsilon_h(obs, vec, Ss2.4)
  }else{
    obs <- mvrnorm(1, Vmu, Mcov)
    e1t <- epsilon_hC(obs, vec, Ss2.4)
  }
  
  qt <- Qt(s1,s2,e1t)
  s12 <- Sh.12(s1,s2,e1t,qt)
  global$s1 <- s12$S1
  global$s2 <- s12$S2
  yt <- qt - k
  
  return(yt)
}

# Funciones para calcular la estadistica de control con datos observados

# est.Yt <- function(e,S1,S2, show.q = F){
#   qt <- Qt(s1,s2,e1t)
#   
#   if(!show.q){
#     return(qt - k)
#   }else{
#     return(list(y = qt - k, q =qt))
#   }
# }

# Función para calcular el vector de epsilon con un conjunto de datos (o dar un conjunto aleatorio)
# esta sirve para calcular el vector de g-distribución de epsilon
M.epsilon <- function(datos, vec_epsi, Ss=NULL){
  # incluir comprobación de dimension columnas Ss == tamaño de vec_epsi
  # comprobar si es tamaño 1 vec entonces SS debe ser null

  antirangos <- apply(datos,1,order) %>% t
  epsi <- antirangos[,vec_epsi]

  if(length(vec_epsi) == 1){
    cat('Epsilon de tamaño uno \n')
    tibble(epsi,indi = 1, id =seq_along(epsi)) %>%
      spread(epsi,indi) %>%
      select(-id) %>%
      mutate_all(~ifelse(is.na(.),0,1)) %>%
      as.matrix() %>%
      return()
  }else if(length(vec_epsi) > 1){
    cat('Epsilon de tamaño ',length(vec_epsi),'\n')
    apply(epsi, 1,epsil_h,Ss) %>% t %>% return()
  }else{
    message('Algún problema con el tamaño de Epsilon \n')
  }
}

# Para revisar la cota del parametro k
val12 <- function(i,g) sum(g[-i])/g[i]
cota.k <- function(g){
  return(sapply(seq_along(g), val12, g))
}

# epsilon = 1 -------------------------------------------------------------
# ORIGINAL
p <- 4
S1.0 <- S2.0 <- numeric(p)
g <- rep(1/p,p)

# k <- 0.5
# L <- 8.029
k <- 1
L <- 6.842

I <- 1e3
pb <- txtProgressBar(max=I,style = 3,char='-')
RLs <- sapply(1:I, RL2.pb)

mean(RLs)

# PROPUESTA
p <- 4
q <- 1
Mcov <- diag(p)
Vmu <- numeric(p)
g <- rep(1/p,p)
vec <- 1:q
S1.0 <- S2.0 <- numeric(p)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 0.5
L <- 19.23

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='-')
RLs <- sapply(1:I, RL2.pb, T)

mean(RLs) # 199.8

# para epsilon = 2 --------------------------------------------------------
# ORIGINAL
p <- 4
q <- 2
dg <- perm(p,q)
S1.0 <- S2.0 <- numeric(dg)
g <- rep(1/dg, dg)

k <- 1
L <- 15.6887

I <- 5e3
pb <- txtProgressBar(max=I,style = 3,char='+')
RLs <- sapply(1:I, RL2.pb)

mean(RLs)
summary(RLs)

# PROPUESTA
p <- 4
q <- 2
Mcov <- diag(p)
Vmu <- numeric(p)
vec <- 1:q
dg <- perm(p,q)
S1.0 <- S2.0 <- numeric(dg)
g <- rep(1/dg, dg)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 1
L <- 36

I <- 5e3
pb <- txtProgressBar(max=I,style = 3,char='-')
RLs <- sapply(1:I, RL2.pb, T)

mean(RLs) 


# para epsilon = 3 --------------------------------------------------------
# ORIGINAL
p <- 4
q <- 3
dg <- perm(p,q)
S1.0 <- S2.0 <- numeric(dg)
g <- rep(1/dg, dg)

k <- 1
L <- 24.845

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='+')
RLs <- sapply(1:I, RL2.pb)

# 164, 167
mean(RLs)
summary(RLs)

# PROPUESTA
p <- 4
q <- 3
Mcov <- diag(p)
Vmu <- numeric(p)
vec <- 1:q
dg <- perm(p,q)
S1.0 <- S2.0 <- numeric(dg)
g <- rep(1/dg, dg)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 1
L <- 50

I <- 5e3
pb <- txtProgressBar(max=I,style = 3,char='-')
RLs <- sapply(1:I, RL2.pb, T)

mean(RLs) 


# calculo del nivel alpha -------------------------------------------------
# ORIGINAL
p <- 4
g <- rep(1/p,p)

k <- 0.5
L <- 8.029

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='+')
global$s1 <- global$s2 <- numeric(p)
ys <- sapply(1:I, yt.pb)

sum(ys >= L)/I # alpha

tibble(y=ys) %>% 
  ggplot(aes(y)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = L) +
  labs(title='Histograma de la estadistica de monitoreo en control',
       subtitle = 'Alternativa original')

sum(ys < 0)

# PROPUESTA
p <- 4
q <- 1
Mcov <- diag(p)
Vmu <- numeric(p)
g <- rep(1/p,p)
vec <- 1:q
# S1.0 <- S2.0 <- numeric(p)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 0.5
L <- 19.23

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='+')
global$s1 <- global$s2 <- numeric(p)
ys <- sapply(1:I, yt2.pb, T)

sum(ys >= L)/I # alpha

tibble(y=ys) %>% 
  ggplot(aes(y)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = L) +
  labs(title='Histograma de la estadistica de monitoreo en control',
       subtitle = 'Alternativa propuesta')



# power -------------------------------------------------------------------
# ORIGINAL
p <- 4
q <- 1
g <- rep(1/p,p)
vec <- 1:q
delta <- 0.5
Mcov <- diag(p)
Vmu <- numeric(p) + delta*c(1,0,0,0)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 0.5
L <- 8.029

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='+')
global$s1 <- global$s2 <- numeric(p)
ys <- sapply(1:I, yt2.pb)

sum(ys >= L)/I # power

tibble(y=ys) %>% 
  ggplot(aes(y)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = L) +
  labs(title='Histograma de la estadistica de monitoreo fuera de control',
       subtitle = 'Alternativa original')

sum(ys < 0)

# PROPUESTA
p <- 4
q <- 1
Mcov <- diag(p)
delta <- -0.5
Vmu <- numeric(p) + delta*c(1,0,0,0)
g <- rep(1/p,p)
vec <- 1:q
S1.0 <- S2.0 <- numeric(p)

Ss2.4 <- gtools::permutations(p,q,1:p)

k <- 0.5
L <- 19.23

I <- 1e4
pb <- txtProgressBar(max=I,style = 3,char='+')
global$s1 <- global$s2 <- numeric(p)
ys <- sapply(1:I, yt2.pb, T)

sum(ys >= L)/I # power

tibble(y=ys) %>% 
  ggplot(aes(y)) +
  geom_histogram(binwidth = 0.1) +
  geom_vline(xintercept = L) +
  labs(title='Histograma de la estadistica de monitoreo en control',
       subtitle = 'Alternativa propuesta')


# datos 2 -----------------------------------------------------------------

tbl2 <- datos2 %>% 
  select(-Obs)

apply(tbl2,1,epsilon_hC, 1:1, Ss2.4) %>% t



# vector g por muestras ---------------------------------------------------
p <- 3
q <- 1

Mcov <- c(1,0.903,0.867,0.903,1,0.864,0.867,0.864,1) %>% 
  matrix(3)
Mcov <- diag(p)

Vmu <- numeric(p)
n <- 1e6

muestra <- MASS::mvrnorm(n,Vmu,Mcov)

# Si estamos intereados en epsilon 1:3
Ss3.4 <- gtools::permutations(p,q,1:p)
epsi.13 <- M.epsilon(muestra,1:3,Ss3.4)
colMeans(epsi.13)
summary(colMeans(epsi.13))

# Si estamos intereados en epsilon 1:2
Ss2.4 <- gtools::permutations(p,q,1:p)
epsi.12 <- M.epsilon(muestra,1:2,Ss2.4)
colMeans(epsi.12)

# Si solo estamos interesados en epsilon 1
epsi.1 <- M.epsilon(muestra,1)
colMeans(epsi.1)
mean(colMeans(epsi.1))
