#---- 
# file: Script922.R
# author: CU
# date:
# description: Replicación de los ejemplos 9.4 y 9.5
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
outPath <- file.path('output')

# read data ---------------------------------------------------------------
# global <- listenv::listenv()

inResid <- file.path(inPath, 'soilResidual.txt')
tbl.alum <- read.table(inResid) %>% 
  tbl_df()

# data(Alsmelterdata) # library(NPMVCP)

# functions ---------------------------------------------------------------
perm <- function(n,k){choose(n,k) * factorial(k)}

# Funciones necesarias para calcular Qh,S1,S2 para el algoritmo
Qt <- function(S1_1,S2_1,et,g1=g) sum(((S1_1-S2_1)+(et-g1))^2/(S2_1+g1))

# S1 observada, S2 esperada
Sh.12 <- function(S1_1,S2_1,et,Qt,g1=g,k1=k){
  if(Qt<=k1){
    S1 <- S2 <- S1_1*0
  }else if(Qt>k1){
    S1 <- ( (S1_1+et)*(Qt-k1) )/(Qt)
    S2 <- ( (S2_1+g1)*(Qt-k1) )/(Qt)
  }else{
    message('Qt revisar')
  }
  return(list(S1=S1,S2=S2))
}

# Función para el calculo del RL de las dos cartas
rl.MACusum <- function(i, L, k, g, hwang=F, ar = NULL, Ss=NULL){
  # Descripción :
  # Esta función realiza el calculo de un RL bajo la carta MA-Cusum.
  # Necesita las funciones Qt y Sh.12. Y epsilon_h cuando hwang=TRUE
  #
  # Parametros :
  # i := Es el i-esimo RL, esto es para implementar map, o barra de progreso.
  # L := es el limite de la carta
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # g := distribución del antirango 
  # hwang := Si es TRUE utiliza la carta de hwang, con una normal-mult.
  #           Necesario que esten definido las variables Vmu, Mcov
  # ar := Cuales son los antirangos tomados por la carta de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # Ss := Es una matriz de con  el espacio muestral del antirango [ar],
  #       necesaria para la carta hwang.
  #
  # Output :
  # Un entero que es el calculo del RL.

  
  # setTxtProgressBar(pb, i)
  RL <- 0
  flag <- T
  s1 <- s2 <- numeric(length(g))
  
  while(flag){
    RL <- RL + 1
    
    if(!hwang){
      e1t <- rmultinom(1, 1, g)[,1]
    }else{ # falta arreglar este pedazo por el conjunto Ss
      if(is.null(a) | is.null(b)) stop('Arguments ar o Ss do not define')
      if(!(exists('Vmu') & exists('Mcov'))) stop('Variables Vmu & Mcov need to be created in the GlobalEnv')
      obs <- MASS::mvrnorm(1, Vmu, Mcov)
      e1t <- epsilon_h(obs, ar, Ss, hwang = T)
    }
    
    qt <- Qt(s1,s2,e1t,g)
    s12 <- Sh.12(s1,s2,e1t,qt,g,k)
    s1 <- s12$S1
    s2 <- s12$S2
    yt <- max(0,qt - k) # Qiu,Hawkins(2001) equación (8)
    flag <- yt < L
  }
  return(RL)
}

pbMACusum.RL <- function(L, k, g, n.iter=1e3, hwang = F, ar = NULL, Ss=NULL){
  # Descripción :
  # Calculo bajo simulación n.iter longitudes de corrida,
  # cuenta con una barra de progreso (pb)
  # Esta función necesita la función rl.MACusum.
  #
  # Parametros :
  # L := es el limite de la carta
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # g := distribución del antirango
  # n.iter := numero de iteraciones
  # hwang := Si es TRUE utiliza la carta de hwang, con una normal-mult.
  #           Necesario que esten definido las variables Vmu, Mcov
  # ar := Cuales son los antirangos tomados por la carta de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # Ss := Es una matriz de con  el espacio muestral del antirango [ar],
  #       necesaria para la carta hwang.
  #
  # Output :
  # Un vector de RL's con longitud n.iter
  
  # pb <- txtProgressBar(max=I,style = 3,char='-')
  RLs <- pbapply::pbsapply(1:n.iter, rl.MACusum,L, k, g,hwang,ar,Ss)
  return(RLs)
}



future_arl <- function(L, k, g, n.iter=1e3, hwang = F, ar = NULL, Ss=NULL){
  # Descripción : 
  # Calculo bajo simulación de un ARL
  # Reparte el calculo de cada ciclo, por los nucleos, (Paralelizar)
  # para tomar menos tiempo en el calculo de todos los RL's
  # Se necesita la función rl.MACusum.
  #
  # Parametros :
  # L := es el limite de la carta
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # g := distribución del antirango
  # n.iter := numero de iteraciones para el calculo del ARL.
  # hwang := Si es TRUE utiliza la carta de hwang, con una normal-mult.
  #           Necesario que esten definido las variables Vmu, Mcov
  # ar := Cuales son los antirangos tomados por la carta de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # Ss := Es una matriz de con  el espacio muestral del antirango [ar],
  #       necesaria para la carta hwang.
  #
  # Output :
  # Un real que es el ARL calculado.
  
  plan(multiprocess)
  # plan(multiprocess, workers = 3)
  RLs <- future_map_dbl(1:n.iter, rl.MACusum,L, k, g, hwang,ar,Ss)
  cat('Finished the ARL calculation with L = ', L, '\n')
  return(mean(RLs))
}



# Funciones usadas para calcular los vectores epsilon de indicadoras para la carta
#epsil_h <- function(uno, Ss) apply(Ss, 1, function(x) ifelse(all(x == uno),1,0))
epsilon_h <- function(obs, ar, Ss=NULL, hwang = F){
  # Descripción : 
  # Retorna el vector epsilon sobre el espacio muestral S,
  # en base de un vector de datos observados o simulado.
  #
  # Parametros : 
  # obs := Vector observado, el cual se obtiene el antirango [ar]
  # ar := Es el vector que dice cuales antirangos se tienen en cuenta pa la carta
  # Ss := Es una matriz de perm(length(obs), length(ar)) por length(ar), con 
  #       el espacio muestral del antirango [ar]
  # hwang := Si es TRUE utiliza la definición de epsilon de hwang.
  #
  # Output : 
  # Un vector indicadora para el calculo de la carta
  
  ep <- order(obs)[ar]
  vep_op <- apply(Ss, 1, function(x) ifelse(all(x == ep),1,0)) # epsil_h(ep, Ss)
  
  if(!hwang){
    return(vep_op)
  }else{ 
    if(1 %in% vec_ep){
      return(vep_op*abs(min(obs)))
    }else if(length(obs) %in% vec_ep){
      return(vep_op*abs(max(obs)))
    }else{
      stop('For hwang chart is necessary that arg(ar) must been defined with the first or last antirank')
    }
  }
}

# Función para calculo de la carta desde historico
carta.MACusum <- function(dataSet, ar, Ss, k, hwang = F, g = NULL){
  # Descripción : 
  # Desde un conjunto de datos calcula la estadistica de la carta MA-Cusum
  # Necesita las funciones epsilon_h, Qt, Sh.12
  # dataSet := Datos observados dimensión nxp, siendo
  #             n la cantidad de muestras y
  #             p cantidad de variables observadas.
  # ar := Cuales son los antirangos tomados por la carta, si es la de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # Ss := Es una matriz con  el espacio muestral del antirango [ar]
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # hwang := Si es TRUE utiliza la carta de hwang.
  # g := distribución del antirango, si es NULL se calcula de los datos

  # Se necesita las funciones Qt y Sh.12
  
  set.epsi <- apply(dataSet, 1, epsilon_h, ar, Ss, hwang) %>% t
  
  if(is.null(g)){ # este if se encarga de generar la distribución del [ar] de los datos
    g <- colSums(set.epsi)/sum(colSums(set.epsi))
    if( !all(g != 0) ){
      # Realizo esta correción de ceros en la distribución del vector [ar]
      # porque considero que es poco probable que sea igual a cero 
      # sino que es poco probable, con una probabilidad pequeña
      # Qué valor dejarle en vez de 0.5?
      warning(paste0('There was a correction in vector g of ', sum(g == 0), 
                     ' elements in 0 by 0.5'))
      conteo <- colSums(set.epsi)
      conteo[conteo == 0] <- 0.5
      g <- conteo/sum(conteo)
    }
  }
  s1 <- s2 <- numeric(length(g))
  
  vec.y <- numeric(nrow(dataSet))
  for (ii in 1:nrow(dataSet)) {
    e1t <- set.epsi[ii,]
    qt <- Qt(s1,s2,e1t,g)
    #if(ii == 1) qt <- 0 # si es como example88.r??
    s12 <- Sh.12(s1,s2,e1t,qt,g,k)
    s1 <- s12$S1
    s2 <- s12$S2
    vec.y[ii] <- max(0,qt - k)
  }
  
  return(vec.y)
}

# --- --- --- --- --- --- --- --- 
OCrl.MACusum <- function(i, L, k, g, Vmu, Mcov, ar, Ss, hwang=F){
  # Descripción : 
  # Calculo de un RL fuera de control, muestreando de una normal multivariada 
  # dado un vector de medias fuera de control Vmu.
  # Necesita las funciones Qt y Sh.12.
  #
  # Parametros :
  # i := Es el i-esimo RL, esto es para implementar map, o barra de progreso.
  # L := es el limite de la carta
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # g := distribución del antirango
  # Vmu := vector de medias fuera de control
  # Mcov := matriz de var-cov bajo control
  # ar := Cuales son los antirangos tomados por la carta de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # Ss := Es una matriz de con  el espacio muestral del antirango [ar],
  #       necesaria para la carta hwang.
  # hwang := Si es TRUE utiliza la carta de hwang.
  #
  # Output :
  # Un entero el cual es la longitud de corrida 
  
  RL <- 0
  flag <- T
  s1 <- s2 <- numeric(length(g))
  
  while(flag){
    RL <- RL + 1
    
    if(!hwang){
      obs <- MASS::mvrnorm(1, Vmu, Mcov)
      e1t <- epsilon_h(obs, ar, Ss)
    }else{
      obs <- MASS::mvrnorm(1, Vmu, Mcov)
      e1t <- epsilon_h(obs, ar, Ss, hwang = T)
    }
    
    qt <- Qt(s1,s2,e1t,g)
    s12 <- Sh.12(s1,s2,e1t,qt,g,k)
    s1 <- s12$S1
    s2 <- s12$S2
    yt <- max(0,qt - k) # 
    flag <- yt < L
  }
  return(RL)
}

future_OCarl.MAC <- function(Vmu,ar,L,Mcov, k, n.iter=1e3, hwang = F){
  # Descripción : 
  # Calculo de un ARL fuera de control, muestreando de una normal multivariada 
  # dado un vector de medias fuera de control Vmu.
  # Necesita la funcion OCrl.MACusum
  #
  # Parametros :
  # Vmu := vector de medias fuera de control
  # ar := Cuales son los antirangos tomados por la carta de hwang,
  #       es necesario que contenga el primer o ultimo antirango.
  # L := es el limite de la carta
  # Mcov := matriz de var-cov bajo control
  # k := es el parametro de la carta \in [0,max(h){\sum{g_i}/g_h})
  # n.iter := numero de iteraciones
  # hwang := Si es TRUE utiliza la carta de hwang.
  #
  # Output :
  # Un real el cual es el ARL fuera de control

  p <- length(Vmu)
  q <- length(ar)
  Ss <- gtools::permutations(p,q,1:p)
  g1 <- numeric(perm(p,q)) + 1/perm(p,q) # se asume intercambiabilidad 
  
  plan(multiprocess)
  # plan(multiprocess, workers = 3)
  RLs <- future_map_dbl(1:n.iter, OCrl.MACusum, L, k, g1, Vmu, Mcov, ar, Ss, hwang)
  cat('Finished the ARL calculation with L = ', Vmu, ' y ', ar, '\n')
  return(mean(RLs))
}

# 9.4 ---------------------------------------------------------------------

# Base 
colnames(tbl.alum) <- paste0('X', 1:5)
nrow(tbl.alum)/2
dat.train <- tbl.alum[1:95,]
dat.test <- tbl.alum[96:nrow(tbl.alum),] %>% 
  mutate(n = 96:nrow(tbl.alum))

tbl.alum %>% 
  mutate(n = seq_along(X1)) %>% 
  gather(X, valor, -n) %>% 
  ggplot(aes(n, valor)) +
  geom_line() +
  facet_grid(X~.)

tbl.alum %>% cor %>% round(.,4)
dat.train %>% cor %>% round(.,4)

tbl.alum %>% 
  gather(X, valor) %>%
  group_by(X) %>% 
  nest() %>% 
  mutate(plots = map2(data, X, ~.x %>% 
                       ggplot(aes(valor)) +
                       geom_density(alpha = 0.2) +
                       xlab(.y) +
                       scale_x_continuous(limits = c(-3, 3))) ) %>% 
  .$plots %>% 
  autoplot()


# carta MA-CUSUM
k1 <- 1

# primer antirango --- --- ---
p <- 5
q <- 1
ar <- 1

Ss <- gtools::permutations(p,q,1:p)
epsi4 <- apply(dat.train,1,epsilon_h, ar, Ss) %>% t

g1 <- colSums(epsi4)/sum(colSums(epsi4))

# Se desea obtener un ARL0 = 200
# MA-CUSUM
rl <- pbMACusum.RL(11.85, k1, g1, 1e3)
mean(rl)

dat.train %>% 
  carta.MACusum(ar, Ss, k1, g1) %>% 
  tibble(Cn = .) %>% 
  mutate(n = seq_along(Cn)) %>% 
  ggplot(aes(n, Cn)) +
  geom_hline(yintercept = 11.85, linetype="dashed") +
  geom_point() +
  geom_line() +
  labs(title='MA-CUSUM chart based on the first antirank',
       subtitle='Train data set')

# Mirar el ordendamiento del espacio de muestreo del vector de ar
# Ss1 <- gtools::permutations(p,q,1:p)
# Ss <- c(Ss1[3:5],Ss1[1:2]) %>% cbind()
# carta.MACusum(dat.train, 1, Ss1, k1, g1)
# carta.MACusum(dat.train, 1, Ss, k1, g1) 

# En busca del L
vec.l <- 11:12
arls <- sapply(vec.l, future_arl, k1, g1, 1e3)
arls

# con los datos de prueba
dat.test %>% 
  select(-n) %>% 
  carta.MACusum(ar, Ss, k1, g1) %>% 
  tibble(Cn = .) %>% 
  mutate(n = dat.test$n) %>% 
  ggplot(aes(n, Cn)) +
  geom_hline(yintercept = 11.85, linetype="dashed") +
  geom_point() +
  geom_line() +
  labs(title='MA-CUSUM chart based on the first antirank',
       subtitle='Test data set')

# El siguiente fragmento no calcula qt en la primera iteracción
# Este fragmento del codigo es del example88.r, no se parece en nada a la grafica
# kP=1
# gnl <- epsi4
# f0 <- g1
# n <- nrow(gnl)
# Unobs=matrix(0,n,10)
# Unexp=matrix(0,n,10)
# CnP = rep(0,n)
# Bn = rep(0,n)
# 
# if(Bn[1] > kP){
#   Unobs[1,]=gnl[1,]*(1-kP/Bn[1])
#   Unexp[1,]=f0*(1-kP/Bn[1])
# }
# if(Bn[1] <= kP){
#   Unobs[1,]=t(rep(0,10))
#   Unexp[1,]=t(rep(0,10))
# }
# 
# CnP[1]=sum((Unobs[1,]-Unexp[1,])^2/Unexp[1,])
# 
# for(i in 2:n){
# 
#   Bn[i]=sum(((Unobs[i-1,]-Unexp[i-1,])+(gnl[i,]-f0))^2/(Unexp[i-1,]+f0))
# 
#   if(Bn[i]>kP){
#     Unobs[i,]=(Unobs[i-1,]+gnl[i,])*(1-kP/Bn[i])
#     Unexp[i,]=(Unexp[i-1,]+f0)*(1-kP/Bn[i])
#   }
#   if(Bn[i]<=kP){
#     Unobs[i,]=t(rep(0,10))
#     Unexp[i,]=t(rep(0,10))
#   }
# 
#   CnP[i]=sum((Unobs[i,]-Unexp[i,])^2/Unexp[i,])
# 
# }


# primer y ultimo antirango --- --- ---
p <- 5
q <- 2
Ss <- gtools::permutations(p,q,1:p)
ar <- c(1,5)
epsi4 <- apply(dat.train,1,epsilon_h, ar, Ss) %>% t

conteo <- colSums(epsi4)
conteo[conteo == 0] <- 0.5
g1 <- conteo/sum(conteo)

# Se desea obtener un ARL0 = 200
# MA-CUSUM
rl <- pbMACusum.RL(34.89, k1, g1)
mean(rl)

dat.train %>% 
  carta.MACusum(ar, Ss, k1, g1) %>% 
  tibble(Cn = .) %>% 
  mutate(n = seq_along(Cn)) %>% 
  ggplot(aes(n, Cn)) +
  geom_hline(yintercept = 34.89, linetype="dashed") +
  geom_point() +
  geom_line() +
  labs(title='MA-CUSUM chart based on the first antirank',
       subtitle='Train data set')

# En busca del L
vec.l <- 34:35
tic()
arls <- sapply(vec.l, future_arl, k1, g1, 1e3)
arls
toc()

# con los datos de prueba
dat.test %>% 
  select(-n) %>% 
  carta.MACusum(ar, Ss, k1, g1) %>% 
  tibble(Cn = .) %>% 
  mutate(n = dat.test$n) %>% 
  ggplot(aes(n, Cn)) +
  geom_hline(yintercept = 34.89, linetype="dashed") +
  geom_point() +
  geom_line() +
  labs(title='MA-CUSUM chart based on the first antirank',
       subtitle='Test data set')


# 9.5 ---------------------------------------------------------------------

k1 <- 1
# ar1-ar4 --- --- ---
p <- 4
ar <- 1
q <- length(ar)

Ss <- gtools::permutations(p,q,1:p)

g1 <- numeric(p) + 1/p

# Se desea obtener un ARL0 = 200
# MA-CUSUM
rl <- pbMACusum.RL(6.842, k1, g1, 1e4)
mean(rl)

# ar12-ar34 --- --- ---
p <- 4
ar <- 1:2
q <- length(ar)

Ss <- gtools::permutations(p,q,1:p)

g1 <- numeric(perm(p,q)) + 1/perm(p,q)

# Se desea obtener un ARL0 = 200
# MA-CUSUM
rl <- pbMACusum.RL(15.6887, k1, g1, 1e4)
mean(rl)

Para la 

# Simulaciones --- --- --- --- ---
k1 <- 1
p <- 4
l <- 15.6887
mu1 <- c(-4,0,0,0)
cov0 <- diag(p)

# ARL OC bajo la carta MA-CUSUM para diferentes cambios de vector de medias
# y distintos ar

# ARs <- list(1,2,3,4,1:2,c(1,3),c(1,4),2:3,c(2,4),3:4)
# ls <- c(rep(6.842, 4 ), rep(15.6887, 6 ))
# mu.list <- list(c(-4,0,0,0), c(-4,-2,0,0), c(-4,-4,0,0), c(-4,-4,-2,0),
#                 c(-4,-2,-2,0), c(-4,-3,-1,0), c(-4,-4,-4,0))
# 
# tbl92 <- tibble(ar = rep(ARs,length(mu.list)),
#                 l = rep(ls, length(mu.list)),
#                 Vmu = rep(mu.list, each = length(ARs))) %>% 
#   # filter(seq_along(ar) %in% c(1)) %>% 
#   mutate( arl = pmap_dbl(list(Vmu,ar,l), future_OCarl.MAC, cov0 ,k1, 1e4))
# 
# save(tbl92, file = file.path(outPath, 'tbl92.Rdata'))

load(file.path(outPath, 'tbl92.Rdata'))

pr_tbl92 <- tbl92 %>% 
  select(-l) %>% 
  mutate(Vmu = map_chr(Vmu, ~paste0('(',paste0(.x, collapse = ','),')')) ,
         ar = map_chr(ar,~paste0('AR',paste0(.x, collapse = ''))) ) %>% 
  spread(Vmu, arl)

