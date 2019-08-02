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
rl.MACusum <- function(i, L, k, g, hwang=F, ar = NULL, Ss=NULL, tbl.distri = NULL){
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
      if(is.null(ar) | is.null(Ss)) stop('Arguments ar o Ss do not define')
      if(is.null(tbl.distri)){
        if(!(exists('Vmu') & exists('Mcov'))) stop('Variables Vmu & Mcov need to be created in the GlobalEnv')
        obs <- MASS::mvrnorm(1, Vmu, Mcov)
        e1t <- epsilon_h(obs, ar, Ss, hwang = T)
      }else{
        e1t <- set.epsi2[sample(1:nrow(set.epsi2),1),]
      }
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



future_arl <- function(L, k, g, n.iter=1e3, hwang = F, ar = NULL, Ss=NULL, tbl.distri = NULL){
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
  RLs <- future_map_dbl(1:n.iter, rl.MACusum,L, k, g, hwang,ar,Ss,tbl.distri)
  # cat('Finished the ARL calculation with L = ', L, '\n')
  #list_return <- list(ARL = mean(RLs), sdRL = sd(RLs))
  #return(list_return)
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
  }else{ # revisar cuando tiene AR1p
    if(1 %in% ar){
      return(vep_op*abs(min(obs)))
    }else if(length(obs) %in% ar){
      return(vep_op*abs(max(obs)))
    }else{
      stop('For hwang chart is necessary that arg(ar) must been defined with the first or last antirank')
    }
  }
}


# Calcula la distribución de g apartir de unos datos
distr.g <- function(dataSet, ar, Ss, hwang = F){
  set.epsi <- apply(dataSet, 1, epsilon_h, ar, Ss, hwang) %>% t
  set.indi <- set.epsi
  set.indi[set.indi != 0] <- 1 # En caso que sea hwang es para dejarlo con unos y ceros
  g <- colSums(set.indi)/sum(colSums(set.indi)) # Solo sumar ceros y unos
  if( !all(g != 0) ){
    # Realizo esta correción de ceros en la distribución del vector [ar]
    # porque considero que es poco probable que sea igual a cero 
    # sino que es poco probable, con una probabilidad pequeña
    # Qué valor dejarle en vez de 0.5?
    warning(paste0('There was a correction in vector g of ', sum(g == 0), 
                   ' elements in 0 by 0.5'))
    conteo <- colSums(set.indi)
    conteo[conteo == 0] <- 0.5
    g <- conteo/sum(conteo)
  }
  return(g)
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
    g <- distr.g(dataSet, ar, Ss, hwang)
    cat('La distribución g = ', paste0(round(g,3),collapse = ', '),'\n')
  }
  s1 <- s2 <- numeric(length(g))
  
  vec.y <- numeric(nrow(dataSet))
  for (ii in 1:nrow(dataSet)) {
    e1t <- set.epsi[ii,]
    qt <- Qt(s1,s2,e1t,g)
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
  # Retorna una lista que contiene el ARL y su desviación estandar
  
  p <- length(Vmu)
  q <- length(ar)
  Ss <- gtools::permutations(p,q,1:p)
  g1 <- numeric(perm(p,q)) + 1/perm(p,q) # se asume intercambiabilidad 
  
  plan(multiprocess)
  # plan(multiprocess, workers = 3)
  RLs <- future_map_dbl(1:n.iter, OCrl.MACusum, L, k, g1, Vmu, Mcov, ar, Ss, hwang)
  cat('Finished the ARL calculation with L = ', Vmu, ' y ', ar, '\n')
  list_return <- list(ARL = mean(RLs), sdRL = sd(RLs))
  return(list_return)
}