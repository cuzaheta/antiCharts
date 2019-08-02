# functions ---------------------------------------------------------------
# mvrnorm <- MASS::mvrnorm
# perm <- function(n,k){choose(n,k) * factorial(k)}

Qt <- function(S1_1,S2_1,et,g1=g){
  A <- diag((S2_1+g1)^(-1))
  a <- (S1_1-S2_1)+(et-g1)
  return(t(a)%*%A%*%a %>% as.numeric())
}

# Yt <- function(S1,S2){
#   A <- diag(S2^(-1))
#   a <- S1-S2
#   return(t(a)%*%A%*%a %>% as.numeric())
# }

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

RL2.pb <- function(i, alter =F){
  setTxtProgressBar(pb, i)
  RL <- 0
  flag <- T
  s1 <- S1.0
  s2 <- S2.0
  
  while(flag){
    RL <- RL + 1
    
    if(!alter){
      e1t <- rmultinom(1, 1, g)[,1]
      # obs <- mvrnorm(1, Vmu, Mcov)
      # e1t <- epsilon_h(obs, vec, Ss2.4)
    }else{
      obs <- mvrnorm(1, Vmu, Mcov)
      e1t <- epsilon_hC(obs, vec, Ss2.4)
    }
    qt <- Qt(s1,s2,e1t)
    # S1_1,S2_1,et,Qt,g1=g,k1=k
    s12 <- Sh.12(s1,s2,e1t,qt)
    s1 <- s12$S1
    s2 <- s12$S2
    yt <- max(0,qt - k) # 
    flag <- yt < L
  }
  return(RL)
}

# Funciones para generar el vector epsilon para la carta 
# este vector esta en base al conjunto S 
epsil_h <- function(uno, Ss) apply(Ss, 1, function(x) ifelse(identical(x,uno),1,0))
epsil_h <- function(uno, Ss) apply(Ss, 1, function(x) ifelse(all(x == uno),1,0))

# se puede borrar la función epsil_h, meterla en las dos de abajo
epsilon_h <- function(obs, vec_ep, Ss=NULL){
  # Original
  # saca el vector epsilon de S, en base a la observación de una muestra
  ep <- order(obs)[vec_ep]
  return(epsil_h(ep, Ss))
}

epsilon_hC <- function(obs, vec_ep, Ss=NULL){
  # saca el vector epsilon de S, en base a la propuesta de hwang
  # ep <- order(obs)[vec_ep]
  # vep_op <- epsil_h(ep, Ss)
  vep_op <- epsilon_h(obs,vec_ep, Ss)
  
  if(1 %in% vec_ep){
    return(vep_op*abs(min(obs)))
  }else if(length(obs) %in% vec_ep){
    return(vep_op*abs(max(obs)))
  }else{
    message('Problemas de dimensión entre la observacion y epsilon')
  }
}

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
























