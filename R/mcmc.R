#' Title
#'
#' @param Mod
#' Modelo a utilizar
#' @param Dat
#' Datos experimentales
#' @param cfg
#' Configuracion
#'
#' @returns
#' Ejecucion
#'
#' @export
#'
#' @examples
#' ejemplo
mcmc <- function(Mod, Dat, cfg) {
  # Desempaqueta cfg
  Ermod               <- cfg$Ermod
  modpar.int          <- cfg$modpar.int
  modpar.ini          <- cfg$modpar.ini
  modpar.multivariada <- cfg$modpar.multivariada
  modpar.covarianza   <- cfg$modpar.covarianza
  modpar.prior        <- cfg$modpar.prior
  ermodpar.int        <- cfg$ermodpar.int
  ermodpar.ini        <- cfg$ermodpar.ini
  n.gibbs             <- cfg$n.gibbs
  n.burnin            <- cfg$n.burnin
  n.thin              <- cfg$n.thin
  mu                  <- cfg$mu
  sigma               <- cfg$sigma

  #nc <- dim(Dat)[1] #no se usa nunca. numero observaciones.
  npar <- Mod$npar+Ermod$npar
  ###Construye las variables de muestreo
  pram.gibbs <- matrix(NA, nrow=n.gibbs, ncol=npar)
  stht.gibbs <- matrix(NA, nrow=n.burnin, ncol=npar) ## step size en metropolis de cada parámetro
  stac.gibbs <- matrix(0, nrow=n.burnin, ncol=npar)  ## guarda la cantidad de aceptados en metropolis
  ###Condiciones iniciales
  validpar <- cbind(modpar.int,ermodpar.int)  ## intervalo de validez para todos los parámetros
  colnames(validpar) <- NULL
  pram.g <- c(modpar.ini,ermodpar.ini)  ## parámetros iniciales
  del <- 1.01               ## delta en el ajuste del paso
  stht.g <- rep(0.01,npar)   ## step size iniciales
  ###Muestreo
  ##Etapa de construcción
  for(g in 1:(n.burnin)){
    ##Guarda los valores actuales
    pram.gibbs[g,] <- pram.g
    stht.gibbs[g,] <- stht.g
    ##Muestrea los valores siguientes
    z <- runif(npar,0,1) # números aleatorios para usar en los metropolis
    for(pp in 1:npar){
      test <- rnorm(1,pram.g[pp],stht.g[pp])
      stht.g[pp] <- del^(-1)*stht.g[pp]
      if(validpar[1,pp] <= test & test < validpar[2,pp]){
        pramtest <- replace(pram.g,pp,test)
        deltest <- exp(Ermod$llikd(Mod,Dat,pramtest, mu, sigma, modpar.covarianza)-Ermod$llikd(Mod,Dat,pram.g, mu, sigma, modpar.covarianza))
        if(z[pp] < deltest){
          pram.g[pp] <- test
          stht.g[pp] <- del^4*stht.g[pp]
          stac.gibbs[g,pp] <- 1
        }
      }
    }
  }

  ##Fija los pazos para el muestreo siguiente
  stht.means <- colMeans(stht.gibbs[(n.burnin*0.9):n.burnin,])
  accept.rate <- colSums(stac.gibbs[(n.burnin*0.9):n.burnin,])/(n.burnin-n.burnin*0.9+1)
  ##Etapa de muestreo propiamente dicho
  for(g in (n.burnin+1):(n.gibbs-1)){
    ##Guarda los valores actuales
    pram.gibbs[g,] <- pram.g
    ##Muestrea los valores siguientes
    z <- runif(npar,0,1) # números aleatoreos para usar en los metropolis
    for(pp in 1:npar){
      test <- rnorm(1,pram.g[pp],stht.means[pp])
      if(validpar[1,pp] <= test & test < validpar[2,pp]){
        pramtest <- replace(pram.g,pp,test)
        deltest <- exp(Ermod$llikd(Mod,Dat,pramtest, mu, sigma, modpar.covarianza)-Ermod$llikd(Mod,Dat,pram.g, mu, sigma, modpar.covarianza))
        if(z[pp] < deltest){
          pram.g[pp] <- test
        }
      }
    }
  }
  ##Guarda el último valor
  pram.gibbs[n.gibbs,] <- pram.g
  ###Muestreo intercalado
  seq.mgt=seq(n.burnin+1,n.gibbs,by=n.thin)
  n.mgt<- length(seq.mgt)
  pram.mgt <- pram.gibbs[seq.mgt,]
  ###Mediana e intervalo del 90% de credibilidad
  #ver si poner:
  #pram.mci.95 <- apply(pram.mgt,2,quantile,probs=c(0.025,0.5,0.975),na.rm=TRUE)


  # ==================================================
  # POSTERIOR MULTIVARIADA
  # ==================================================

  # Seleccionar solo los parámetros del modelo.
  # Los parámetros del error quedan afuera.
  modpar.mgt <- pram.mgt[
    ,
    1:Mod$npar
  ]

  posterior.multivariada <- list(
    # Vector de medias posterior
    media = colMeans(
      modpar.mgt,
      na.rm = TRUE
    ),

    # Matriz de covarianza posterior
    covarianza = cov(
      modpar.mgt,
      use = "complete.obs"
    ),

    # Matriz de correlación posterior
    correlacion = cor(
      modpar.mgt,
      use = "complete.obs"
    ),

    # Muestras conjuntas de la posterior
    muestras = modpar.mgt
  )

  ### Resultado
  ret <- list(
    param.ps = pram.mgt,
    nm.ps = n.mgt,
    stht.means = stht.means,
    accept.rate = accept.rate,
    posterior.multivariada = posterior.multivariada,
    data = Dat,
    mod = Mod,
    ermod = Ermod
  )
  return(ret)
}
