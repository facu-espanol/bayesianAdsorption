#calcula prediccion simulada con error normal
PredNorm <- function(tt,Mod,par){
  val <- rnorm(1,Mod$Cr(tt,par[1:Mod$npar]),par[(Mod$npar+1)]) ##par[(Mod$npar+1)]
  return(val)
}
#modelo de error
Ermod.Norm <- list(llikd=LogLikNorm.dens,pred=PredNorm,npar=1,parint=c(0,0.1),parini=0.05)

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
#' Ejecucions
#'
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




#' Title
#'
#' @param Mod Modelo
#' @param Dat Datos
#' @param config Objeto configuracion
#'
#' @returns curva
#' @export
#'
#' @examples ejemplo1
PostSamp <- function(
    Mod,
    Dat,
    config = list()
) {

  cfg <- modifyList(PostSampDefaults, config)

  if (is.null(cfg$modpar.int)) {
    cfg$modpar.int <- Mod$parint
  }

  if (is.null(cfg$ermodpar.int)) {
    cfg$ermodpar.int <- cfg$Ermod$parint
  }

  if (is.null(cfg$ermodpar.ini)) {
    cfg$ermodpar.ini <- cfg$Ermod$parini
  }


  # PRIOR MULTIVARIADA

  if (!is.null(cfg$modpar.multivariada)) {

    # La multivariada contiene el vector de medias
    # y la matriz de covarianza.
    mu <- cfg$modpar.multivariada$media
    cfg$modpar.covarianza <- cfg$modpar.multivariada$covarianza

    # Se obtienen los desvíos individuales de la diagonal.
    # Se mantienen por compatibilidad con la likelihood.
    sigma <- sqrt(diag(cfg$modpar.covarianza))

    # Si no se proporciona un valor inicial,
    # el algoritmo comienza desde la media de la multivariada.
    if (is.null(cfg$modpar.ini)) {
      cfg$modpar.ini <- mu
    }

  } else {

    # PRIOR INDEPENDIENTE O SIN PRIOR

    if (is.null(cfg$modpar.ini)) {
      cfg$modpar.ini <- Mod$parini
    }

    if (
      is.null(cfg$modpar.ini) &&
      !is.null(cfg$modpar.prior)
    ) {
      cfg$modpar.ini <- sapply(
        cfg$modpar.prior,
        function(x) x$mean
      )
    }

    if (!is.null(cfg$modpar.prior)) {

      mu <- sapply(
        cfg$modpar.prior,
        function(x) x$mean
      )

      sigma <- sapply(
        cfg$modpar.prior,
        function(x) x$sd
      )

    } else {

      mu <- NULL
      sigma <- NULL
      cfg$modpar.covarianza <- NULL
    }
  }

  # validar que modpar.multivariada tenga
  # los elementos "media" y "covarianza".
  #
  # validar las dimensiones de media y covarianza.
  #
  # validar que la covarianza sea simétrica.
  #
  # validar que la covarianza sea definida positiva.
  #
  # validar que modpar.ini esté dentro
  # de los intervalos permitidos.

  cfg$mu    <- mu
  cfg$sigma <- sigma

  ret <- mcmc(Mod, Dat, cfg)
  return(ret)
}
