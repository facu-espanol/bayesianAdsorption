# Modelo de error normal: calcula la log-verosimilitud y genera predicciones
# agregando ruido gaussiano (σ constante) a la salida del modelo determinístico.
#' Title
#'
#' @param Mod Modelo
#' @param Dat Datos de entrada
#' @param par Par (a,b)
#' @param mu opcional, si hay distribucion a priori
#' @param sigma opcional, si hay distribucion a priori
#' @param covarianza opcional, si hay distribucion a priori
#'
#' @returns log verosimilitud modificada
#' @export
#'
#' @examples ejemplo1
LogLikNorm.dens <- function(Mod, Dat, par, mu = NULL, sigma = NULL, covarianza=NULL){

  Cf <- Mod$Cr(Dat$tt, par[1:Mod$npar])

  #Calcula la log-verosimilitud
  val <- sum(dnorm(Dat$Cr, Cf, par[(Mod$npar+1)], log = TRUE),na.rm = TRUE);

  #si se entro a mcmc con distribucion a priori, entonces se suma el ln de la probabilidad multivariada
  if(!is.null(mu) && !is.null(sigma)){
    theta <- par[1:Mod$npar]  #primer parametro, los parametros a y b, no el error
    if (is.null(covarianza)){
      covarianza <- diag(sigma^2) #matriz de covarianza normal multivariada
    }

    val <- val + mvtnorm::dmvnorm(
      theta,
      mean = mu,
      sigma = covarianza,
      log = TRUE
    ) #como era log sumo en vez de multiplicar log(a*b)=log(a)+log(b)
    #ln(sumatoria*probabilidad)=ln(sumatoria)+ln(probabilidad)
  }

  return(val)
}


# Modelo de error normal: calcula la log-verosimilitud y genera predicciones
# agregando ruido gaussiano (σ constante) a la salida del modelo determinístico.
#' Title
#'
#' @param Mod Modelo
#' @param Dat Datos
#' @param par Par
#' @param mu opcional, si hay distribucion a priori
#' @param sigma opcional, si hay distribucion a priori
#' @param covarianza opcional, si hay distribucion a priori
#'
#'
#' @returns log verosimilitud
#' @export
#'
#' @examples ejemplo1
LogLikNorm.dens <- function(Mod, Dat, par, mu = NULL, sigma = NULL, covarianza=NULL){

  Cf <- Mod$Cr(Dat$tt, par[1:Mod$npar])

  #Calcula la log-verosimilitud
  val <- sum(dnorm(Dat$Cr, Cf, par[(Mod$npar+1)], log = TRUE),na.rm = TRUE);

  return(val)
}

#calcula prediccion simulada con error normal
PredNorm <- function(tt,Mod,par){
  val <- rnorm(1,Mod$Cr(tt,par[1:Mod$npar]),par[(Mod$npar+1)]) ##par[(Mod$npar+1)]
  return(val)
}
#modelo de error
Ermod.Norm <- list(llikd=LogLikNorm.dens,pred=PredNorm,npar=1,parint=c(0,0.1),parini=0.05)
