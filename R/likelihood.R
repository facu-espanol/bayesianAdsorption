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


# FALTA AGREGAR LA LOG VEROSIMILITUD PERO DIFERENCIADA DE LA OTRA. NO SON LO MISMO
