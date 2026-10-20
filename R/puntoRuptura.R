#Recibe
#ps: ejecucion de postsamp
#ttseq: secuencia de tiempo
#objetivo: punto a identificar
#tolerancia: diferencia permitida entre objetivo y valor real
#' Title
#'
#' @param ps ejecucion
#' @param ttseq tiempo
#' @param objetivo que punto
#' @param tolerancia +- cuanto
#'
#' @returns distribucion de probabilidad de objetivo
#' @export
#'
#' @examples ejemplo1
punto_ruptura <- function(
    ps,
    ttseq = NULL,
    objetivo = 0.1,
    tolerancia = 0.005
) {

  Mod <- ps$mod

  #toma la secuencia de tiempo de la ejecucion
  if (is.null(ttseq)) {
    ttseq <- seq(
      min(ps$data$tt, na.rm = TRUE),
      max(ps$data$tt, na.rm = TRUE),
      by = 1
    )
  }

  # Curvas para todas las muestras posteriores. Para cada combinacion de a y b hace la curva
  pscr <- matrix(
    NA,
    nrow = ps$nm.ps,
    ncol = length(ttseq)
  )

  #para cada conjunto de parametros, calcula la curva.
  for (ii in seq_len(ps$nm.ps)) {
    pscr[ii, ] <- Mod$Cr(
      ttseq,
      ps$param.ps[ii, 1:Mod$npar]
    )
  }

  # Curva mediana de la posterior. cada tiempo queda asociado con una concentracion
  curva <- apply(
    pscr,
    2,
    median,
    na.rm = TRUE
  )

  # Buscar dónde Cr está cerca del objetivo
  idx <- which(
    abs(curva - objetivo) <= tolerancia
  )

  if (length(idx) == 0) {
    return(NA_real_)
  }

  #devuelve el primer punto que llego a 0.1+-0.005
  return(ttseq[idx[1]])
}
