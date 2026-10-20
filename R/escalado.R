#entradas
#Salida de postsamp (estimados)
#op.origen y destino pueden contener: vd, c0, A, u, M, h
#funcion de relacion entre parametros operativos

#salida
#la misma que si ejecuto el modelo normal (postsamp). deberia ser como ejecutar el modelo

#' Title
#'
#' @param ps corrida
#' @param op.origen operativas origen
#' @param op.objetivo operativas destino
#' @param funcion.escalado funcion de escalado
#'
#' @returns muestreo escalado
#' @export
#'
#' @examples ejemplo1
EscalarPostSamp <- function(
    ps,
    op.origen,
    op.objetivo,
    funcion.escalado
) {
  Mod <- ps$mod

  #no es necesario, los parametros se pasan por valor. es solo por legibilidad
  resultado <- ps

  matriz_muestras <- ps$param.ps

  #recorre todas las muestras
  for (ii in seq_len(ps$nm.ps)) {

    #extrae parametros del modelo. el npar hace que sean todos menos el de error.
    parametros.modelo <- ps$param.ps[
      ii,
      seq_len(Mod$npar)
    ]

    #se llama a la funcion de escalado para esa fila de la matriz
    parametros.escalados <- funcion.escalado(
      parametros.modelo,
      op.origen,
      op.objetivo
    )

    #construye la matriz de muestras escaladas
    matriz_muestras[
      ii,
      seq_len(Mod$npar)
    ] <- parametros.escalados
  }

  resultado$param.ps <- matriz_muestras

  modpar.mgt <- resultado$param.ps[
    ,
    seq_len(Mod$npar)
  ]

  #media de cada parametro
  resultado$posterior.multivariada <- list(
    media = colMeans(
      modpar.mgt,
      na.rm = TRUE
    ),
    #matriz de covarianza
    covarianza = cov(
      modpar.mgt,
      use = "complete.obs"
    ),
    #matriz de correlacion
    correlacion = cor(
      modpar.mgt,
      use = "complete.obs"
    ),
    muestras = modpar.mgt
  )

  return(resultado)
}


#' EscalarMdr
#'
#' @param parametros par de parametros a escalar
#' @param op.origen operativas origen
#' @param op.objetivo operativas objetivo
#'
#' @returns par (a, b) escalados
#' @export
#'
#' @examples ejemplo1
EscalarMdr <- function(
    parametros,
    op.origen,
    op.objetivo
) {
  # Parámetros estimados en la corrida anterior
  a.origen <- parametros[1]
  b.origen <- parametros[2]

  # Recupero parámetros físicos desde la corrida anterior
  kyv <- b.origen
  q0 <- (op.origen$c0 * op.origen$A * op.origen$u)/(exp(a.origen/kyv)*op.origen$M)

  # Calculo los nuevos parámetros a y b para la condición objetivo
  a.nuevo <- kyv * log((op.objetivo$c0 * op.objetivo$A * op.objetivo$u)/(q0 * op.objetivo$M))
  b.nuevo <- kyv

  return(c(a.nuevo, b.nuevo))
}

#' EscalarLogistico
#'
#' @param parametros par de parametros a escalar
#' @param op.origen operativas origen
#' @param op.objetivo operativas objetivo
#'
#' @returns par (a, b) escalados
#' @export
#'
#' @examples ejemplo1
EscalarLogistico <- function(
    parametros,
    op.origen,
    op.objetivo
) {
  # Parámetros estimados en la corrida anterior
  a.origen <- parametros[1]
  b.origen <- parametros[2]

  # Recupero parámetros físicos desde la corrida anterior
  kba <- b.origen / op.origen$c0
  n0 <- (a.origen * op.origen$u) /(op.origen$h * kba)

  # Calculo los nuevos parámetros a y b para la condición objetivo
  a.nuevo <- kba * n0 * (op.objetivo$h / op.objetivo$u)
  b.nuevo <- kba * op.objetivo$c0

  return(c(a.nuevo, b.nuevo))
}

#' EscalarNormal
#'
#' @param parametros par de parametros a escalar
#' @param op.origen operativas origen
#' @param op.objetivo operativas objetivo
#'
#' @returns par (a, b) escalados
#' @export
#'
#' @examples ejemplo1
EscalarNormal <- function(
    parametros,
    op.origen,
    op.objetivo
) {
  # Parámetros estimados en la corrida anterior
  a.origen <- parametros[1]
  b.origen <- parametros[2]

  # Recupero parámetros físicos desde la corrida anterior
  sigma.a <- op.origen$h / a.origen
  vd.origen <- b.origen * sigma.a

  # Calculo los nuevos parámetros a y b para la condición objetivo
  a.nuevo <- op.objetivo$h / sigma.a
  b.nuevo <- vd.origen / sigma.a

  return(c(a.nuevo, b.nuevo))
}
