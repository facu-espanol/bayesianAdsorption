#Parametro de entrada: salida de PostSamp
#Salida graficos
#' Title
#'
#' @param ps ejecucion
#'
#' @returns histograma
#' @export
#'
#' @examples ejemplo1
PlotTrace <- function(ps){
  npar <- ncol(ps$param.ps)

  for(i in 1:npar){
    plot(
      ps$param.ps[,i],
      type = "l",
      main = paste("Traceplot parametro", i),
      xlab = "Iteracion",
      ylab = paste("Parametro", i)
    )
  }
}



#Entradas
#ps : salida de postsamp
#mci: matriz de cuantiles
#i: numero del parametro
#' Title
#'
#' @param ps ejecucion
#' @param mci parametros
#' @param i numero de parametro
#'
#' @returns distribucion de un parametro
#' @export
#'
#' @examples ejemplo1
GraficarParametro <- function(ps, mci, i) {
  plot(
    density(ps$param.ps[, i]),
    main = paste("Densidad del parámetro", i),
    xlab = paste("Parámetro", i)
  )

  abline(v = mci[1, i], lty = 2)
  abline(v = mci[2, i], lty = 1)
  abline(v = mci[3, i], lty = 2)
}


#Entradas
#ps: ejecucion del modelo
#Dat: si se quieren graficar los puntos de un set de datos
#ttseq: tiempo, determina la escala del eje x
#probs: que lineas queremos marcar
#titulo: titulo del grafico
#xlab: nombre del eje x
#ylab: nombre del eje y
#' Title
#'
#' @param ps Ejecucion
#' @param Dat datos
#' @param ttseq secuencia de tiempo
#' @param probs cuantiles
#' @param titulo titulo del grafico
#' @param xlab eje x
#' @param ylab eje y
#'
#' @returns grafico de curva de ruptura
#' @export
#'
#' @examples ejemplo1
GraficarAjusteModelo <- function(
    ps,
    Dat=NULL,
    ttseq = NULL,
    probs = c(0.025, 0.5, 0.975),
    titulo = "Ajuste del modelo",
    xlab = "Tiempo",
    ylab = "Cr"
) {
  if (is.null(Dat)){
    Dat <- ps$data
  }
  Mod <- ps$mod
  if (is.null(ttseq)){
    ttseq <- seq(
      min(Dat$tt, na.rm = TRUE),
      max(Dat$tt, na.rm = TRUE),
      by = 1
    )
  }
  pscr <- matrix(
    NA,
    nrow = ps$nm.ps,
    ncol = length(ttseq)
  )

  for (ii in seq_len(ps$nm.ps)) {
    pscr[ii, ] <- Mod$Cr(
      ttseq,
      ps$param.ps[ii, 1:Mod$npar]
    )
  }
  #cuantiles de las curvas
  pscr.mci <- apply(
    pscr,
    2,
    quantile,
    probs = probs,
    na.rm = TRUE
  )

  #cuantiles de los parametros
  pspar.mci <- apply(
    ps$param.ps,
    2,
    quantile,
    probs = probs,
    na.rm = TRUE
  )

  #10% saturacion
  t.ruptura <- punto_ruptura(
    ps = ps,
    ttseq = ttseq
  )

  plot(
    Dat$tt,
    Dat$Cr,
    xlim = range(c(Dat$tt, ttseq), na.rm = TRUE),
    ylim = range(c(0, 1)),
    xlab = xlab,
    ylab = ylab,
    main = titulo
  )

  #curva central
  lines(
    ttseq,
    Mod$Cr(
      ttseq,
      pspar.mci[2, 1:Mod$npar]
    )
  )

  lines(ttseq, pscr.mci[1, ], lty = 2)
  lines(ttseq, pscr.mci[2, ], lty = 2)
  lines(ttseq, pscr.mci[3, ], lty = 2)
  # Punto de ruptura
  abline(
    h = 0.1,
    lty = 3
  )

  abline(
    v = t.ruptura,
    lty = 3
  )

}

#' Title
#'
#' @param ps ejecucion
#' @param Dat datos
#' @param ttseq secuencia de tiempo
#' @param probs cuantiles
#' @param saturacion hasta donde quiero llegar con la grafica
#' @param titulo titulo del grafico
#' @param xlab eje x
#' @param ylab eje y
#'
#' @returns grafico curva de ruptura
#' @export
#'
#' @examples ejemplo1
GraficarEscalado <- function(
    ps,
    Dat = NULL,
    ttseq = NULL,
    probs = c(0.025, 0.5, 0.975),
    saturacion = 0.99,
    titulo = "Modelo escalado",
    xlab = "Tiempo",
    ylab = "Cr"
) {

  Mod <- ps$mod


  # ==================================================
  # PARAMETROS CENTRALES
  # ==================================================

  parametros.mediana <- apply(
    ps$param.ps[, 1:Mod$npar, drop = FALSE],
    2,
    median,
    na.rm = TRUE
  )


  # ==================================================
  # DETERMINAR HASTA QUE TIEMPO GRAFICAR
  # ==================================================
  if (is.null(ttseq)) {

    t.fin <- 0
    paso <- 10

    while (
      Mod$Cr(t.fin, parametros.mediana) < saturacion
    ) {
      t.fin <- t.fin + paso
    }

    ttseq <- seq(
      0,
      t.fin,
      by = 1
    )
  }


  # ==================================================
  # CURVAS POSTERIORES
  # ==================================================

  pscr <- matrix(
    NA,
    nrow = ps$nm.ps,
    ncol = length(ttseq)
  )

  for (ii in seq_len(ps$nm.ps)) {

    pscr[ii, ] <- Mod$Cr(
      ttseq,
      ps$param.ps[ii, 1:Mod$npar]
    )
  }


  # ==================================================
  # INTERVALO DE CREDIBILIDAD
  # ==================================================

  pscr.mci <- apply(
    pscr,
    2,
    quantile,
    probs = probs,
    na.rm = TRUE
  )


  # ==================================================
  # GRAFICO
  # ==================================================

  if (is.null(Dat)) {

    plot(
      ttseq,
      pscr.mci[2, ],
      type = "n",
      xlab = xlab,
      ylab = ylab,
      main = titulo,
      ylim = range(pscr.mci, na.rm = TRUE)
    )

  } else {

    plot(
      Dat$tt,
      Dat$Cr,
      xlim = range(c(Dat$tt, ttseq), na.rm = TRUE),
      ylim = range(c(Dat$Cr, pscr.mci), na.rm = TRUE),
      xlab = xlab,
      ylab = ylab,
      main = titulo
    )
  }


  # Curva central
  #lines(
  # ttseq,
  #Mod$Cr(
  # ttseq,
  #parametros.mediana
  #)
  #)
  lines(
    ttseq,
    pscr.mci[2, ]
  )


  # Intervalo de credibilidad
  lines(
    ttseq,
    pscr.mci[1, ],
    lty = 2
  )

  lines(
    ttseq,
    pscr.mci[2, ],
    lty = 2
  )

  lines(
    ttseq,
    pscr.mci[3, ],
    lty = 2
  )


  # ==================================================
  # PUNTO DE RUPTURA
  # ==================================================

  t.ruptura <- punto_ruptura(
    ps = ps,
    ttseq = ttseq
  )

  abline(
    h = 0.1,
    lty = 3
  )

  abline(
    v = t.ruptura,
    lty = 3
  )
}
