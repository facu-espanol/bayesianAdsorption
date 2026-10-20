MinCuadrados <- function(parametros, Dat, Mod) {

  prediccion <- Mod$Cr(Dat$tt, parametros)

  #suma de (datos reales - prediccion)^2
  return (sum((Dat$Cr - prediccion)^2))
}

#' Title
#'
#' @param Dat datos
#' @param Mod modelo
#' @param parametros_iniciales parametros iniciales
#'
#' @returns par (a, b) optimizado
#' @export
#'
#' @examples ejemplo1
AjustarMinCuadrados <- function(Dat, Mod, parametros_iniciales) {
  #optim busca los valores optimpos para los parametros a y b, usando la funcion minimos cuadrados
  resultado <- optim(
    par = parametros_iniciales,
    fn = MinCuadrados,
    Dat = Dat,
    Mod = Mod)
}
