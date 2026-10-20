#' Title
#'
#' @param ps ejecucion
#'
#' @returns cuantiles
#' @export
#'
#' @examples ejemplo1
PostSummary <- function(ps){
  rbind(
    media   = colMeans(ps$param.ps, na.rm = TRUE),
    sd      = apply(ps$param.ps, 2, sd, na.rm = TRUE),
    q025    = apply(ps$param.ps, 2, quantile, probs = 0.025, na.rm = TRUE),
    mediana = apply(ps$param.ps, 2, quantile, probs = 0.5, na.rm = TRUE),
    q975    = apply(ps$param.ps, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}
