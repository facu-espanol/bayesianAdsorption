#Model.ps: salida de PostSamp
#' Title
#'
#' @param Model.ps ejecucion
#'
#' @returns dic
#' @export
#'
#' @examples ejemplo1
dic.modelps <- function(Model.ps){
  ## desviación relativa media E(D[par])+C
  devc.mgt <- rep(NA,Model.ps$nm)
  for(ii in 1:Model.ps$nm){
    devc.mgt[ii] <- (-2)*Model.ps$ermod$llikd(Model.ps$mod,Model.ps$data,Model.ps$param.ps[ii,])
  }
  devc.mean <- mean(devc.mgt)
  ## desviación relativa en los parámetros medios D(E[par])+C
  psparmean <- apply(Model.ps$param.ps,2,mean,na.rm=TRUE)
  devc.psparmean <- (-2)*Model.ps$ermod$llikd(Model.ps$mod,Model.ps$data,psparmean)
  ## número efectivo de parámetros
  pd <- devc.mean-devc.psparmean
  dic <- devc.mean+pd
  return(list(pd=pd,dic=dic))
}

#' Title
#'
#' @param Model.ps ejecucion
#'
#' @returns pplc
#' @export
#'
#' @examples ejemplo1
pplc.modelps <- function(Model.ps){
  nm <- Model.ps$nm.ps
  nt <- length(Model.ps$data$tt)
  crpred<-matrix(NA,nrow=nm,ncol=nt)
  for(ii in 1:nm){
    for(it in 1:nt){
      crpred[ii,it] <- Model.ps$ermod$pred(Model.ps$data$tt[it],Model.ps$mod,Model.ps$param.ps[ii,])
      ## ver como corregir para descartar valores negativos
      ## crpred[ii,it] <- max(Model.ps$ermod$pred(Model.ps$data$tt[it],Model.ps$mod,Model.ps$param.ps[ii,]),0)
    }
  }
  crpred.mean <- apply(crpred,2,mean,na.rm=TRUE)
  crpred.var <- apply(crpred,2,var,na.rm=TRUE)
  pploss.g <- sum((crpred.mean - Model.ps$data$Cr)^2)
  pploss.p <- sum(crpred.var)
  pploss.pplc <- pploss.g+pploss.p
  return(list(g=pploss.g,p=pploss.p,pplc=pploss.pplc))
}
