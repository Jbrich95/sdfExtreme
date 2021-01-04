#' Theoretical triple-wise \eqn{\chi}  (\eqn{\chi_q}) for a (Inverted) Brown-Resnick process
#'
#' Calculates the theoretical triple-wise \eqn{\chi(s_i,s_j,s_k)} or \eqn{\chi_q(s_i,s_j,s_k)} for a given 3 by 3 matrix of semivariogram. values.

#' @param v_H 3 by 3 matrix of semivariogram. values.
#' @param q Exceedance threshold for \eqn{\chi_q(s_i,s_j,s_k)}.
#' @return Theoretical \eqn{\chi(s_i,s_j,s_k)} or \eqn{\chi_q(s_i,s_j,s_k)} measure for the corresponding matrix of semivariogram values.
#'
#'
#' @examples
#'
#' data(Aus_Heat)
#' data(Aus_Heat_Output)
#' Z<-Aus_Heat$Temp.
#' Gcoords<-Aus_Heat$coords
#' sdf<-Aus_Heat_Output$sdf
#' likD.MSP<-Aus_Heat_Output$likD.MSP
#' likD.IMSP<-Aus_Heat_Output$likD.IMSP
#'
#' Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
#'
#'  ind.triple<-c(1,2,3) ##Denotes indices of triple for which triple-wise chi calculated.
#'
#'  #MSP
#'  v_H<-(rdist.earth(Dcoords[ind.triple,],miles=F)/likD.MSP$par[2])^likD.MSP$par[1]

#'  print(chi3(v_H))
#'
#'   #IMSP - likD.IMSP only contains range parameter as smoothing parameter is 2.
#'  v_H<-(rdist.earth(Dcoords[ind.triple,],miles=F)/likD.IMSP$par[1])^2

#'  print(chi3q(v_H,q=0.95))
#'
#' @rdname chi3
#' @export
chi3=function(v_H){
  return(3-theta2(v_H[1,2])-theta2(v_H[1,3])-theta2(v_H[2,3])+theta3(v_H))
}
#' @rdname chi3
#' @export
chi3q=function(v_H,q){
  return((1-q)^(theta3(v_H)-1))
}

theta3=function(v_H){
  corr=rep(0,3)
  corr[1]=0.5*(v_H[1,2]+v_H[1,3]-v_H[2,3])/(v_H[1,2]*v_H[1,3])^0.5
  corr[2]=0.5*(v_H[1,2]+v_H[2,3]-v_H[1,3])/(v_H[1,2]*v_H[2,3])^0.5
  corr[3]=0.5*(v_H[3,2]+v_H[1,3]-v_H[2,1])/(v_H[3,2]*v_H[1,3])^0.5

  theta_H=pmvnorm(lower=-Inf,upper=c((v_H[1,2]/2)^0.5,(v_H[1,3]/2)^0.5),mean=c(0,0),corr=matrix(c(1,corr[1],corr[1],1),nrow=2,ncol=2))[1]+
    pmvnorm(lower=-Inf,upper=c((v_H[1,2]/2)^0.5,(v_H[2,3]/2)^0.5),mean=c(0,0),corr=matrix(c(1,corr[2],corr[2],1),nrow=2,ncol=2))[1]+
    pmvnorm(lower=-Inf,upper=c((v_H[1,3]/2)^0.5,(v_H[2,3]/2)^0.5),mean=c(0,0),corr=matrix(c(1,corr[3],corr[3],1),nrow=2,ncol=2))[1]
  return(theta_H)
}
theta2=function(v_H){
  2*pnorm(sqrt(v_H/2))
}

