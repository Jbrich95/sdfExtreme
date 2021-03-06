#' Find spline parameters (Extremal)
#'
#'Find spline parameters for deformation using extremal dependence methods. Provides Frobenuis norm of difference between pairwise empircal and theoretical \eqn{\chi} measures.  Output to be minimised to estimate spline parameters.
#'
#' @param par Parameter values \eqn{\phi=(b_{1},b_{2},\rho,\kappa,\delta^{(1)}_{4}...\delta^{(1)}_{m},\delta^{(2)}_{4}...\delta^{(2)}_{m})}. If m < 4, \eqn{\delta} parameters are not needed.
#' @param Gcoords A \code{d} by 2 matrix of G-plane coordinates.
#' @param m.ind A vector of length \code{m < d} giving the indices of the anchor points in \code{Gcoords}.
#' @param emp.dep A \code{d} by \code{d} matrix of pairwise empirical chi values.
#' @param type \code{"CHI_BR"} for theoretical \eqn{\chi(h^{*}_{ij})} from a Brown-Resnick model. \code{"CHI_q_IBR"} for theoretical \eqn{\chi_q(h^{*}_{ij})} from an inverted Brown-Resnick model.
#' @param q Threhold for \eqn{\chi_q(h^{*}_{ij})}. Only needed if \code{type=="CHI_q_IBR"}.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?

#' @return Frobenius norm of difference between theoretical \eqn{\chi(h^{*}_{ij})} or \eqn{\chi_{q}(h^{*}_{ij})} matrix and \code{emp.chi}.
#' @examples
#' data("Aus_Heat")
#'
#' Z<-Aus_Heat$Temp.
#' Gcoords<-Aus_Heat$coords
#'
#' Z_U<-Z
#'
#' unif<-function(x){rank(x)/(length(x)+1)}
#'
#'#Transform to uniform margins
#' for(i in 1:dim(Z_U)[2]){
#'
#'   Z_U[,i]<-unif(Z[,i])
#' }
#'
#' chi.emp<-function(U,V,z) sum((U>=z)&(V>=z))/sum(U>=z)
#' q<-0.98
#'
#'#Calculate pairwise empirical chi
#' emp.chi<-matrix(rep(0,dim(Z_U)[2]^2),nrow=dim(Z_U)[2],ncol=dim(Z_U)[2])
#' for(i in 1:dim(Z_U)[2]){
#'   for(j in 1:i){
#'
#'     emp.chi[i,j]<-chi.emp(Z_U[,i],Z_U[,j],q)
#'   }
#'
#' }
#' emp.chi<-emp.chi+t(emp.chi)
#' diag(emp.chi)<-diag(emp.chi)/2
#'
#'# Set number of anchor points
#'m<-10
#'# Sample anchor points
#'m.ind<-sample(1:dim(Gcoords)[1],m,replace=FALSE)
#'
#'#Transform to D-plane using Brown-Resnick theoretical chi function
#' sdf<-optim(fn=FindSplineParFNEXTR,par=c(0.05,0.05,0,1,rep(0,2*m-6)),type="CHI_BR",sphere.dis=TRUE,
#'            Gcoords=Gcoords,m.ind=m.ind,control=list(maxit=2000), emp.dep=emp.chi,
#'	      method = "Nelder-Mead")
#'
#' sdf<-optim(fn=FindSplineParFNEXTR,par=sdf$par,sphere.dis=TRUE,type="CHI_BR",
#'            Gcoords=Gcoords,m.ind=m.ind,
#'            control=list(maxit=2000), emp.dep=emp.chi,method = "BFGS")
#' sdf$m.ind<-m.ind
#'
#'#Plot Dcoords
#'
#'Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
#'plot(Dcoords)
#' @export
FindSplineParFNEXTR=function(par,Gcoords,m.ind,emp.dep,type=c("CHI_BR","CHI_q_IBR"),q=0,sphere.dis=FALSE){
  if(type!="CHI_BR"&&type!="CHI_q_IBR"){
    stop("No type")
  }
  b1 = par[1]
  b2 = par[2]
  rho = par[3]
  theta_2=par[4]

  #Different constraints on theta_2 for SGP and IMSP
  if(type=="CHI_BR"||type=="CHI_q_IBR"){
    if(b1<0||b2<0||theta_2>2||theta_2<0.01){return(10e10)}
  }

  m=length(m.ind)
  coordm<-Gcoords[m.ind,]
  delta1<-delta2<-rep(0,3)
  if(m>3){
    delta1<-par[5:(length(par)-m+3)]
    delta2<-par[(length(par)-m+4):length(par)]


    dmat<-matrix(c(1,1,1,coordm[c(m-2,m-1,m),1],coordm[c(m-2,m-1,m),2]),3,3,byrow = TRUE)
    rhsd1<- -c(sum(delta1),sum(delta1*coordm[-c(m-2,m-1,m),1]), sum(delta1*coordm[-c(m-2,m-1,m),2]))
    rhsd2<- -c(sum(delta2),sum(delta2*coordm[-c(m-2,m-1,m),1]), sum(delta2*coordm[-c(m-2,m-1,m),2]))

    delta1final3<- solve(dmat)%*%rhsd1
    delta2final3<- solve(dmat)%*%rhsd2

    delta1<-c(delta1,delta1final3)
    delta2<-c(delta2,delta2final3)
  }

  etafunc<-function(Gcoordvec,s1,s2,sphere.dis=FALSE)
  {
    if(sphere.dis==FALSE){
      r<-sqrt((s1-Gcoordvec[1])^2+(s2-Gcoordvec[2])^2)
    }else{
      r=rdist.earth(rbind(c(s1,s2),Gcoordvec),miles=FALSE)[1,2]

    }
    if(r==0){return(0)}
    return(r^2 * log(r))
  }

  f1<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b1^2*s1 + rho*b1*b2*s2 + sum(delta1*eta)
  }
  f2<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b2^2*s2 + rho*b1*b2*s1 + sum(delta2*eta)
  }
  f<-function(s1s2)
  {
    s1<-s1s2[1]
    s2<-s1s2[2]
    return(c(f1(s1,s2),f2(s1,s2)))
  }

  Dcoord<-t(apply(Gcoords,1,f))
  if(sphere.dis==FALSE){
    vario_H=(rdist(Dcoord))^theta_2
  }else{
    vario_H=(rdist.earth(Dcoord,miles=FALSE))^theta_2

  }
  theta_h=2*pnorm(sqrt(vario_H/2))

  if(type =="CHI_BR"){
    dep.mat=2-theta_h
  }else if(type=="CHI_q_IBR"){
    dep.mat=(1-q)^(theta_h-1)
  }


  dis = norm(dep.mat-emp.dep,type="F")
  return(dis)
}
