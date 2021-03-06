#' Find spline parameters (Correlation)
#'
#'Find spline parameters for deformation based on correlation methods. Either use the Frobenuis norm via \code{type=="F-norm"} or the original Smith (1996) method via \code{type=="Smith"}. Output to be minimised to estimate spline parameters.
#'
#' @param par Parameter values \eqn{\phi=(b_{1},b_{2},\rho,\kappa,\delta^{(1)}_{4}...\delta^{(1)}_{m},\delta^{(2)}_{4}...\delta^{(2)}_{m})}. If \eqn{m < 4}, \eqn{\delta} parameters are not needed.
#' @param Gcoords A \code{d} by 2 matrix of G-plane coordinates.
#' @param m.ind A vector of length \code{m<d} giving the indices of the anchor points in \code{Gcoords}.
#' @param emp.cor A \code{d} by \code{d} matrix of pairwise empirical correlation values.
#' @param type \code{"F-norm"} for Frobenius norm method using theoretical \eqn{\rho(h^{*}_{ij})} from a stationary Matérn correlation model. \code{"Smith"} for original Smith (1996) method; also uses stationary Matérn correlation model.
#' @param n Number of data points used to estimate \code{emp.cor}. Only neccesary for \code{type=="Smith"}.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?

#' @return Frobenius norm of difference between theoretical \eqn{\rho(h^{*}_{ij})} matrix and \code{emp.cor}, or negative log-likelihood for a stationary Gaussian process fit using D-plane coordinates.
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
#'#Transform to Gaussian margins
#'
#' Z_N<-qnorm(Z_U)
#'
#'#Calculate pairwise empirical correlation
#' emp.cor<-matrix(rep(0,dim(Z_N)[2]^2),nrow=dim(Z_N)[2],ncol=dim(Z_N)[2])
#' for(i in 1:dim(Z_N)[2]){
#'   for(j in 1:i){
#'
#'     emp.cor[i,j]<-cor(Z_N[,i],Z_N[,j])
#'   }
#'
#' }
#' emp.cor<-emp.cor+t(emp.cor)
#' diag(emp.cor)<-diag(emp.cor)/2
#'
#'# Set number of anchor points
#'m<-10
#'# Sample anchor points
#'m.ind<-sample(1:dim(Gcoords)[1],m,replace=FALSE)
#'
#'#Transform to D-plane using Frobenius norm method
#' sdf<-optim(fn=FindSplineParFNCOR,par=c(0.05,0.05,0,1,rep(0,2*m-6)),
#'            type="F-norm",sphere.dis=TRUE,Gcoords=Gcoords,m.ind=m.ind,
#'            control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead")
#'
#' sdf<-try(optim(fn=FindSplineParFNCOR,par=sdf$par,sphere.dis=TRUE,type="F-norm",
#'          Gcoords=Gcoords,m.ind=m.ind, control=list(maxit=2000), emp.cor=emp.cor,
#'	    method = "BFGS"))
#'sdf$m.ind<-m.ind
#'#Plot Dcoords
#'
#'Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
#'plot(Dcoords)
#' @export
FindSplineParFNCOR=function(par,Gcoords,m.ind,n=0,emp.cor,type=c("F-norm","Smith"),sphere.dis=FALSE){
  if(type!="Smith"&&type!="F-norm"){
    stop("No type")
  }
  b1 = par[1]
  b2 = par[2]
  rho = par[3]
  theta_2=par[4]

  if(b1<0||b2<0||theta_2>3.5||theta_2<0.01){return(10e10)}

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
      r<-rdist.earth(rbind(c(s1,s2),Gcoordvec),miles=FALSE)[1,2]

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
    cor.mat<-Matern(rdist(Dcoord),nu=theta_2,phi=1)
  }else{
    cor.mat<-Matern(rdist.earth(Dcoord,miles=FALSE),nu=theta_2,phi=1)


  }
  if(type =="F-norm"){
    dis = norm(cor.mat-emp.cor,type="F")
    return(dis)}else{
      if(det(cor.mat)==0){return(1e10)}

      nll<- (n/2)*log(det(cor.mat))+((n-1)/2)*sum(diag(solve(cor.mat)%*%emp.cor))
      return(nll)
    }

}
