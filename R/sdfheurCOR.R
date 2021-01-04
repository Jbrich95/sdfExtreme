#' Heuristic for creating bijective deformations (Correlation)
#'
#'This algorithm uses \code{FindSplineParFNCOR} to create a deformation with \code{m.init} inital anchor points.
#'The initial anchor points are the last \code{m.init} entries of \code{Full.m.ind}. After creating this deformation,
#'the spline values are used as initial parameters in creating a deformation with \code{m.init+1} anchor points.
#'This iterative procedure repeats until a deformation with all \code{Full.m.ind} anchor points is created.
#'
#'
#' @param m.init Number of inital anchor points. Must have \code{m.init < length(Full.m.ind)}.
#' @param Full.m.ind Full vector of indices for anchor points in \code{Gcoords}.
#' @param Gcoords A \code{d} by 2 matrix of G-plane coordinates.
#' @param emp.cor A \code{d} by \code{d} matrix of pairwise empirical correlation values.
#' @param type \code{"F-norm"} for Frobenius norm method using theoretical \eqn{\rho(h^{*}_{ij})} from a Matérn correlation model. \code{"Smith"} for original Smith (1996) method, also using Matérn correlation model.
#' @param n Number of data points used to estimate \code{emp.cor}. Only neccesary for \code{type=="Smith"}.
#' @param par Initial parameters for first deformation. If not stated, initial parameters are given.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?
#' @return List with three elements:\describe{
#' \item{par}{Spline parameter values.}
#' \item{value}{Objective value from final optimisation.}
#' \item{m.ind}{Vector of indices for full set of anchor points.}
#' }
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
#'
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
#'
#'m.init<-6
#'Full.m.ind<-sample(1:dim(Gcoords)[1],m.init+2)
#'
#'#Transform to D-plane using Smith (1996) method
#'
#'##WARNING: This may take a while to run.
#' sdf<-sdf.heur.Cor(m.init,Full.m.ind,Gcoords,emp.cor,type="Smith",n=dim(Z)[1],sphere.dis=TRUE)
#'
#'#Plot Dcoords
#'
#'Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
#'plot(Dcoords,main="D-plane",ylab="",xlab="")
#' @export
sdf.heur.Cor=function(m.init,Full.m.ind,Gcoords,emp.cor,type=c("F-norm","Smith"),n=0,par=NULL,sphere.dis=F){
  if(type!="F-norm" & type!="Smith"){
    stop("Invalid type")
  }
  if(m.init > (length(Full.m.ind)-1)){
    stop("Initial m too large")
  }

  m.final=length(Full.m.ind)
  m.ind=Full.m.ind[(m.final-m.init+1):m.final]

  reject=m.reject(m.ind,Gcoords,sphere.dis=sphere.dis)
  if(reject==1){
    stop("Initial anchor points rejected - May be too close together")
  }else print("Initial anchor points accepted")
  if(is.null(par)){
    par=c(0.5,0.2,1,0.5,rep(0,2*m.init-6))
  }
  if(type=="F-norm"){
    sdf<-optim(fn=FindSplineParFNCOR,par=par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead")

    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "BFGS")
    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead")
    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000),emp.cor=emp.cor,method = "BFGS")
  }else if(type =="Smith"){
    sdf<-optim(fn=FindSplineParFNCOR,par=par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000),emp.cor=emp.cor,method = "Nelder-Mead",n=n)

    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "BFGS",n=n)
    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead",n=n)
    sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
               control=list(maxit=2000), emp.cor=emp.cor,method = "BFGS",n=n)

  }
  print(paste("sdf with m = ", m.init," finished"))
  for(it in (m.init+1):m.final){

    m.ind=Full.m.ind[(m.final-it+1):m.final]
    tempnd=(length(sdf$par)-4)/2
    temppar=c(sdf$par[1:4],0,sdf$par[(5:(5+tempnd-1))],0,sdf$par[-c(1:(5+tempnd-1))])
    if(type =="F-norm"){
      sdf<-optim(fn=FindSplineParFNCOR,par=temppar,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead")

      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "BFGS")
      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead")
      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="F-norm",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000),emp.cor=emp.cor,method = "BFGS")
    }else if(type=="Smith"){

      sdf<-optim(fn=FindSplineParFNCOR,par=temppar,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead",n=n)

      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000),emp.cor=emp.cor,method = "BFGS",n=n)
      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "Nelder-Mead",n=n)
      sdf<-optim(fn=FindSplineParFNCOR,par=sdf$par,type="Smith",Gcoords=Gcoords,m.ind=m.ind,sphere.dis=sphere.dis,
                 control=list(maxit=2000), emp.cor=emp.cor,method = "BFGS",n=n)

    }
    print(paste("sdf with m = ", it," finished"))

  }
  sdf$m.ind=Full.m.ind
  return(sdf)
}
