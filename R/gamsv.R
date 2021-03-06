#' Non-stationary 'centered' semivariogram.
#'
#' Calculates the Fouedijo et al. (2015) theoretical semivariogram. This is a function of some central
#' location \eqn{o=} \code{centre} and distance between locations \code{s1} and \code{s2}.
#' The semivariogram. has the form \deqn{\gamma^*(s_2,s_1)=\gamma(||\psi(s_2)-\psi(s_1)||),} where
#' \deqn{\psi(s)=o+(s-o)||s-o||} and \deqn{\gamma(x)=(x/\lambda)^\kappa} for \eqn{\lambda > 0} and \eqn{\kappa \in (0,2]}.
#'
#' @param s1 Vector of length 2 giving coordinates of first location.
#' @param s2 Vector of length 2 giving coordinates of second location.
#' @param centre Vector of length 2 giving coordinates of centre of non-stationarity. If \code{NULL}, returns stationary semivariogram.
#' @param lam Value of \eqn{\lambda}.
#' @param kap Value of \eqn{\kappa}.
#' @return Non-stationary semivariogram. between locations \code{s1} and \code{s2}
#'
#'@references Fouedijo et al. (2015) Spatial Statistics, 13:45-61,
#' (\href{https://doi.org/10.1016/j.spasta.2015.05.001}{doi})
#'
#' @examples
#'
#' ##Creating correlation values to simulate non-stationary Brown-Resnick process. See help(brnsims).
#'
#'lambda<-2
#'centre<-c(0,0)
#'kappa<-0.8
#'
#'n.grid<-8
#'sim.coords<-as.matrix(expand.grid(seq(-1,1,length=n.grid),seq(-1,1,length=n.grid)))
#'
#'p<-dim(sim.coords)[1]
#'tau<-matrix(NA,nrow=p,ncol=p)
#'
#'for(i in 1:p){
#'for(j in 1:p){
#'tau[i,j]<-gamsv(s1=sim.coords[i,],s2=c(0,0),lam=lambda,kap=kappa,centre=NULL)+
#'gamsv(s1=sim.coords[j,],s2=c(0,0),lam=lambda,kap=kappa,centre=NULL)-
#'gamsv(s1=sim.coords[i,],s2=sim.coords[j,],lam=lambda,kap=kappa,centre=centre)
#'}
#'}
#'
#'##Simulates 10 realisations of non-stationary BR process. 
#'
#'Sim<-brnsims(reps=10,locs=sim.coords,kappa=kappa,lambda=lambda,centre=centre,tau=tau)
#'
#'


gamsv<-function(s1,s2,centre=NULL,lam,kap){
  Euc.dist<-function(s1,s2=NULL){
    if(is.null(s2)){s2=c(0,0)}
    return( sqrt((s1[1]-s2[1])^2+(s1[2]-s2[2])^2))

  }
  if(is.null(centre)){
    h=Euc.dist(s1,s2)
  }else{
    fx=Euc.dist(s1,centre)*(s1-centre)
    fy=Euc.dist(s2,centre)*(s2-centre)
    h=Euc.dist(fx,fy)
  }
  (h/lam)^kap

}
