#' Simulate non-stationary Brown-Resnick process
#'
#'Simulate a non-stationary (or stationary) Brown-Resnick process using the 'stopping rule'
#'methodology of Dieker and Mikosch (2015). Non-stationary Brown-Resnick processes are simulated
#'using the non-stationary semivariogram. detailed in the help file for \code{gamsv}.
#'
#'
#' @param reps Number of realisations, \eqn{n}.
#' @param locs A \eqn{p} by 2 matrix of coordinates.
#' @param lambda Value of \eqn{\lambda}. See \code{help(gamsv)}.
#' @param kappa Value of \eqn{\kappa}. See \code{help(gamsv)}.
#' @param centre Vector of length 2 giving coordinates of centre of non-stationarity for semivariogram. If \code{NULL}, stationary semivariogram used. See \code{help(gamsv)}.
#' @param tau A \eqn{p} by \eqn{p} matrix of correlation values.
#' @return \eqn{n} by \eqn{p} matrix
#'
#'@references Dieker and Mikosch (2015) Extremes, 18(2):301-314,
#' (\href{https://doi.org/10.1007/s10687-015-0214-4}{doi})
#'
#' @examples
#' ##Creating correlation values to simulate non-stationary Brown-Resnick process.
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
#' gamsv(s1=sim.coords[j,],s2=c(0,0),lam=lambda,kap=kappa,centre=NULL)-
#' gamsv(s1=sim.coords[i,],s2=sim.coords[j,],lam=lambda,kap=kappa,centre=centre)
#'}
#'}
#'
#'##Simulates 10 realisations of non-stationary BR process. 
#'
#'Sim<-brnsims(reps=10,locs=sim.coords,kappa=kappa,lambda=lambda,centre=centre,tau=tau)
#'
#'
#' @rdname brnsims
#' @export
brnsims<-function(reps,locs,kappa,lambda,tau,centre=NULL){

 brsim<-function(kappa,lambda,B,locs,mu=1,tau,centre=NULL){

brspec<-function(kappa,lambda=lambda,locs,tau,centre=NULL){

  p<-dim(locs)[1]

  W<-gamma<-Y<-Z<-rep(NA,p)

  i<-sample(seq(1:p),size=1)   ## Sample index of a random site


  mu<-rep(0,p)  ## Mean vector for GP to simulate from.

  Z<-mvrnorm(1,mu,Sigma=tau)  ## Simulate a GP

  for(j in 1:p){

    gamma[j]<-gamsv(s1=locs[j,],s2=locs[i,], lam=lambda,kap=kappa,centre=centre)

    Y[j]<-Z[j]-gamma[j]  ##Subtract \gamma_j from Y_j

  }

  for(j in 1:p){

    W[j]<-exp(Y[j])/((1/p)*sum(exp(Y)))  ##Take mean

  }

  list(Ws=W) ## Outputs needed for input into next part of algorithm


}


rho<-0 ##Initialise rho

Z<-rep(0,dim(locs)[1])


W_i<-rep(NA,dim(locs)[1])



for(c in 1:100000){ ## Arbitrary loop

  rho<-rho+rexp(1)

  W_i<-brspec(kappa=kappa,lambda=lambda,locs=locs,tau=tau,centre=centre)$Ws ## Calling spectral function to obtain W.


  Z<-pmax(Z,W_i/rho)

  if((B/rho)<min(Z)){

    break

  }

}

Z<-(1/mu)*Z

list(Z=Z)

}


  brsims<-matrix(NA,nrow=reps,ncol=dim(locs)[1])

  p<-dim(locs)[1]

  B=length(locs)
  for(n in 1:reps){

    brsims[n,]<-brsim(kappa=kappa,lambda=lambda,B=B,locs=locs,tau=tau,centre=centre)$Z


  }

  return(brsims)

}




