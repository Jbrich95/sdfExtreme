#' Negative composite censored log-likelihood for the stationary inverted Brown-Resnick process
#' 
#' Calculate the negative composite censored log-likelihood for the stationary inverted Brown-Resnick process on standard exponential margins.
#'
#' @param par Stationary semivariogram parameters \eqn{(\kappa,\lambda)}.
#' @param ZBE A list of length \code{d} with each element a matrix with 2 columns.
#' @param ZXE A list of length \code{d} with each element a matrix with 2 columns.
#' @param ZYE A list of length \code{d} with each element a matrix with 2 columns.
#' @param ZNE A list of length \code{d} with each element a matrix with 2 columns.
#' @param coord A \code{d} by 2 matrix of coordinates.
#' @param n.a A vector with length determined by the number of unique pairs in \code{ZNE}.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?

#' @return Negative composite censored log-likelihood for inverted Brown-Resnick process.
#' @examples
#' # For a N by d matrix of data "Z" and d by 2 matrix of coordinates "Gcoords".
#' # We use a very small subset of data(Aus_Heat) as an example.
#'##THIS WILL TAKE A LONG TIME TO RUN WITH THE FULL DATASET##
#'
#' library(fields)
#' data(Aus_Heat) 
#' Z<-Aus_Heat$Temp.[,1:3]
#' Gcoords<-Aus_Heat$coords[1:3,]
#'
#' unif<-function(x) rank(x)/(length(x)+1)
#' Z_U<-Z
#' for(i in 1:dim(Z_U)[2]) Z_U[,i]<-unif(Z[,i]) # Transform to uniform margins
#' Z_Exp<-qexp(Z_U) #Transform to exponential margins
#'
#' q<-0.98
#' u<-quantile(Z_Exp,prob=q) # Censoring threshold
#' # Create a list of length 4 of pairwise exceedances and index in coordinate matrix
#'Zpair<-makepairs(Z_Exp,u=u) 
#'
#'#Identify number of non-exceedances per pair
#'#Speeds up likelihood computation
#' dist<-rdist(Gcoords+runif(dim(Gcoords)[1]*2,0,1))
#' dist2<-apply(Zpair[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
#' unique.d<-as.matrix(unique(dist2))
#' n.a<-rep(0, dim(unique.d)[1])
#' for(i in 1:length(unique.d)){
#'   n.a[i]<-sum(dist2==unique.d[i,1])
#' }
#'
#'#WARNING- pairwise likelihoods take some time to run
#'# Inverted Brown-Resnick process fit to G-plane coordinate system
#'likG.IMSP<-optim(fn=nllIMSPexp,par=c(1,500),ZBE=as.matrix(Zpair[[1]]),
#'            ZXE=as.matrix(Zpair[[2]]),
#'            ZYE=as.matrix(Zpair[[3]]),ZNE=as.matrix(Zpair[[4]]),
#'            n.a=n.a,sphere.dis=TRUE,coord=Gcoords,
#'            control=list(maxit=2000),
#'            method = "Nelder-Mead",hessian=TRUE)

nllIMSPexp=function(par,ZBE,ZXE,ZYE,ZNE,coord,n.a,sphere.dis=F){


  smooth=par[1]
  range=par[2]

  #Different constraints on theta_2 for SGP and IMSP

  if(smooth>2||smooth<0.01||range<0.01){return(10e10)}
  if(sphere.dis==F){
    vario_H=(rdist(coord)/range)^smooth
  }else{
    vario_H=(rdist.earth(coord,miles=F)/range)^smooth


  }
  #Both Exceed
  if(dim(ZBE)[1]!=0){
    var=apply(matrix(ZBE[,3:4],ncol=2,byrow=F),1,function(x){vario_H[x[1],x[2]]})

    a=(2*var)^0.5
    ZBE=ZBE[,-c(3,4)]
    ZBE=cbind(ZBE,a)
    pBE=apply(ZBE,1,function(x){

      x[1]=1/x[1]
      x[2]=1/x[2]
      V=Vterm(x[1],x[2],x[3])
      V1=Vpart(x[1],x[2],x[3])
      V2=Vpart(x[2],x[1],x[3])
      V12=Vpart2(x[1],x[2],x[3])
      2*log(x[1]*x[2])+log(V1*V2-V12)-V
    }
    )
  }else{pBE=0}

  if(dim(ZXE)[1]!=0){
    #X Exceed
    var=apply(matrix(ZXE[,3:4],ncol=2,byrow=F),1,function(x){vario_H[x[1],x[2]]})

    a=(2*var)^0.5
    ZXE=ZXE[,-c(3,4)]
    ZXE=cbind(ZXE,a)
    pXE=apply(ZXE,1,function(x){

      V=Vterm(1/x[1],1/x[2],x[3])
      V1=Vpart(1/x[1],1/x[2],x[3])
      log(exp(-x[1])+x[1]^(-2)*(V1)*exp(-V))
    }
    )
  }else{pXE=0}

  if(dim(ZYE)[1]!=0){
    #Y Exceed
    var=apply(matrix(ZYE[,3:4],ncol=2,byrow=F),1,function(x){vario_H[x[1],x[2]]})

    a=(2*var)^0.5
    ZYE=ZYE[,-c(3,4)]
    ZYE=cbind(ZYE,a)
    pYE=apply(ZYE,1,function(x){


      V=Vterm(1/x[1],1/x[2],x[3])
      V1=Vpart(1/x[2],1/x[1],x[3])
      log(exp(-x[2])+x[2]^(-2)*(V1)*exp(-V))
    }
    )
  }else{pYE=0}
  #No exceed

  # EDITED
  if(dim(ZNE)[1]!=0){
    var=apply(unique(matrix(ZNE[,c(3,4)],ncol=2,byrow=F)),1,function(x){vario_H[x[1],x[2]]})

    a=(2*var)^0.5
    unique.a=matrix(a,ncol=1)

    u=as.numeric(ZNE[1,1])
    pNE=apply(unique.a,1,function(x){
      V=Vterm(1/u,1/u,x[1])
      log(1-2*(exp(-u))+exp(-V))
    })
    pNE=n.a*pNE
  }else{pNE=0}

  loglik=sum(pBE)+sum(pXE)+sum(pYE)+sum(pNE)

  if(is.finite(loglik)){
    return(-loglik)
  }else{return(1e10)}
}

