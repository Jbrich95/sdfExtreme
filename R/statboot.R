#' Stationary Bootstrap
#'
#' Creates one bootstrap sample of a \eqn{N} by \eqn{d} matrix, \code{X}. Here \eqn{N} denotes the number of time points
#' and \eqn{d} the number of sampling locations. The stationary bootstrap (Politis and Romano, 1994) generates a bootstrap sample
#' by repeated sampling of random blocks with expected value \code{block.size} until a sample of length \eqn{N} has been created.
#'
#' @param X Matrix with \eqn{N} rows corresponding to time points and \eqn{d} columns corresponding to spatial locations.
#' @param mean.block.size Expected value of temporal block size.
#' 
# @return
#'
#'@references Politis and Romano (1994), JASA,5(4):303-336,
#' (\href{https://doi.org/10.1080/01621459.1994.10476870}{doi})
#'
#'
#' @examples \eqn{N} by \eqn{d} matrix
#'
#' data(Aus_Heat)
#' Z<-Aus_Heat$Temp.
#'
#' unif<-function(x) rank(x)/(length(x)+1)
#' Z_U<-Z
#' for(i in 1:dim(Z_U)[2]) Z_U[,i]<-unif(Z[,i]) # Transform to uniform margins
#'
#' q<-0.95
#'
#' ind.triple<-c(1,2,3) ##Denotes indices of triple for which triple-wise chi calculated.
#'
#' block.mean<-14 #Mean block size for random block choice - Here a fortnight
#'
#'boot <- stat.boot(Z_U[,ind.triple],block.mean)
#'
#'#Create one bootstrap estimate of triple-wise chi
#'chi3.emp(boot[,1],boot[,2],boot[,3],q)
#'

stat.boot=function(X,mean.block.size){
  N=dim(X)[1]
  
  
    #Wraps around
  tempX=rbind(X,X)
  
  ind=sample(1:N,1)
  b=0
  while(b==0){
   ##Ensure block size is not 0
   b=rgeom(1,1/(mean.block.size+1))

  }
  boot=matrix(tempX[(ind):(ind+b-1),],nrow=b)
 while(dim(boot)[1]<N){
    ind=sample(1:N,1)
    b=0
    while(b==0){
      ##Ensure block size is not 0
      b=rgeom(1,1/(mean.block.size+1))
      
    }
    boot=rbind(boot,tempX[(ind):(ind+b-1),])

}
  
boot=boot[1:N,]
  return(boot)
}
