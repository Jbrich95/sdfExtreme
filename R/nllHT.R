#' Negative log-likelihood for bivariate Heffernan and Tawn (2004) model
#'
#'Calculates the negative log-likelihood for \eqn{Y|X=x}. This is a profile log-likelihood
#'with free parameters \eqn{\alpha} and \eqn{\beta}. We use the same normalising functions as in the original paper.
#' @param X Vector of Laplace variables. This is the conditioning variable i.e. \eqn{X=x}.
#' @param Y Vector of Laplace variables, \eqn{Y}.
#' @param par Vector \eqn{(\alpha,\beta)}.
#' @return Negative log-likelihood for Heffernan and Tawn model fit to \eqn{Y|X=x}.
#'
#' @references Heffernan and Tawn (2004), Journal of the Royal Statistical Society: Series B, 66:497-546,
#' (\href{https://doi.org/10.1111/j.1467-9868.2004.02050.x}{doi})
#' @examples
#' #For a N by d matrix of data "Z" and d by 2 matrix of coordinates "Gcoords".
#' # We use data(Aus_Heat) as an example.
#'
#' library(rmutil)
#' library(fields)
#'data(Aus_Heat) 
#' Z<-Aus_Heat$Temp.
#' Gcoords<-Aus_Heat$coords
#'
#' unif<-function(x) rank(x)/(length(x)+1)
#' Z_U<-Z
#' for(i in 1:dim(Z_U)[2]) Z_U[,i]<-unif(Z[,i]) # Transform to uniform margins
#' Z_LP<-qlaplace(Z_U) # Transform to Laplace margins
#'

#'
#'p<-dim(Gcoords)[1]
#'ConExp<-matrix(0,nrow=p,ncol=p)
#'
#'u<-quantile(Z_LP,0.98) #Quantile for estimating conditional expectation
#'
#'#Calculate conditional expectation for each pair
#'for(i in 1:p){
#'  for(j in 1:i){
#'    Exceedances<-cbind(Z_LP[,i],Z_LP[,j])[which(Z_LP[,i]>=u),]
#'    opt<-optim(nllHT,X=Exceedances[,1],Y=Exceedances[,2],par=c(0.3,0.8))
#'    #print(opt)
#'    alpha<-opt$par[1]
#'    beta<-opt$par[2]
 #'   mu<-mean((Exceedances[,2]-alpha*Exceedances[,1])/Exceedances[,1]^beta)
 #'   ConExp[i,j]<-alpha*u+u^beta*mu
 #' }
#' 
#'}
#'ConExp<-ConExp+t(ConExp)  ##Symmetry assumed
#'diag(ConExp)<-diag(ConExp)/2
#'
#'plot(rdist.earth(Gcoords,miles=F),ConExp, ylab="Conditional Expectation",
#'     xlab="Distance (Km)", main="")

nllHT=function(X,Y,par){

  alpha=par[1]
  beta=par[2]


  if(alpha < 0 || alpha > (1-10^(-6))  || beta > 1-10^(-6) || beta < 0){return(10e10)}

  mu=mean((Y-alpha*X)/X^beta)
  sig=(mean((Y-mu*X^beta-alpha*X)^2/X^(2*beta)))^0.5
  if(sig<=0){return(10e10)}
  negloglik <- -sum(dnorm(Y,alpha*X+mu*X^beta,sig*X^beta,log=TRUE))
  if(is.finite(negloglik)){
    return(negloglik)
  }else{return(1e10)}
}
