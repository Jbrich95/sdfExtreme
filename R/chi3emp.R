#' Empirical estimate of triple-wise chi
#'
#' Calculates empirical estimate of triple-wise chi for \code{(U,V,W)}.
#'
#' @param U Vector of length \eqn{n} of standard uniform variables.
#' @param V Vector of length \eqn{n} of standard uniform variables.
#' @param W Vector of length \eqn{n} of standard uniform variables.
#' @param q Value of exceedance probability used in estimation.
#'
#' @return Estimate of triple-wise chi for \eqn{(U,V,W)}.
#'
#'
#' @examples
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
#'chi3.emp(Z_U[,ind.triple[1]],Z_U[,ind.triple[2]],Z_U[,ind.triple[3]],q)
#'

chi3.emp=function(U,V,W,q){
sum((U>=q)&(V>=q)&(W>=q))/(length(U)*(1-q))
}
