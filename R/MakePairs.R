#' Make pairs
#' 
#'   Produce a list of all configurations of pairwise exceedances above some threshold \code{u}. Used in censored likelihood estimation.
#'
#' @param Z A \code{N} by \code{d} matrix.
#' @param u Censoring threshold.

#' @return List of four data frames. Each data frame has 4 columns and the rows sum to \code{N*d*(d-1)/2}. For each data frame the first two columns are the values and the last two are the indices for the location of these values in the columns of \code{Z}.
#' \describe{
#' \item{Component 1}{Pairs of observations where both components exceed \code{u}. }
#' \item{Component 2}{Pairs of observations where the first component exceeds \code{u} and the second equals \code{u}. }
#' \item{Component 3}{Pairs of observations where the second component exceeds \code{u} and the first equals \code{u}. }
#' \item{Component 4}{Pairs of observations where both components equal \code{u}. }
#'}
#' @examples
#' #For a N by d matrix of data "Z" and d by 2 matrix of coordinates "Gcoords".
#' # We use a very small subset of data(Aus_Heat) as an example.
#'##THIS WILL TAKE A LONG TIME TO RUN WITH THE FULL DATASET##
#'
#'data(Aus_Heat) 
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

makepairs=function(Z,u){
  out=list()
  N=dim(Z)[1]
  M=dim(Z)[2]
  allz=c()
  for(i in 1:N){
    allz=rbind(allz,cbind(expand.grid(Z[i,],Z[i,])))
  }
  allz=cbind(allz,rep(1:M,M*N),rep(1:M,each=M,N))
  names(allz)=c()
  allz=subset(allz,allz[,3]<allz[,4])


  allzBE=subset(allz, allz[,1]>u & allz[,2]>u)
  allzYE=subset(allz, allz[,1]<=u & allz[,2]>u)
  if(!is.na(allzYE[1,1])){
    allzYE[,1]=u
  }
  allzXE=subset(allz, allz[,1]>u & allz[,2]<=u)
  if(!is.na(allzXE[1,1])){
    allzXE[,2]=u
  }



  allzNE=subset(allz, allz[,1]<=u & allz[,2]<=u)
  if(!is.na(allzNE[1,1])){

    allzNE[,1]=u
    allzNE[,2]=u
  }
  out[[1]]=allzBE
  out[[2]]=allzXE
  out[[3]]=allzYE
  out[[4]]=allzNE
  return(out)

}
