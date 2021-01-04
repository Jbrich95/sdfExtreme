#' Score function estimation
#'
#'For \code{Score_Est}, the score for the likelihood of a Brown-Resnick or inverted Brown-Resnick model
#'is estimated at each time point. This can then be used to calculate the CLAIC, see examples.
#'\code{Score_Est_Smith} estimates the score for the likelihood of a Smith or inverted Smith model.
#'
#' @param Z_Exp A \code{N} by \code{d} matrix of exponential random variables.
#' @param u Censoring threshold.
#' @param coord A \code{d} by 2 matrix of coordinates.
#' @param lik Optim output from model fitting. See help(\code{nllMSPexp}) or help(\code{nllMSPexpSmith}).
#' @param type Either Brown-Resnick (\code{"BR"}) or inverted Brown-Resnick (\code{"IBR"}) for \code{Score_Est}. Smith or inverted Smith for \code{Score_Est_Smith}.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?
#' @return An \code{N} by 2 matrix of score estimates for semivariogram parameters \eqn{(\kappa,\lambda)} (just \eqn{\lambda} for \code{Score_Est_Smith}).
#' @examples
#'# Z_Exp: N by d matrix of data on exponential margins
#'# Gcoords and Dcoords: d by 2 matrix of sampling locations in the G-plane
#'# and D-plane respectively. See help(returnDcoord) for obtaining Dcoords.
#'# likG.MSP and likD.MSP: optim outputs from model fitting. See help(nllMSPexp).
#'
#'#We use data(Aus_Heat) as an example.
#'##THIS WILL TAKE A LONG TIME TO RUN WITH THE FULL DATASET##
#'
#'data(Aus_Heat) 
#' Z<-Aus_Heat$Temp.
#' Gcoords<-Aus_Heat$coords
#' 
#'data(Aus_Heat_Output)
#'sdf<-Aus_Heat_Output$sdf
#'likG.MSP<-Aus_Heat_Output$likG.MSP
#'likD.MSP<-Aus_Heat_Output$likD.MSP
#'
#'Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
#'
#'  unif<-function(x) rank(x)/(length(x)+1)
#' Z_U<-Z
#' for(i in 1:dim(Z_U)[2]) Z_U[,i]<-unif(Z[,i]) # Transform to uniform margins
#' Z_Exp<-qexp(Z_U) #Transform to exponential margins
#' 
#' q <- 0.98 # 98% quantile used as threshold in composite likelihood
#' u <- quantile(Z_Exp,prob=q)
#' 
#' 
#' s_G_MSP <- Score_Est(Z_Exp,u,coord=Gcoords,lik=likG.MSP,type="BR",sphere.dis=T)
#' s_D_MSP <- Score_Est(Z_Exp,u,coord=Dcoords,lik=likD.MSP,type="BR",sphere.dis=T)
#'
#' #Estimate variance of score - This is specfic to the Australian summer temperatures data.
#' # Set block.size. Here we take 90 and 91, corresponding to a regular season
#' block.sizes <- c(91,90) 
#' #and a season with a leap year
#' years <- 1957:2014
#' k <- length(years)
#'
#'temp <- matrix(0,nrow=k,ncol=2)
#'temp2 <- matrix(0,nrow=k,ncol=2)
#'int <- 0
#'for(l in 1:k){
#'  if(years[l]%%4==0) block.size=block.sizes[1] else block.size=block.sizes[2]
#'
#'  int <- int + block.size
#'  temp[l,] <- colSums(s_G_MSP[(int-block.size+1):(int),])
#'  temp2[l,] <- colSums(s_D_MSP[(int-block.size+1):(int),])
#'
#'}
#'
#'#Estimate variance of score
#'
#'varS_G_MSP <- var(temp)
#'varS_D_MSP <- var(temp2)
#'
#'
#'#Estimate CLAIC
#'CLAIC_G_MSP <-2*likG.MSP$value+2*sum(diag(varS_G_MSP%*%solve(likG.MSP$hessian)))
#'CLAIC_D_MSP <- 2*likD.MSP$value+2*sum(diag(varS_D_MSP%*%solve(likD.MSP$hessian)))
#'
#'#Estimate CLAIC for Inverted Smith model
#'#'# likG.IMSP and likD.IMSP: optim outputs from model fitting. See help(nllMSPexp).
#'
#'
#'likG.IMSP<-Aus_Heat_Output$likG.IMSP
#'likD.IMSP<-Aus_Heat_Output$likD.IMSP
#'
#' q <- 0.98 # 98% quantile used as threshold in composite likelihood
#' u <- quantile(Z_Exp,prob=q)
#' s_G_IMSP <- Score_Est_Smith(Z_Exp,u,coord=Gcoords,lik=likG.IMSP,type="InvSmith",sphere.dis=T)
#' s_D_IMSP <- Score_Est_Smith(Z_Exp,u,coord=Dcoords,lik=likD.IMSP,type="InvSmith",sphere.dis=T)
#'
#' #Estimate variance of score - This is specfic to the Australian summer temperatures data.
#' block.sizes <- c(91,90) #Set block.size # Here we take 90 and 91, corresponding to a regular season
#' #and a season with a leap year
#' years <- 1957:2014
#' k <- length(years)
#'
#'temp <- matrix(0,nrow=k,ncol=1)
#'temp2 <- matrix(0,nrow=k,ncol=1)
#'int <- 0
#'for(l in 1:k){
#'  if(years[l]%%4==0) block.size=block.sizes[1] else block.size=block.sizes[2]
#'
#'  int <- int + block.size
#'  temp[l] <- sum(s_G_IMSP[(int-block.size+1):(int)])
#'  temp2[l] <- sum(s_D_IMSP[(int-block.size+1):(int)])
#'
#'}
#'
#'#'#Estimate variance of score
#'
#'varS_G_IMSP <- var(temp)
#'varS_D_IMSP <- var(temp2)
#'
#'
#'#Estimate CLAIC
#'CLAIC_G_IMSP <-2*likG.IMSP$value+2*sum(diag(varS_G_IMSP%*%solve(likG.IMSP$hessian)))
#'CLAIC_D_IMSP <- 2*likD.IMSP$value+2*sum(diag(varS_D_IMSP%*%solve(likD.IMSP$hessian)))
#'
#' @rdname Score_Est
#' @export
Score_Est=function(Z_Exp,u,coord,lik,type=c("BR,IBR"),sphere.dis=F){

  N=dim(Z_Exp)[1]
  s=array(0,dim=c(N,2))
  for( i in 1:N){
    one.timepoint=t(as.matrix(Z_Exp[i,]))
    temppairs=makepairs(Z=one.timepoint,u=u)
    dist=rdist(coord+runif(dim(coord)[1],0,1))
    dist2=apply(temppairs[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
    unique.d=as.matrix(unique(dist2))
    n.a=rep(0, dim(unique.d)[1])
    if(length(n.a)!=0){
      for(j in 1:length(unique.d)){
        n.a[j]=sum(dist2==unique.d[j,1])
      }
    }
    ZBE=as.matrix(temppairs[[1]])
    if(dim(ZBE)[2]==1){
      ZBE=t(ZBE)
    }
    ZXE=as.matrix(temppairs[[2]])
    if(dim(ZXE)[2]==1){
      ZXE=t(ZXE)
    }
    ZYE=as.matrix(temppairs[[3]])
    if(dim(ZYE)[2]==1){
      ZYE=t(ZYE)
    }
    ZNE=as.matrix(temppairs[[4]])
    if(dim(ZNE)[2]==1){
      ZNE=t(ZNE)
    }

    if(type=="BR"){
      func=function(x){

        -nllMSPexp(ZBE,ZXE,ZYE,ZNE,par=x,coord,n.a,sphere.dis)
      }
    }else if(type=="IBR"){
      func=function(x){

        -nllIMSPexp(ZBE,ZXE,ZYE,ZNE,par=x,coord,n.a,sphere.dis)
      }
    }
    s[i,]=grad(func,x=lik$par)

  }
  return(s)
}

#' @rdname Score_Est
#' @export
Score_Est_Smith=function(Z_Exp,u,coord,lik,type=c("Smith,InvSmith"),sphere.dis=F){

  N=dim(Z_Exp)[1]
  s=array(0,dim=c(N,1))
  for( i in 1:N){
    one.timepoint=t(as.matrix(Z_Exp[i,]))
    temppairs=makepairs(Z=one.timepoint,u=u)
    dist=rdist(coord+runif(dim(coord)[1],0,1))
    dist2=apply(temppairs[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
    unique.d=as.matrix(unique(dist2))
    n.a=rep(0, dim(unique.d)[1])
    if(length(n.a)!=0){
      for(j in 1:length(unique.d)){
        n.a[j]=sum(dist2==unique.d[j,1])
      }
    }
    ZBE=as.matrix(temppairs[[1]])
    if(dim(ZBE)[2]==1){
      ZBE=t(ZBE)
    }
    ZXE=as.matrix(temppairs[[2]])
    if(dim(ZXE)[2]==1){
      ZXE=t(ZXE)
    }
    ZYE=as.matrix(temppairs[[3]])
    if(dim(ZYE)[2]==1){
      ZYE=t(ZYE)
    }
    ZNE=as.matrix(temppairs[[4]])
    if(dim(ZNE)[2]==1){
      ZNE=t(ZNE)
    }

    if(type=="Smith"){
      func=function(x){

        -nllMSPexpSmith(ZBE,ZXE,ZYE,ZNE,par=x,coord,n.a,sphere.dis)
      }
    }else if(type=="InvSmith"){
      func=function(x){

        -nllIMSPexpSmith(ZBE,ZXE,ZYE,ZNE,par=x,coord,n.a,sphere.dis)
      }
    }
    s[i]=grad(func,x=lik$par)

  }
  return(s)
}
