#' A heuristic for deciding if initial anchor points have enough spread.
#'
#' This heuristic ensures that chosen anchor points are not along a straight line and,
#' if \code{length(m.ind)<5}, that the maximum distance between the initial anchor points in the
#' G-plane is sufficiently large i.e. over \eqn{60\%} of the maximum pairwise distance of all sampling locations in the G-plane.
#'
#'
#' @param m.ind A vector of indices denoting the anchor points in \code{Gcoords}.
#' @param Gcoords A \code{d} by 2 matrix of G-plane sampling locations.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?
#' @return \code{1} if anchor points are rejected, \code{0} otherwise.
#' @examples
#' data("Aus_Heat")
#'
#'Gcoords<-Aus_Heat$coords
#'# Set number of anchor points
#'m<-10
#'# Sample anchor points
#'m.ind<-sample(1:dim(Gcoords)[1],m,replace=FALSE)
#'
#'reject<-m.reject(m.ind,Gcoords,sphere.dis=TRUE)
#'print(reject)


m.reject=function(m.ind,Gcoords,sphere.dis=FALSE){
  m.coord=Gcoords[m.ind,]

  if(length(m.ind)<=5){
    if(sphere.dis==FALSE){
      max.dist=0.6*max(rdist(Gcoords))
      max.dis.comp=max(rdist(m.coord))
    }else{
      max.dist=0.6*max(rdist.earth(Gcoords,miles=FALSE))

      max.dis.comp=max(rdist.earth(m.coord,miles=FALSE))

    }
    angle = function(x,y){
      prod=sum(x*y)/(sqrt(sum(x*x))*sqrt(sum(y*y)))
      prod=round(prod,3)
      theta=acos(prod)
      round(as.numeric(theta),2)
    }
    min.angle=2*pi+1
    for(i in 3:(length(m.ind))){
      for(j in 2:(i-1)){
        for(k in 1:(j-1)){
          angle.compare=angle(m.coord[i,]-m.coord[j,],m.coord[i,]-m.coord[k,])
          if(angle.compare<=min.angle){
            min.angle=angle.compare
          }

        }
      }
    }
    if(min.angle>0 && max.dis.comp>=max.dist*0.6){
      reject=0

    }else{
      reject=1
    }
  }else{
    reject=0
  }
  m=length(m.ind)
  dmat<-matrix(c(1,1,1,m.coord[c(m-2,m-1,m),1],m.coord[c(m-2,m-1,m),2]),3,3,byrow = TRUE)
  test=try(solve(dmat),silent=TRUE)
  if(length(test)==1){
    reject=1
  }
  if(reject==1){print("Anchor Points Rejected")}else if(reject==0){print("Anchor Points Accepted")}
  return(reject)
}
