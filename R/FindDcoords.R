#' Return D-plane coordinates
#'
#' Returns the D-plane coordinates given spline parameter values calculated using \code{FindSplineParFNEXTR} or \code{FindSplineParFNCOR}.
#'
#' @param par Spline parameter values calculated using either \code{FindSplineParFNEXTR} or \code{FindSplineParFNCOR} or \code{sdf.heur.Cor} or \code{sdf.heur.EXTR}.
#' @param Gcoords A \code{d} by 2 matrix of G-plane sampling locations.
#' @param gridcoord A \code{K} by 2 matrix of any coordinates in the G-plane.
#' @param whichm A vector of length \code{m < d} giving the indices of the anchor points in \code{Gcoords}. Must be the same as used to find \code{par}.
#' @param  sphere.dis Is Spherical distance or Euclidean distance used?

#' @return \code{returnDcoord} returns a \code{d} by 2 matrix of coordinates. \code{returnDGridcoord} returns a \code{K} by 2 matrix of coordinates.

#' @rdname ReturnDcoords
#' @export
returnDcoord<-function(par,Gcoords,whichm,sphere.dis=F){
  m<-length(whichm)

  b1<-par[1]
  b2<-par[2]
  rho<-par[3]
  nu<-par[4]

  coordm<-Gcoords[whichm,]
  delta1<-delta2<-rep(0,3)
  if(m>3){
    delta1<-par[5:(length(par)-m+3)]
    delta2<-par[(length(par)-m+4):length(par)]


    dmat<-matrix(c(1,1,1,coordm[c(m-2,m-1,m),1],coordm[c(m-2,m-1,m),2]),3,3,byrow = T)
    rhsd1<- -c(sum(delta1),sum(delta1*coordm[-c(m-2,m-1,m),1]), sum(delta1*coordm[-c(m-2,m-1,m),2]))
    rhsd2<- -c(sum(delta2),sum(delta2*coordm[-c(m-2,m-1,m),1]), sum(delta2*coordm[-c(m-2,m-1,m),2]))

    delta1final3<- solve(dmat)%*%rhsd1
    delta2final3<- solve(dmat)%*%rhsd2

    delta1<-c(delta1,delta1final3)
    delta2<-c(delta2,delta2final3)
  }

  etafunc<-function(Gcoordvec,s1,s2,sphere.dis=FALSE)
  {
    if(sphere.dis==FALSE){
      r<-sqrt((s1-Gcoordvec[1])^2+(s2-Gcoordvec[2])^2)
    }else{
      r<-rdist.earth(rbind(c(s1,s2),Gcoordvec),miles=F)[1,2]

    }
    if(r==0){return(0)}
    return(r^2 * log(r))
  }

  f1<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b1^2*s1 + rho*b1*b2*s2 + sum(delta1*eta)
  }
  f2<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b2^2*s2 + rho*b1*b2*s1 + sum(delta2*eta)
  }
  f<-function(s1s2)
  {
    s1<-s1s2[1]
    s2<-s1s2[2]
    return(c(f1(s1,s2),f2(s1,s2)))
  }

  Dcoord<-t(apply(Gcoords,1,f))
  return(Dcoord)
}
#' @rdname ReturnDcoords
#' @export
returnDGridcoord<-function(par,gridcoord,Gcoords,whichm,sphere.dis=F){
  m<-length(whichm)

  b1<-par[1]
  b2<-par[2]
  rho<-par[3]
  nu<-par[4]

  coordm<-Gcoords[whichm,]
  if(m>3){
    delta1<-par[5:(length(par)-m+3)]
    delta2<-par[(length(par)-m+4):length(par)]


    dmat<-matrix(c(1,1,1,coordm[c(m-2,m-1,m),1],coordm[c(m-2,m-1,m),2]),3,3,byrow = T)
    rhsd1<- -c(sum(delta1),sum(delta1*coordm[-c(m-2,m-1,m),1]), sum(delta1*coordm[-c(m-2,m-1,m),2]))
    rhsd2<- -c(sum(delta2),sum(delta2*coordm[-c(m-2,m-1,m),1]), sum(delta2*coordm[-c(m-2,m-1,m),2]))

    delta1final3<- solve(dmat)%*%rhsd1
    delta2final3<- solve(dmat)%*%rhsd2

    delta1<-c(delta1,delta1final3)
    delta2<-c(delta2,delta2final3)
  }

  etafunc<-function(Gcoordvec,s1,s2,sphere.dis=FALSE)
  {
    if(sphere.dis==FALSE){
      r<-sqrt((s1-Gcoordvec[1])^2+(s2-Gcoordvec[2])^2)
    }else{
      r<-rdist.earth(rbind(c(s1,s2),Gcoordvec),miles=F)[1,2]

    }
    if(r==0){return(0)}
    return(r^2 * log(r))
  }

  f1<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b1^2*s1 + rho*b1*b2*s2 + sum(delta1*eta)
  }
  f2<-function(s1,s2)
  {
    eta<-apply(coordm,1,etafunc,s1=s1,s2=s2,sphere.dis=sphere.dis)
    b2^2*s2 + rho*b1*b2*s1 + sum(delta2*eta)
  }
  f<-function(s1s2)
  {
    s1<-s1s2[1]
    s2<-s1s2[2]
    return(c(f1(s1,s2),f2(s1,s2)))
  }

  Dcoord<-t(apply(gridcoord,1,f))
  return(Dcoord)
}
