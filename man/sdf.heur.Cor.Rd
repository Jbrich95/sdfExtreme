% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdfheurCOR.R
\name{sdf.heur.Cor}
\alias{sdf.heur.Cor}
\title{Heuristic for creating bijective deformations (Correlation)}
\usage{
sdf.heur.Cor(
  m.init,
  Full.m.ind,
  Gcoords,
  emp.cor,
  type = c("F-norm", "Smith"),
  n = 0,
  par = NULL,
  sphere.dis = F
)
}
\arguments{
\item{m.init}{Number of inital anchor points. Must have \code{m.init < length(Full.m.ind)}.}

\item{Full.m.ind}{Full vector of indices for anchor points in \code{Gcoords}.}

\item{Gcoords}{A \code{d} by 2 matrix of G-plane coordinates.}

\item{emp.cor}{A \code{d} by \code{d} matrix of pairwise empirical correlation values.}

\item{type}{\code{"F-norm"} for Frobenius norm method using theoretical \eqn{\rho(h^{*}_{ij})} from a Matérn correlation model. \code{"Smith"} for original Smith (1996) method, also using Matérn correlation model.}

\item{n}{Number of data points used to estimate \code{emp.cor}. Only neccesary for \code{type=="Smith"}.}

\item{par}{Initial parameters for first deformation. If not stated, initial parameters are given.}

\item{sphere.dis}{Is Spherical distance or Euclidean distance used?}
}
\value{
List with three elements:\describe{
\item{par}{Spline parameter values.}
\item{value}{Objective value from final optimisation.}
\item{m.ind}{Vector of indices for full set of anchor points.}
}
}
\description{
This algorithm uses \code{FindSplineParFNCOR} to create a deformation with \code{m.init} inital anchor points.
The initial anchor points are the last \code{m.init} entries of \code{Full.m.ind}. After creating this deformation,
the spline values are used as initial parameters in creating a deformation with \code{m.init+1} anchor points.
This iterative procedure repeats until a deformation with all \code{Full.m.ind} anchor points is created.
}
\examples{
data("Aus_Heat")

Z<-Aus_Heat$Temp.
Gcoords<-Aus_Heat$coords

Z_U<-Z

unif<-function(x){rank(x)/(length(x)+1)}

#Transform to uniform margins
for(i in 1:dim(Z_U)[2]){

  Z_U[,i]<-unif(Z[,i])
}


#Transform to Gaussian margins

Z_N<-qnorm(Z_U)

#Calculate pairwise empirical correlation
emp.cor<-matrix(rep(0,dim(Z_N)[2]^2),nrow=dim(Z_N)[2],ncol=dim(Z_N)[2])
for(i in 1:dim(Z_N)[2]){
  for(j in 1:i){

    emp.cor[i,j]<-cor(Z_N[,i],Z_N[,j])
  }

}
emp.cor<-emp.cor+t(emp.cor)
diag(emp.cor)<-diag(emp.cor)/2


m.init<-6
Full.m.ind<-sample(1:dim(Gcoords)[1],m.init+2)

#Transform to D-plane using Smith (1996) method

##WARNING: This may take a while to run.
sdf<-sdf.heur.Cor(m.init,Full.m.ind,Gcoords,emp.cor,type="Smith",n=dim(Z)[1],sphere.dis=TRUE)

#Plot Dcoords

Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
plot(Dcoords,main="D-plane",ylab="",xlab="")
}
