% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sdfheurEXTR.R
\name{sdf.heur.EXTR}
\alias{sdf.heur.EXTR}
\title{Heuristic for creating bijective deformations (Extremal)}
\usage{
sdf.heur.EXTR(
  m.init,
  Full.m.ind,
  Gcoords,
  emp.chi,
  type = c("CHI_BR", "CHI_q_IBR"),
  q = 0,
  par = NULL,
  sphere.dis = F
)
}
\arguments{
\item{m.init}{Number of inital anchor points. Must have \code{m.init < length(Full.m.ind)}.}

\item{Full.m.ind}{Full vector of indices for anchor points in \code{Gcoords}.}

\item{Gcoords}{A \code{d} by 2 matrix of G-plane coordinates.}

\item{emp.chi}{A \code{d} by \code{d} matrix of pairwise empirical chi values.}

\item{type}{\code{"CHI_BR"} for theoretical \eqn{\chi(h^{*}_{ij})} from a Brown-Resnick model. \code{"CHI_q_IBR"} for theoretical \eqn{\chi_q(h^{*}_{ij})} from an inverted Brown-Resnick model.}

\item{q}{Threhold for \eqn{\chi_q(h^*_ij)}. Only needed if \code{type=="CHI_q_IBR"}.}

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
This algorithm uses \code{FindSplineParFNEXTR} to create a deformation with \code{m.init} inital anchor points.
The initial anchor points are the last \code{m.init} entries of \code{Full.m.ind}. After creating this deformation,
the spline values are used as initial parameters in creating a deformation with \code{m.init + 1} anchor points.
This iterative procedure repeats until a deformation with all \code{Full.m.ind} anchor points is created.
}
\examples{
# Note that the given code will not produce the Australian Summer Temperature deformation
# used in the paper. See Data_Example.R and help(Aus_Heat_Output) for the deformation
# used in the paper. This deformation will use much fewer anchor points as the full
# deformation takes a long time to run.

data("Aus_Heat")

Z<-Aus_Heat$Temp.
Gcoords<-Aus_Heat$coords

Z_U<-Z

unif<-function(x){rank(x)/(length(x)+1)}

#Transform to uniform margins
for(i in 1:dim(Z_U)[2]){

  Z_U[,i]<-unif(Z[,i])
}

chi.emp<-function(U,V,z) sum((U>=z)&(V>=z))/sum(U>=z)
q<-0.98

#Calculate pairwise empirical chi
emp.chi<-matrix(rep(0,dim(Z_U)[2]^2),nrow=dim(Z_U)[2],ncol=dim(Z_U)[2])
for(i in 1:dim(Z_U)[2]){
  for(j in 1:i){

    emp.chi[i,j]<-chi.emp(Z_U[,i],Z_U[,j],q)
  }

}
emp.chi<-emp.chi+t(emp.chi)
diag(emp.chi)<-diag(emp.chi)/2

m.init<-6
Full.m.ind<-sample(1:dim(Gcoords)[1],m.init+2)

#Transform to D-plane using Brown-Resnick theoretical chi function

##WARNING: This may take a while to run. 
sdf<-sdf.heur.EXTR(m.init,Full.m.ind,Gcoords,emp.chi,type="CHI_BR",
                     par=c(.05,.05,0.2,1.6,rep(0,2*m.init-6)),sphere.dis=TRUE)

#Plot Dcoords

Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=TRUE)
plot(Dcoords,main="D-plane",ylab="",xlab="")
}
