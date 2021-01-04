library(sdfEXTREME)
##Simulation Study##
set.seed(1)
#This is repeated for 50 different seeds. We then take each of the CLAIC values for each of the model
#fits from each seed and calculate the proportion of seeds for which the particular deformation method
#provides the lowest.

##Choose one process -
#NSBR.true - Non-stationary Brown-Resnick and inverted Brown-Resnick processes
#maxmix.true - Max-mixture and inverted max-mixture
#Gaussmix.true - Gaussian-Mixture

#Note that if Gaussmix.true==TRUE, we only consider an asymptotically independent process
NSBR.true<-TRUE
maxmix.true<-FALSE
Gaussmix.true<-FALSE

if(NSBR.true == T | maxmix.true == T){
lambda<-2
centre<-c(0,0)
kappa<-0.8

n.grid<-8
sim.coords<-as.matrix(expand.grid(seq(-1,1,length=n.grid),seq(-1,1,length=n.grid)))

p<-dim(sim.coords)[1]
tau<-matrix(NA,nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    tau[i,j]<-gamsv(s1=sim.coords[i,],s2=c(0,0),lam=lambda,kap=kappa,centre=NULL)+gamsv(s1=sim.coords[j,],s2=c(0,0),lam=lambda,kap=kappa,centre=NULL)-gamsv(s1=sim.coords[i,],s2=sim.coords[j,],lam=lambda,kap=kappa,centre=centre)
  }
}

##Simulate Non-stationary Brown-Resnick process
##Simulates N=1000 replicates - Will take some time

N=1000
Sim.AD<-brnsims(reps=N,locs=sim.coords,kappa=kappa,lambda=lambda,centre=centre,tau=tau)


if(maxmix.true == T){
cor.mat<-Matern(rdist(sim.coords),range=1,smoothness = 1.2)
Y<-qgev(pnorm(rmvnorm(N,mean=rep(0,p),sigma=cor.mat)),1,1,1)
omega<-0.3

for(i in 1:p){
  Sim.AD[,i]<-apply(cbind(omega*Sim.AD[,i],(1-omega)*Y[,i]),1,max)
}
}

Sim.AI<-1/Sim.AD #Creates asymptotically independent version of process


}
if(Gaussmix.true==T){
  require("condMVNorm")
RadialDis<-function(s1,s2,centre){
  Euc.dist<-function(s1,s2=NULL){
    if(is.null(s2)){s2=c(0,0)}
    return( sqrt((s1[1]-s2[1])^2+(s1[2]-s2[2])^2))

  }
  fx=Euc.dist(s1,centre)*(s1-centre)
  fy=Euc.dist(s2,centre)*(s2-centre)
  h=Euc.dist(fx,fy)
}
n.grid<-9
sim.coords<-as.matrix(expand.grid(seq(-1,1,length=n.grid),seq(-1,1,length=n.grid)))
p<-dim(sim.coords)[1]
lambda<-2
centre<-c(0,0)
kappa<-0.8
#Make non-stationary correlation matrix
cor.NS<-matrix(NA,nrow=p,ncol=p)
for(i in 1:p){
  for( j in 1:p){
    cor.NS[i,j]<-Matern(RadialDis(sim.coords[i,],sim.coords[j,],centre),range=lambda,smoothness =kappa)
  }
}
#Make stationary correlation matrix

lambda.S<-2
kappa.S<-1
cor.S<-Matern(rdist(sim.coords),lambda.S,kappa.S)

N <- 1000# sample size
Sim.AI<-matrix(0,nrow=N,ncol=p)

for(i in 1:N){
  initial<-rnorm(1)
  centre.ind=ceiling(p/2)

  Sim.AI[i,centre.ind]<- initial #Iniate at centre

  if(initial > qnorm(0.9)){
    #Draw from non-stationary GP
    Sim.AI[i,-centre.ind]<- rcmvnorm(1,rep(0,p),sigma=cor.NS,dep=c(1:p)[-centre.ind],given.ind=centre.ind, X.given = initial)
  }else{
    #Draw from non-stationary GP
    Sim.AI[i,-centre.ind]<- rcmvnorm(1,rep(0,p),sigma=cor.S,dep=c(1:p)[-centre.ind],given.ind=centre.ind, X.given = initial)

  }
}
}

unif<-function(x){rank(x)/(length(x)+1)}
chi.emp<-function(U,V,z) sum((U>=z)&(V>=z))/sum(U>=z)

#Calculate pairwise chi and correlation estimates - Chi estimates above threshold q

q<-0.9

Sim.AI_U<-Sim.AI
#Transform to uniform margins
for(i in 1:dim(Sim.AI_U)[2]){

  Sim.AI_U[,i]<-unif(Sim.AI[,i])
}
#Transform to Normal margins and exponential margins (for likelihood)
Sim.AI_N<-qnorm(Sim.AI_U)
Sim.AI_Exp<-qexp(Sim.AI_U)

#Calculate pairwise empirical chi
emp.chi.AI<-matrix(rep(0,dim(Sim.AI_U)[2]^2),nrow=dim(Sim.AI_U)[2],ncol=dim(Sim.AI_U)[2])
for(i in 1:dim(Sim.AI_U)[2]){
  for(j in 1:i){

    emp.chi.AI[i,j]<-chi.emp(Sim.AI_U[,i],Sim.AI_U[,j],q)
  }

}
emp.chi.AI<-emp.chi.AI+t(emp.chi.AI)
diag(emp.chi.AI)<-diag(emp.chi.AI)/2

#Calculate pairwise empirical rho
emp.cor.AI<-matrix(rep(0,dim(Sim.AI_N)[2]^2),nrow=dim(Sim.AI_N)[2],ncol=dim(Sim.AI_N)[2])
for(i in 1:dim(Sim.AI_N)[2]){
  for(j in 1:i){

    emp.cor.AI[i,j]<-cor(Sim.AI_N[,i],Sim.AI_N[,j])
  }

}
emp.cor.AI<-emp.cor.AI+t(emp.cor.AI)
diag(emp.cor.AI)<-diag(emp.cor.AI)/2

#Asymptotically dependent process only considered if using NS Brown-Resnick or max-mixture
if(NSBR.true == T | maxmix.true == T){
Sim.AD_U<-Sim.AD
#Transform to uniform margins
for(i in 1:dim(Sim.AD_U)[2]){

  Sim.AD_U[,i]<-unif(Sim.AD[,i])
}
#Transform to Normal margins and exponential margins (for likelihood)
Sim.AD_N<-qnorm(Sim.AD_U)
Sim.AD_Exp<-qexp(Sim.AD_U)

#Calculate pairwise empirical chi
emp.chi.AD<-matrix(rep(0,dim(Sim.AD_U)[2]^2),nrow=dim(Sim.AD_U)[2],ncol=dim(Sim.AD_U)[2])
for(i in 1:dim(Sim.AD_U)[2]){
  for(j in 1:i){

    emp.chi.AD[i,j]<-chi.emp(Sim.AD_U[,i],Sim.AD_U[,j],q)
  }

}
emp.chi.AD<-emp.chi.AD+t(emp.chi.AD)
diag(emp.chi.AD)<-diag(emp.chi.AD)/2

#Calculate pairwise empirical rho
emp.cor.AD<-matrix(rep(0,dim(Sim.AD_N)[2]^2),nrow=dim(Sim.AD_N)[2],ncol=dim(Sim.AD_N)[2])
for(i in 1:dim(Sim.AD_N)[2]){
  for(j in 1:i){

    emp.cor.AD[i,j]<-cor(Sim.AD_N[,i],Sim.AD_N[,j])
  }

}
emp.cor.AD<-emp.cor.AD+t(emp.cor.AD)
diag(emp.cor.AD)<-diag(emp.cor.AD)/2
}

#Set an initial number of anchor points m.init

# Set indices for all anchor points used in final deformation. In the paper these were the same for each
# deformation and each random seed. These were chosen as they gave generally good deformations
# regardless of the seed. Note that we use more anchor points for the Gaussian mixture as there are more sampling locations

if(NSBR.true == T | maxmix.true == T){
  m.init<-10
  Full.m.ind<-c(37, 47, 55, 59, 29, 12, 18, 1, 32, 63,
                4, 57, 6, 17, 50, 7)
  }else if(Gaussmix.true == T ){
    m.init<-11
    Full.m.ind<-c(74, 19, 17, 8, 29, 27, 23, 39, 55, 51,
                  31, 10, 56, 71, 32, 53, 59, 14, 9, 52)
}


#Create deformations using heuristic. This will take a long time to run. We create 4 deformations for the sampled processes.
#sdf_chi: This uses a theoretical stationary chi function.
#sdf_chiq: This uses a theoretical stationary chi_q function.
#sdf_rho: This uses a theoretical stationary correlation function.
#sdf_smith: This uses the original Smith (1996) method, which uses the full stationary GP likelihood.



sdf.AI_chi<-sdf.heur.EXTR(m.init=m.init,Full.m.ind=Full.m.ind,
                          Gcoords=sim.coords,emp.chi=emp.chi.AI,type="CHI_BR")
sdf.AI_chiq<-sdf.heur.EXTR(m.init=m.init,Full.m.ind=Full.m.ind,
                           Gcoords=sim.coords,emp.chi=emp.chi.AI,type="CHI_q_IBR",q=q)
sdf.AI_rho<-sdf.heur.Cor(m.init=m.init,Full.m.ind=Full.m.ind,
                         Gcoords=sim.coords,emp.cor=emp.cor.AI,type="F-norm")
sdf.AI_smith<-sdf.heur.Cor(m.init=m.init,Full.m.ind=Full.m.ind,
                           Gcoords=sim.coords,emp.cor=emp.cor.AI,type="Smith",N=N)
if(NSBR.true == T | maxmix.true == T){
  sdf.AD_chi<-sdf.heur.EXTR(m.init=m.init,Full.m.ind=Full.m.ind,
                            Gcoords=sim.coords,emp.chi=emp.chi.AD,type="CHI_BR")
  sdf.AD_chiq<-sdf.heur.EXTR(m.init=m.init,Full.m.ind=Full.m.ind,
                             Gcoords=sim.coords,emp.chi=emp.chi.AD,type="CHI_q_IBR",q=q)
  sdf.AD_rho<-sdf.heur.Cor(m.init=m.init,Full.m.ind=Full.m.ind,
                           Gcoords=sim.coords,emp.cor=emp.cor.AD,type="F-norm")
  sdf.AD_smith<-sdf.heur.Cor(m.init=m.init,Full.m.ind=Full.m.ind,
                             Gcoords=sim.coords,emp.cor=emp.cor.AD,type="Smith",N=N)
}
#Create each set of Dcoords

Dcoords.AI_chi<-returnDcoord(sdf.AI_chi$par,sim.coords,sdf.AI_chi$m.ind,sphere.dis=F)
Dcoords.AI_chiq<-returnDcoord(sdf.AI_chiq$par,sim.coords,sdf.AI_chiq$m.ind,sphere.dis=F)
Dcoords.AI_rho<-returnDcoord(sdf.AI_rho$par,sim.coords,sdf.AI_rho$m.ind,sphere.dis=F)
Dcoords.AI_smith<-returnDcoord(sdf.AI_smith$par,sim.coords,sdf.AI_smith$m.ind,sphere.dis=F)
if(NSBR.true == T | maxmix.true == T){
  Dcoords.AD_chi<-returnDcoord(sdf.AD_chi$par,sim.coords,sdf.AD_chi$m.ind,sphere.dis=F)
  Dcoords.AD_chiq<-returnDcoord(sdf.AD_chiq$par,sim.coords,sdf.AD_chiq$m.ind,sphere.dis=F)
  Dcoords.AD_rho<-returnDcoord(sdf.AD_rho$par,sim.coords,sdf.AD_rho$m.ind,sphere.dis=F)
  Dcoords.AD_smith<-returnDcoord(sdf.AD_smith$par,sim.coords,sdf.AD_smith$m.ind,sphere.dis=F)
}
######Model Fitting######

#We fit 5 models for each process - One to the original Gcoords and one for each set of Dcoords

#Create pairs for likelihood

##AI fits
#u = threshold in Composite likelihood estimation
u<-quantile(Sim.AI_Exp,prob=q)

##Warning - Takes a while
Sim.AIpair<-makepairs(Sim.AI_Exp,u=q)

#Identify number of non-exceedances per pair
#Speeds up likelihood computation
dist<-rdist(sim.coords+runif(dim(sim.coords)[1]*2,0,1))
dist2<-apply(Sim.AIpair[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
unique.d<-as.matrix(unique(dist2))
n.a<-rep(0, dim(unique.d)[1])
for(i in 1:length(unique.d)){
  n.a[i]<-sum(dist2==unique.d[i,1])
}

#WARNING- pairwise likelihoods take some time to run.

#####IMSP fit on G-plane
# likAI.G<-optim(fn=nllIMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.AIpair[[1]]),
#             ZXE=as.matrix(Sim.AIpair[[2]]),
#             ZYE=as.matrix(Sim.AIpair[[3]]),ZNE=as.matrix(Sim.AIpair[[4]]),
#             n.a=n.a,sphere.dis=F,coord=sim.coords,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

####IMSP fit on chi D-plane
# likAI.Dchi<-optim(fn=nllIMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.AIpair[[1]]),
#             ZXE=as.matrix(Sim.AIpair[[2]]),
#             ZYE=as.matrix(Sim.AIpair[[3]]),ZNE=as.matrix(Sim.AIpair[[4]]),
#             n.a=n.a,sphere.dis=F,coord=Dcoords.AI_chi,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

####IMSP fit on chiq D-plane
# likAI.Dchiq<-optim(fn=nllIMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.AIpair[[1]]),
#             ZXE=as.matrix(Sim.AIpair[[2]]),
#             ZYE=as.matrix(Sim.AIpair[[3]]),ZNE=as.matrix(Sim.AIpair[[4]]),
#             n.a=n.a,sphere.dis=F,coord=Dcoords.AI_chiq,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

####IMSP fit on rho D-plane
# likAI.Drho<-optim(fn=nllIMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.AIpair[[1]]),
#             ZXE=as.matrix(Sim.AIpair[[2]]),
#             ZYE=as.matrix(Sim.AIpair[[3]]),ZNE=as.matrix(Sim.AIpair[[4]]),
#             n.a=n.a,sphere.dis=F,coord=Dcoords.AI_rho,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

####IMSP fit on Smith D-plane
# likAI.Dsmith<-optim(fn=nllIMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.AIpair[[1]]),
#             ZXE=as.matrix(Sim.AIpair[[2]]),
#             ZYE=as.matrix(Sim.AIpair[[3]]),ZNE=as.matrix(Sim.AIpair[[4]]),
#             n.a=n.a,sphere.dis=F,coord=Dcoords.AI_smith,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

if(NSBR.true == T | maxmix.true == T){
  ##AD fits
  #u = threshold in Composite likelihood estimation
  u<-quantile(Sim.AD_Exp,prob=q)

  ##Warning - Takes a while
  Sim.ADpair<-makepairs(Sim.AD_Exp,u=q)

  #Identify number of non-exceedances per pair
  #Speeds up likelihood computation
  dist<-rdist(sim.coords+runif(dim(sim.coords)[1]*2,0,1))
  dist2<-apply(Sim.ADpair[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
  unique.d<-as.matrix(unique(dist2))
  n.a<-rep(0, dim(unique.d)[1])
  for(i in 1:length(unique.d)){
    n.a[i]<-sum(dist2==unique.d[i,1])
  }

  #WARNING- pairwise likelihoods take some time to run.

  #####MSP fit on G-plane
  # likAD.G<-optim(fn=nllMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.ADpair[[1]]),
  #             ZXE=as.matrix(Sim.ADpair[[2]]),
  #             ZYE=as.matrix(Sim.ADpair[[3]]),ZNE=as.matrix(Sim.ADpair[[4]]),
  #             n.a=n.a,sphere.dis=F,coord=sim.coords,
  #             control=list(maxit=2000),
  #             method = "Nelder-Mead",hessian=T)

  ####MSP fit on chi D-plane
  # likAD.Dchi<-optim(fn=nllMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.ADpair[[1]]),
  #             ZXE=as.matrix(Sim.ADpair[[2]]),
  #             ZYE=as.matrix(Sim.ADpair[[3]]),ZNE=as.matrix(Sim.ADpair[[4]]),
  #             n.a=n.a,sphere.dis=F,coord=Dcoords.AD_chi,
  #             control=list(maxit=2000),
  #             method = "Nelder-Mead",hessian=T)

  ####MSP fit on chiq D-plane
  # likAD.Dchiq<-optim(fn=nllMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.ADpair[[1]]),
  #             ZXE=as.matrix(Sim.ADpair[[2]]),
  #             ZYE=as.matrix(Sim.ADpair[[3]]),ZNE=as.matrix(Sim.ADpair[[4]]),
  #             n.a=n.a,sphere.dis=F,coord=Dcoords.AD_chiq,
  #             control=list(maxit=2000),
  #             method = "Nelder-Mead",hessian=T)

  ####MSP fit on rho D-plane
  # likAD.Drho<-optim(fn=nllMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.ADpair[[1]]),
  #             ZXE=as.matrix(Sim.ADpair[[2]]),
  #             ZYE=as.matrix(Sim.ADpair[[3]]),ZNE=as.matrix(Sim.ADpair[[4]]),
  #             n.a=n.a,sphere.dis=F,coord=Dcoords.AD_rho,
  #             control=list(maxit=2000),
  #             method = "Nelder-Mead",hessian=T)

  ####MSP fit on Smith D-plane
  # likAD.Dsmith<-optim(fn=nllMSPexp,par=c(0.5,1),ZBE=as.matrix(Sim.ADpair[[1]]),
  #             ZXE=as.matrix(Sim.ADpair[[2]]),
  #             ZYE=as.matrix(Sim.ADpair[[3]]),ZNE=as.matrix(Sim.ADpair[[4]]),
  #             n.a=n.a,sphere.dis=F,coord=Dcoords.AD_smith,
  #             control=list(maxit=2000),
  #             method = "Nelder-Mead",hessian=T)

}

#CLAIC estimation

#Estimate score at each time point - Takes a while and must be done for each individual likelihood
u<-quantile(Sim.AI_Exp,prob=q)
s.AI_G= Score_Est(Sim.AI_Exp,u,coord=sim.coords,lik=likAI.G,type="IBR",sphere.dis=F)
s.AI_Dchi= Score_Est(Sim.AI_Exp,u,coord=Dcoords.AI_chi,lik=likAI.Dchi,type="IBR",sphere.dis=F)
s.AI_Dchiq= Score_Est(Sim.AI_Exp,u,coord=Dcoords.AI_chiq,lik=likAI.Dchiq,type="IBR",sphere.dis=F)
s.AI_Drho= Score_Est(Sim.AI_Exp,u,coord=Dcoords.AI_rho,lik=likAI.Drho,type="IBR",sphere.dis=F)
s.AI_Dsmith= Score_Est(Sim.AI_Exp,u,coord=Dcoords.AI_smith,lik=likAI.Dsmith,type="IBR",sphere.dis=F)
if(NSBR.true == T | maxmix.true == T){
  u<-quantile(Sim.AD_Exp,prob=q)
  s.AD_G= Score_Est(Sim.AD_Exp,u,coord=sim.coords,lik=likAD.G,type="BR",sphere.dis=F)
  s.AD_Dchi= Score_Est(Sim.AD_Exp,u,coord=Dcoords.AD_chi,lik=likAD.Dchi,type="BR",sphere.dis=F)
  s.AD_Dchiq= Score_Est(Sim.AD_Exp,u,coord=Dcoords.AD_chiq,lik=likAD.Dchiq,type="BR",sphere.dis=F)
  s.AD_Drho= Score_Est(Sim.AD_Exp,u,coord=Dcoords.AD_rho,lik=likAD.Drho,type="BR",sphere.dis=F)
  s.AD_Dsmith= Score_Est(Sim.AD_Exp,u,coord=Dcoords.AD_smith,lik=likAD.Dsmith,type="BR",sphere.dis=F)
}

#Estimate variance of score - We use block sizes of 1 as the simulations are all temporally independent

varS.AI_G=var(s.AI_G)
varS.AI_Dchi=var(s.AI_Dchi)
varS.AI_Dchiq=var(s.AI_Dchiq)
varS.AI_Drho=var(s.AI_Drho)
varS.AI_Dsmith=var(s.AI_Dsmith)
if(NSBR.true == T | maxmix.true == T){
  varS.AD_G=var(s.AD_G)
  varS.AD_Dchi=var(s.AD_Dchi)
  varS.AD_Dchiq=var(s.AD_Dchiq)
  varS.AD_Drho=var(s.AD_Drho)
  varS.AD_Dsmith=var(s.AD_Dsmith)
}
#Estimate CLAIC values
CLAIC.AI_G=2*likAI.G$value+2*sum(diag(varS.AI_G%*%solve(likAI.G$hessian)))
CLAIC.AI_Dchi=2*likAI.Dchi$value+2*sum(diag(varS.AI_Dchi%*%solve(likAI.Dchi$hessian)))
CLAIC.AI_Dchiq=2*likAI.Dchiq$value+2*sum(diag(varS.AI_Dchiq%*%solve(likAI.Dchiq$hessian)))
CLAIC.AI_Drho=2*likAI.Drho$value+2*sum(diag(varS.AI_Drho%*%solve(likAI.Drho$hessian)))
CLAIC.AI_Dsmith=2*likAI.Dsmith$value+2*sum(diag(varS.AI_Dsmith%*%solve(likAI.Dsmith$hessian)))
if(NSBR.true == T | maxmix.true == T){
  CLAIC.AD_G=2*likAD.G$value+2*sum(diag(varS.AD_G%*%solve(likAD.G$hessian)))
  CLAIC.AD_Dchi=2*likAD.Dchi$value+2*sum(diag(varS.AD_Dchi%*%solve(likAD.Dchi$hessian)))
  CLAIC.AD_Dchiq=2*likAD.Dchiq$value+2*sum(diag(varS.AD_Dchiq%*%solve(likAD.Dchiq$hessian)))
  CLAIC.AD_Drho=2*likAD.Drho$value+2*sum(diag(varS.AD_Drho%*%solve(likAD.Drho$hessian)))
  CLAIC.AD_Dsmith=2*likAD.Dsmith$value+2*sum(diag(varS.AD_Dsmith%*%solve(likAD.Dsmith$hessian)))
}
