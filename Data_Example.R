library(sdfEXTREME)
##Code has been written for Australian Sumemr temperatures case study, but is easily adaptable for precipitation case studies,
data(Aus_Heat)
##data(Snow_Precip)
##data(High_Precip)


Z<-Aus_Heat$Temp.
Gcoords<-Aus_Heat$coords

Z_U<-Z

unif<-function(x){rank(x)/(length(x)+1)}

#Transform to uniform margins
for(i in 1:dim(Z_U)[2]){

  Z_U[,i]<-unif(Z[,i])
}

chi.emp<-function(U,V,z) sum((U>=z)&(V>=z))/sum(U>=z)
q<-0.98 ##Australian Heatwaves. For precipitation, use q <- 0.95

#Calculate pairwise empirical chi
emp.chi<-matrix(rep(0,dim(Z_U)[2]^2),nrow=dim(Z_U)[2],ncol=dim(Z_U)[2])
for(i in 1:dim(Z_U)[2]){
  for(j in 1:i){

    emp.chi[i,j]<-chi.emp(Z_U[,i],Z_U[,j],q)
  }

}
emp.chi<-emp.chi+t(emp.chi)
diag(emp.chi)<-diag(emp.chi)/2

#These are the anchor points, correctly order, that were used to create
#Australian summer temperature deformation in the paper. Deformations are ready made for precipitation data, see help(Snow_Precip_Output) and help(High_Precip_Output).
m.init<-10
Full.m.ind<-c(32, 61, 7, 64, 72, 66, 39, 55, 9, 11,
              43, 46, 52, 26, 3, 1, 42, 23)

#Transform to D-plane using Brown-Resnick theoretical chi function

##WARNING: This may take a while to run. Use help(Aus_Heat_Output) for sdf output
#sdf<-sdf.heur.EXTR(m.init,Full.m.ind,Gcoords,emp.chi,type="CHI_BR",par=c(.05,.05,0.2,1.6,rep(0,2*m.init-6)),sphere.dis=T)
data("Aus_Heat_Output")
sdf<-Aus_Heat_Output$sdf

#Plot Dcoords

Dcoords<-returnDcoord(sdf$par,Gcoords,sdf$m.ind,sphere.dis=T)
plot(Dcoords)
# For details of plotting grid-lines, see help(returnDGridcoord)

#plot Chi vs. distance

par(mfrow=c(2,1))
plot(rdist.earth(Gcoords,miles=F),emp.chi,ylab="",main="G-Plane",xlab=expression(h[ij]))
mtext(side=2, expression(hat(chi)(h[ij])),cex = 1, line=2.2)
plot(rdist.earth(Dcoords,miles=F),emp.chi,ylab="",main="D-Plane",xlab=expression(h[ij]))
mtext(side=2, expression(hat(chi)(h[ij])),cex = 1, line=2.2)
#Theoretical Chi for deformation
v_H=(seq(0,max(rdist.earth(Dcoords,miles=F)),length=2000))^sdf$par[4]
t_H=2*pnorm(sqrt(v_H/2))
chi=2-t_H
points(y=chi,x=seq(0,max(rdist.earth(Dcoords,miles=F)),length=2000),type="l",col="red")

######Model Fitting######

#Transform to exponential margins
Z_Exp<-qexp(Z_U)
#Create pairs
#u = threshold in Composite likelihood estimation
u<-quantile(Z_Exp,prob=q)
##Warnings - Takes a while
Zpair<-makepairs(Z_Exp,u=u)

#Identify number of non-exceedances per pair
#Speeds up likelihood computation
dist<-rdist(Gcoords+runif(dim(Gcoords)[1]*2,0,1))
dist2<-apply(Zpair[[4]][,3:4],1,function(x){dist[x[1],x[2]]})
unique.d<-as.matrix(unique(dist2))
n.a<-rep(0, dim(unique.d)[1])
for(i in 1:length(unique.d)){
  n.a[i]<-sum(dist2==unique.d[i,1])
}

#WARNING- pairwise likelihoods take some time to run. Use data(Aus_Heat_Output).
data("Aus_Heat_Output")
likG.MSP<-Aus_Heat_Output$likG.MSP
likD.MSP<-Aus_Heat_Output$likD.MSP

#####MSP fit on G-plane (code for IMSP fit also available)
# likG.MSP<-optim(fn=nllMSPexp,par=c(1.5,500),ZBE=as.matrix(Zpair[[1]]),
#             ZXE=as.matrix(Zpair[[2]]),
#             ZYE=as.matrix(Zpair[[3]]),ZNE=as.matrix(Zpair[[4]]),
#             n.a=n.a,sphere.dis=T,coord=Gcoords,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

####MSP fit on D-plane
# likD.MSP<-optim(fn=nllMSPexp,par=c(1.5,1),ZBE=as.matrix(Zpair[[1]]),
#             ZXE=as.matrix(Zpair[[2]]),
#             ZYE=as.matrix(Zpair[[3]]),ZNE=as.matrix(Zpair[[4]]),
#             n.a=n.a,sphere.dis=T,coord=Dcoords,
#             control=list(maxit=2000),
#             method = "Nelder-Mead",hessian=T)

##See help("nllIMSPexpSmith") for IMSP fits
#CLAIC estimation

#Estimate score at each time point - Take a while
s_G<- Score_Est(Z_Exp,u,coord=Gcoords,lik=likG.MSP,type="BR",sphere.dis=T)
s_D<- Score_Est(Z_Exp,u,coord=Dcoords,lik=likD.MSP,type="BR",sphere.dis=T)

#Estimate variance of score
block.sizes<-c(91,90) #Set block.size # Here we take 90 and 91, corresponding to a regular season
#and a season with a leap year. For Precipitation data, all block sizes are of length 180 ( 90*2 half days)
years<-1957:2014
k<-length(years)

temp<-matrix(0,nrow=k,ncol=2)
temp2<-matrix(0,nrow=k,ncol=2)
int<-0
for(l in 1:k){
  if(years[l]%%4==0) block.size<-block.sizes[1] else block.size<-block.sizes[2]

  int <- int + block.size
  temp[l,]<-colSums(s_G[(int-block.size+1):(int),])
  temp2[l,]<-colSums(s_D[(int-block.size+1):(int),])

}

varS_G<-var(temp)
varS_D<-var(temp2)


#Estimate CLAIC
CLAIC_G<-2*likG.MSP$value+2*sum(diag(varS_G%*%solve(likG.MSP$hessian)))
CLAIC_D<-2*likD.MSP$value+2*sum(diag(varS_D%*%solve(likD.MSP$hessian)))


#Diagnostics

#Stationary Bootstrap - Here we give for the max-stable fit only, but chi3 can be replaced by chi3q for IMSP fit

#Define set of possible indices - Here we use North-West/South-East transect only. Index.set will need changing for precipitation data.
Index.set<-vector("list",10)
Index.set[[1]]<-1:5
Index.set[[2]]<-6:10
Index.set[[3]]<-11:17
Index.set[[4]]<-18:24
Index.set[[5]]<-25:31
Index.set[[6]]<-32:40
Index.set[[7]]<-41:48
Index.set[[8]]<-49:58
Index.set[[9]]<-59:66
Index.set[[10]]<-67:72
block.mean<-14 #Mean block size for random block choice - Here a fortnight
B<-1000 # Number of bootstrap samples for each instance

M<-10
emp.chi3<-chi3.UQ<-chi3.LQ<-the.chi3<-the.chi3G<-rep(0,M) # M - number of triples considered
q<-0.98 #0.95 for precipitation
#Calculate theoretical chi(h_ijk) and bootstrap estimates + 5% and 95% quantiles
for(it in 1:M){
  ind.sub<-sample(1:length(Index.set),1) # Pick possible subset
  inds<-Index.set[[ind.sub]]
  start.ind<-sample(1:(length(inds)-2),1)
  ind.triple<-inds[start.ind:(start.ind+2)] # Will not pick indices that are not in that list component/ will not wrap around

  emp.chi3[it]<-chi3.emp(Z_U[,ind.triple[1]],Z_U[,ind.triple[2]],Z_U[,ind.triple[3]],q)

  chi3.emp.boot<-rep(0,B)
  for(i in 1:B){
    b<-0
    while(b==0){

      b<-rgeom(1,1/(block.mean+1))

    }
    boot<-stat.boot(Z_U[,ind.triple],b)
    chi3.emp.boot[i]<-chi3.emp(boot[,1],boot[,2],boot[,3],q)
  }

  chi3.UQ[it]<-quantile(chi3.emp.boot,prob=0.975)
  chi3.LQ[it]<-quantile(chi3.emp.boot,prob=0.025)

  v_H<-(rdist(Dcoords[ind.triple,])/likD.MSP$par[2])^likD.MSP$par[1]

  the.chi3[it]<-chi3(v_H)
}

library(plotrix)

plotCI(1:M,emp.chi3,ui=chi3.UQ,li=chi3.LQ,ylab="",xlab="Index",main="Model Fit Diagnostic")
mtext(side=2, expression(hat(chi)(h[ijk])),cex = 1, line=2.2)

points(the.chi3,col="red")

#Conditional expectation (Heffernan & Tawn, 2004)

Z_LP<-qlaplace(Z_U) # Change to Laplace margins

p<-dim(Gcoords)[1]
ConExp<-matrix(0,nrow=p,ncol=p)

u<-quantile(Z_LP,q) #Quantile for estimating conditional expectation

#Calculate conditional expectation measure for each pair
for(i in 1:p){
  for(j in 1:i){
    Exceedances<-cbind(Z_LP[,i],Z_LP[,j])[which(Z_LP[,i]>=u),]
    opt<-optim(nllHT,X=Exceedances[,1],Y=Exceedances[,2],par=c(0.3,0.8))
    #print(opt)
    alpha<-opt$par[1]
    beta<-opt$par[2]
    mu<-mean((Exceedances[,2]-alpha*Exceedances[,1])/Exceedances[,1]^beta)
    ConExp[i,j]<-alpha*u+u^beta*mu
  }
  print(i)
}
ConExp<-ConExp+t(ConExp) #Symmetry assumed
diag(ConExp)<-diag(ConExp)/2

#plot on both planes
par(mfrow=c(2,1))
plot(rdist.earth(Gcoords,miles=F),ConExp,ylab="Conditional Expectation",xlab="Distance",main="G-Plane")
plot(rdist.earth(Dcoords,miles=F),ConExp,ylab="Conditional Expectation",xlab="Distance",main="D-Plane")

#overlay with mean distances normalised
par(mfrow=c(1,1))
plot(rdist.earth(Gcoords,miles=F),ConExp,xlab="Pairwise distance",ylab="Conditional Expectation",main="Deformation Diagnostic")
Dcoordnorm<-Dcoords*mean(rdist.earth(Gcoords,miles=F))/mean(rdist.earth(Dcoords,miles=F))
points(rdist.earth(Dcoordnorm,miles=F),ConExp,col="red")

