library(foreach)
library(GORCure)
source("sim_cure.R")


f1<-function(x){sqrt(2)*sin(2*pi*x)}
f2<-function(x){sqrt(2)*cos(2*pi*x)}

simu<-function(true_theta_0,true_theta_xi,true_beta_0,true_beta_xi){

timegrid = seq(0,1,by=0.01)

ssize=200
a1=rnorm(ssize,0,2)
a2=rnorm(ssize,0,1)

observed = lapply(1:ssize,function(i){
mu = 0
pc1add=sapply(timegrid,f1)*a1[i]
pc2add= sapply(timegrid,f2)*a2[i]
err=rnorm(1,mean=0,sd=0.01)
as.numeric(mu+pc1add+pc2add+err)
})

timepoints = lapply(1:ssize,function(i){
timegrid
})

library(fdapace)
library(fda)
res_pace<- FPCA(observed, timepoints,list(dataType='Dense',error=TRUE, kernel='epan', verbose=TRUE,methodBwCov="GCV",methodBwMu="GCV",methodSelectK = 2))


spline_basis=create.bspline.basis(rangeval=c(0,1),nbasis=10,norder=4)

pacemu = Data2fd(argvals = res_pace$workGrid,y=res_pace$mu,spline_basis)
pacefpc1 = Data2fd(argvals = res_pace$workGrid,y=res_pace$phi[,1],spline_basis)
pacefpc2 = Data2fd(argvals = res_pace$workGrid,y=res_pace$phi[,2],spline_basis)

imse = integrate(function(x){(eval.fd(x, pacefpc1) - f1(x))^2}, lower = 0, upper = 1)

scoresmat = cbind(a1,a2)

  N=ssize
  intercept <- 1
  z1 <- stats::rbinom(N,1,0.5)
  z2 <- stats::rnorm(N)
  cov_theta_0 <- data.frame(intercept, z1, z2)
  cov_beta_0 <- data.frame(z1, z2)
  
  
  cov_theta_xi = as.data.frame(scoresmat)
  cov_beta_xi= as.data.frame(scoresmat)
  
  

  theta_0 = true_theta_0
  theta_xi=  true_theta_xi
  beta_0 =   true_beta_0
  beta_xi= true_beta_xi
  
  truepar = c(theta_0,theta_xi,beta_0,beta_xi)
  
  #cureset=sim_cure(ssize, cov_theta_0,cov_beta_0,cov_theta_xi,cov_beta_xi,theta_0 ,theta_xi,beta_0 ,beta_xi,A = 5, B = 15) 
  cureset=sim_cure_GOR(ssize, cov_theta_0,cov_beta_0,cov_theta_xi,cov_beta_xi,theta_0 ,theta_xi,beta_0 ,beta_xi,A = 50, B = 15) 
  
  
  
  scoreest=res_pace$xiEst[,1:2]

  signmat = sign(scoresmat)==sign(scoreest)
  ss=apply(signmat,2,sum)
  if(ss[1]<100){
  scoreest[,1]=scoreest[,1]*(-1)
  }
   if(ss[2]<100){
  scoreest[,2]=scoreest[,2]*(-1)
  }
  
  
cureset[,c("theta_xi1", "theta_xi2")]=scoreest
cureset[,c("beta_xi1", "beta_xi2")]=scoreest
  
  
fulldata=as.data.frame(cureset)
fulldata$R[which(fulldata$R==Inf)]<-NA

fit<-GORMC(survfun=Surv(L,R)~theta_z1+theta_z2+theta_xi1+theta_xi2,curefun=~beta_z1+beta_z2+beta_xi1+beta_xi2,data=fulldata,r=0)
#summary(fit)

xx=c(fit$ParEst$Eta,fit$ParEst$Beta)
yy=sqrt(diag(fit$ParVcov$vcov.bg))
return(list(xx,yy))

}





						

true_theta_0 = c(0.05,0,-1)
true_theta_xi=c(-1,0)
true_beta_0 = c(0.2,0.3)
true_beta_xi=c(0.0,2)
truepar = c(true_theta_0,true_theta_xi,true_beta_0,true_beta_xi)

numofsimu=500
olist=list()
for(h in 1:numofsimu){			
	olist[[h]]=try(simu(true_theta_0,true_theta_xi,true_beta_0,true_beta_xi))

	while(( 'try-error' %in% class(olist[[h]] ) || is.nan(sum(olist[[h]][[2]])) )){
	olist[[h]] = try(simu(true_theta_0,true_theta_xi,true_beta_0,true_beta_xi))
	}
}



			
		
est = lapply(olist,function(o){
est = o[[1]]
std = o[[2]]
return(est)
})

std = lapply(olist,function(o){
est = o[[1]]
std = o[[2]]
return(std)
})


stdmat = Reduce(rbind,std)
stdmat2=stdmat[sapply(std,function(x)!is.nan(sum(x))),]
#average SD
stdall=apply(stdmat2,2,mean)

#average est
estall=apply(Reduce(rbind,est),2,mean)

#average bias
apply(Reduce(rbind,est),2,mean) - truepar

index=which(sapply(std,function(x)!is.nan(sum(x))))

			
coveragelist=c()			
for(kk in 1:9){
coverage = sapply(1:numofsimu,function(i){

o=olist[[i]]
est = o[[1]]
std = o[[2]]

if (est[kk]+1.95*std[kk] > truepar[kk] && est[kk]-1.95*std[kk] < truepar[kk] ){
return(1)
}else{
return(0)
}

})
coveragelist = c(coveragelist,sum(coverage)/length(coverage))
}

tblo = rbind(truepar,estall,stdall,coveragelist)


			
  

  
  
  
  