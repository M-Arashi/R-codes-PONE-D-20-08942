rm(list=ls())
library(MASS)
library(KernSmooth)
library(snpMatrix)
library(leaps)
library(robustbase)
library(mvtnorm)
library(rrcov)
library(boot)
require(graphics)
library(Matrix)
library(glmnet)
library(linprog)
library(lpSolve)
library(plus)
library(parallel)
library(scalreg)
library(Rfit)



n=200;p=240;p1=floor(0.1*p);s=0.44;T=1000;rho=0.50
b=c(rep(rnorm(p1,mean=0.,sd=1)),rep(0,p-p1))
alpha=0.60
h=floor(n*alpha);h

R22=matrix(0,T,1);P1=matrix(0,T,1);P2=matrix(0,T,1)
Bhk=matrix(0,p,T);Bhrk=matrix(0,p,T);Bhrks=matrix(0,p,T)
msek=c();mserk=c();mserks=c()

for(tt in 1:T){

Z=X=matrix(0,n,p)
for(j in 1:p) Z[,j]=rnorm(n, j, 1)
for(j in 1:p) X[,j]=sqrt(1-rho^2)*Z[,j]+rho*Z[,p1]
kappa(t(X)%*%X)

V=matrix(0,h,h)
for(i in 1:h) 
for(j in 1:h)
V[i,j]=exp(-9*(abs(i-j)))
V=s*V

ep1=mvrnorm(1,rep(0,h),V)
ep2=rchisq(n-h, df=1, ncp=8)
ep=c(ep1,ep2)
y=X%*%b+ep

#bestlasso <- cv.glmnet(X,y,alpha=1,nfolds=10)
#Coef_lasso <- coef(bestlasso, s="lambda.min")[-1]
#X1=x1_las<- X[,which(!Coef_lasso %in% 0)]
#dim(x1_las)



C=t(X)%*%X#eigen(C)$values
Ck=function(k){
a=C+k*diag(p)
return(a)
}


###Ridge Estimator#################################################
bhk=function(k) solve(C+k*diag(p))%*%t(X)%*%y

L1=function(k) X%*%solve(C+k*diag(p))%*%t(X)
GCV1=function(k) (t(y-L1(k)%*%y)%*%(y-L1(k)%*%y)/n)/(1-sum(diag(L1(k)))/n)^2

gcv1=function(kd){
k=kd[1]
a=GCV1(k)
return(a)}


kopt=optim(c(0.3),gcv1,method="L-BFGS-B",lower=c(.001),upper=c(900.9))$par

Bhk[,tt]=bhk(kopt)
r1=function(k) t(bhk(k)-b)%*%(bhk(k)-b)
msek[tt]=r1(kopt)

#Rank-based Ridge Estimator######################################################

fit<-rfit(y~-1+X[,1:p1],scores=wscores)
yhatr=fit$fitted.values
res=fit$residuals
tau=fit$tauhat

r=rank(y)/(n+1)
a=getScores(wscores,r)
sigma_a=1/n*sum(a^2)
Xk=function(k) solve(C/n+k*diag(p))-k*solve(C/n+k*diag(p))%*%solve(C/n+k*diag(p))
Rn=function(k) 1/sigma_a*t(a)%*%X%*%solve(C/n+k*diag(p))%*%solve(Xk(k))%*%solve(C/n+k*diag(p))%*%t(X)%*%a

bhrk=function(k) solve(C+k*diag(p))%*%t(X)%*%yhatr
L2=function(k) 2*tau*X%*%solve(C+k*diag(p))%*%t(X)
GCV2=function(k) (t(y-L2(k)%*%y)%*%(y-L2(k)%*%y)/n)/(1-sum(diag(L2(k)))/n)^2

gcv2=function(kd){
k=kd[1]
a=GCV2(k)
return(a)}


kdopt2=optim(c(0.5),gcv2,method="L-BFGS-B",lower=c(0.0001),upper=c(1000))$par

Bhrk[,tt]=bhrk(kdopt2)
r2=function(k) t(bhrk(k)-b)%*%(bhrk(k)-b)
mserk[tt]=r2(kdopt2)


#Stein-type Shrinkage Estimator######################################################

bhrks=function(k,d) (1-d/Rn(k))[1,1]*bhrk(k)


L=function(k,d) (1-d/Rn(k))[1,1]*2*tau*X%*%solve(C+k*diag(p))%*%t(X)
GCV=function(k,d) (t(y-L(k,d)%*%y)%*%(y-L(k,d)%*%y)/n)/(1-sum(diag(L(k,d)))/n)^2

gcv=function(kd){
k=kd[1]
d=kd[2]
a=GCV(k,d)
return(a)}


kdopt=optim(c(0.5,0.5),gcv,method="L-BFGS-B",lower=c(0.09,.09),upper=c(1100,98.9))$par
kdopt

Bhrks[,tt]=bhrks(kdopt[1],kdopt[2])
r3=function(k,d) t(bhrks(k,d)-b)%*%(bhrks(k,d)-b)
mserks[tt]=r3(kdopt[1],kdopt[2])
 

##############OLSE
C=t(X)%*%X
Ck=function(k){
a=C+k*diag(p)
return(a)
} 
ols=summary(lm(y~X));
F=ols$fstatistic;
f=qf(0.95,p,n-p,ncp=0,lower.tail=TRUE)
F1=F[[1]];
P1[tt,]=I(F1>f)



 

xi=qchisq(0.95,p,ncp=0,lower.tail=TRUE,log.p=FALSE)
P2[tt]=I(Rn(kdopt2)>xi)

}
date()
R22;R2=sort(abs(R22))
c(mean(P1),mean(P2))

round(c(mean(msek),mean(mserk),mean(mserks)),4)
round(mean(msek)/c(mean(msek),mean(mserk),mean(mserks)),4)






######################################
