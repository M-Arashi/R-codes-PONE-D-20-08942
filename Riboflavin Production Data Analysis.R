rm(list=ls())
library(MASS)
library(Matrix)
library(glmnet)
library(linprog)
library(lpSolve)
library(plus)
library(parallel)
library(scalreg)
library(hdi)



######### READ DATA
data(riboflavin)
data <- na.omit(riboflavin);
#design matrix
X <- as.matrix(data$x);  #if it is necessary#  X<- scale(X);
#response variable
y<-  as.matrix(data$y);


###############Bivariate boxplot for effective genes in the riboflavin production data set 
bestlasso <- cv.glmnet(X,y,alpha=1,nfolds=20)
plot(bestlasso)
Coef_lasso <- coef(bestlasso, s="lambda.min")[-1]
X1=x1_las<- X[,which(!Coef_lasso %in% 0)]
X2=x2_las<- X[,which(Coef_lasso %in% 0)]

dim(x1_las)
dim(x2_las)

n=71;p=4088;p1=41;p2=p-p1
X=cbind(X1,X2);dim(X1)


require(grid) 
require(lattice)
require(MVA)

u1<-cbind(y,X1[,1]);u2<-cbind(y,X1[,2]);u3<-cbind(y,X1[,3]);u4<-cbind(y,X1[,4]);u5<-cbind(y,X1[,5]);u6<-cbind(y,X1[,6]);u7<-cbind(y,X1[,7]);u8<-cbind(y,X1[,8]);u9<-cbind(y,X1[,9]);u10<-cbind(y,X1[,10]);
u11<-cbind(y,X1[,11]);u12<-cbind(y,X1[,12]);u13<-cbind(y,X1[,13]);u14<-cbind(y,X1[,14]);u15<-cbind(y,X1[,15]);u16<-cbind(y,X1[,16]);u17<-cbind(y,X1[,17]);u18<-cbind(y,X1[,18]);u19<-cbind(y,X1[,19]);u20<-cbind(y,X1[,20]);
u21<-cbind(y,X1[,21]);u22<-cbind(y,X1[,22]);u23<-cbind(y,X1[,23]);u24<-cbind(y,X1[,24]);u25<-cbind(y,X1[,25]);u26<-cbind(y,X1[,26]);u27<-cbind(y,X1[,27]);u28<-cbind(y,X1[,28]);u29<-cbind(y,X1[,29]);u30<-cbind(y,X1[,30]);
u31<-cbind(y,X1[,31]);u32<-cbind(y,X1[,32]);u33<-cbind(y,X1[,33]);u34<-cbind(y,X1[,34]);u35<-cbind(y,X1[,35]);u36<-cbind(y,X1[,36]);u37<-cbind(y,X1[,37]);u38<-cbind(y,X1[,38]);u39<-cbind(y,X1[,39]);u40<-cbind(y,X1[,40]);
u41<-cbind(y,X1[,41])


#Fig3
par(mfrow = c(3, 5))
bvbox(u1, mtitle = "", xlab ="ARGF_at", ylab ="y")
bvbox(u2, mtitle = "", xlab ="DNAJ_at", ylab ="y")
bvbox(u3, mtitle = "", xlab ="GAPB_at", ylab ="y")
bvbox(u4, mtitle = "", xlab ="LYSC_at", ylab ="y")
bvbox(u5, mtitle = "", xlab ="ARGF_at", ylab ="y")
bvbox(u6, mtitle = "", xlab ="PKSA_at", ylab ="y")
bvbox(u7, mtitle = "", xlab ="PRIA_at", ylab ="y")
bvbox(u8, mtitle = "", xlab ="SPOIIAA_at", ylab ="y")
bvbox(u9, mtitle = "", xlab ="SPOVAA_at", ylab ="y")
bvbox(u10, mtitle = "", xlab ="THIK_at", ylab ="y")
bvbox(u12, mtitle = "", xlab ="XHLB_at", ylab ="y")
bvbox(u13, mtitle = "", xlab ="YACN_at", ylab ="y")
bvbox(u14, mtitle = "", xlab ="YBFI_at", ylab ="y")
bvbox(u15, mtitle = "", xlab ="YCDH_at", ylab ="y")
bvbox(u16, mtitle = "", xlab ="YCGO_at", ylab ="y")
bvbox(u17, mtitle = "", xlab ="YCKE_at", ylab ="y")
bvbox(u18, mtitle = "", xlab ="YCLB_at", ylab ="y")
bvbox(u19, mtitle = "", xlab ="YCLF_at", ylab ="y")
bvbox(u20, mtitle = "", xlab ="YDDH_at", ylab ="y")
bvbox(u21, mtitle = "", xlab ="YDDK_at", ylab ="y")
bvbox(u22, mtitle = "", xlab ="YEBC_at", ylab ="y")
bvbox(u23, mtitle = "", xlab ="YFHE_r_at", ylab ="y")
bvbox(u24, mtitle = "", xlab ="YFII_at", ylab ="y")
bvbox(u25, mtitle = "", xlab ="YFIO_at", ylab ="y")
bvbox(u26, mtitle = "", xlab ="YFIR_at", ylab ="y")
bvbox(u27, mtitle = "", xlab ="YHDS_r_at", ylab ="y")
bvbox(u28, mtitle = "", xlab ="YKBA_at", ylab ="y")
bvbox(u29, mtitle = "", xlab ="YKVJ_at", ylab ="y")
bvbox(u30, mtitle = "", xlab = "YLXW_at", ylab ="y")
bvbox(u31, mtitle = "", xlab = "YMFE_at", ylab ="y")
bvbox(u32, mtitle = "", xlab = "YOAB_at", ylab ="y")
bvbox(u33, mtitle = "", xlab = "YPGA_at", ylab ="y")
bvbox(u34, mtitle = "", xlab = "YQJT_at", ylab ="y")
bvbox(u35, mtitle = "", xlab = "YQJU_at", ylab ="y")
bvbox(u36, mtitle = "", xlab = "YRVJ_at", ylab ="y")
bvbox(u37, mtitle = "", xlab = "YTGB_at", ylab ="y")
bvbox(u38, mtitle = "", xlab = "YUID_at", ylab ="y")
bvbox(u39, mtitle = "", xlab = "YURQ_at", ylab ="y")
bvbox(u40, mtitle = "", xlab = "YXLE_at", ylab ="y")
bvbox(u41, mtitle = "", xlab = "YYBG_at", ylab ="y")
bvbox(u42, mtitle = "", xlab = "YYDA_at", ylab ="y")



m00=lm(y~X1)
#Fig1
par(mfrow=c(2,3),mar=c(4,6,4,4))
mlab <- "Observations"
plab <- "Standardized Residuals"
plot(rstudent(m00),xlab=mlab,ylab=plab)
abline(h=c(-2,2),lty=3,col="darkblue",lwd=2)
identify(rstandard(m00)+.25)
plot(m00)





library(robustHD)
frac <- seq(0.2, 0.05,by =-0.05)
fitSparseLTS=sparseLTS(X,y,lambda=0.05,mode="fraction")
a<-plot(fitSparseLTS,method =  "diagnostic")
a[1]
a[2]
a[3]
a[4]


diagnosticPlot(fitSparseLTS, id.n = NULL)
bhr=as.vector(coef(fitSparseLTS, zeros = FALSE))
g=(fitSparseLTS$best)
x1=X[g,]
dim(x1)






library(robustbase)
mRobust=ltsReg(y~X1[,1]+X1[,2]+X1[,3]+X1[,4]+X1[,5]+X1[,6]+X1[,7]+X1[,8]+X1[,9]+X1[,10]+X1[,11]+X1[,12]+X1[,13]+X1[,14]+X1[,15]+X1[,16]+X1[,17]+X1[,18]+X1[,19]+X1[,20]+X1[,21]+X1[,22]+X1[,23]+X1[,24]+X1[,25]+X1[,26]+X1[,27]+X1[,28]+X1[,29]
+X1[,30]+X1[,31]+X1[,32]+X1[,33]+X1[,34],model=TRUE)
summary(mRobust)
plot(mRobust, which = "rqq")
plot(mRobust, which = "all")



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


kopt=optim(c(0.3),gcv1,method="L-BFGS-B",lower=c(.001),upper=c(90.9))$par
kdopt
gc1=gcv1(kopt)

bh1=bhk(kopt)
RSS1=t(y-X1%*%bh1)%*%(y-X1%*%bh1);RSS1
ssy=sum(y^2)-n*mean(y)^2
R12=1-gcv1(kopt)/ssy;R12


#Rank-based Ridge Estimator######################################################


library(Rfit)
fit<-rfit(y~X1,scores=wscores)
#summary(fit)
#Bhr=fit$coef
yhatr=fit$fitted.values
res=fit$residuals
tau=fit$tauhat
plot(y~yhatr)

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


kdopt2=optim(c(.5),gcv2,method="L-BFGS-B",lower=c(.001),upper=c(1000))$par
kdopt2
gcv2(kdopt2)
bh2=bhk(kdopt)
RSS2=t(y-X%*%bh2)%*%(y-X%*%bh2);RSS2
ssy=sum(y^2)-n*mean(y)^2
R22=1-gcv2(kdopt2)/ssy;R22
#Fig2
rr=as.vector(y-X%*%bh2)
qqnorm(rr)
qqline(rr)


#Stein-type Shrinkage Estimator######################################################

bhrks=function(k,d) (1-d/Rn(k))[1,1]*bhrk(k)


L=function(k,d) (1-d/Rn(k))[1,1]*2*tau*X%*%solve(C+k*diag(p))%*%t(X)
GCV=function(k,d) (t(y-L(k,d)%*%y)%*%(y-L(k,d)%*%y)/n)/(1-sum(diag(L(k,d)))/n)^2

gcv=function(kd){
k=kd[1]
d=kd[2]
a=GCV(k,d)
return(a)}


kdopt=optim(c(.5,.5),gcv,method="L-BFGS-B",lower=c(.001,.002),upper=c(1000,98.9))$par
kdopt
gcv(kdopt)

R32=1-gcv(kdopt)/ssy;R32

round(c(R12,R22,R32),6)


#Fig4
options(digits=15)
t1=(4600:4799)/10000;t1
t2=(100:399)/100000;t2
z=matrix(ncol=300,nrow=200)
for(i in 1:200){
      for(j in 1:300){
           z[i,j]=GCV((i+4599)/10000,(j+99)/100000)}}

par(cex.axis=0.5,cex.lab=.9)
persp(x=t1[80:94],y=t2[110:150],z[80:94,110:150],xlab="Ridge parameter",ylab="Shrinkage parameter",zlab="GCV",ticktype="detailed"
,theta=50,phi=20,expand=1,axes=TRUE,box=TRUE,col="lightgrey")

contour(x=t1,y=t2,z,drawlabels=FALSE,
xlab="Ridge parameter",ylab="Shrinkage parameter")

u3=c()
t3=(4600:4799)/10000
for(i in 1:200)
{
u3[i]=GCV((i+4599)/10000,kdopt[2])
}
u3;min(u3)

u4=c()
t4=(100:399)/100000
for(j in 1:300)
{
u4[j]=GCV(kdopt[1],(j+99)/100000)
}
u4;min(u4)

mlab <- expression(Ridge~parameter~k)
plab <- expression(GCV)
plot(t3[80:94],u3[80:94],type="l",xlab= mlab,ylab = plab,ylim=c(min(u3[80:94]),max(u3[80:94])+.000000001),pch=".")
ex.cs<-expression(paste(italic(d==0.002332)))
utils::str(legend(.46885,u3[80]+.000000001, ex.cs,col=c(1)))

mlab <- expression(Shrinkage~parameter~d)
plab <- expression(GCV)
plot(t4[110:150],u4[110:150],type="l",xlab= mlab,ylab = plab,pch="*")
ex.cs<-expression(paste(italic(k==0.468759)))
utils::str(legend(0.00236,u4[110], ex.cs,col=c(1)))

#######################################################


