# insurance <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/insurance.txt",sep=",")
# kidney <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/kidney.txt",sep=",",header=T)
# ncog <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/ncog.txt",sep=" ",header=T)
# 
# 
# write.table(insurance,file="/Users/marco/Desktop/insurance.txt",sep=" ")
# write.table(kidney,file="/Users/marco/Desktop/kidney.txt",sep=" ")
# write.table(cancer,file="/Users/marco/Desktop/cancer.txt",sep=" ")
# write.table(sleep,file="/Users/marco/Desktop/sleep.txt",sep=" ")
rm(list=ls())


setwd("/Users/marco/Library/Mobile Documents/com~apple~CloudDocs/Columbia/teaching/Statistical Infererence and Modelling/Homeworks/HW2")



##################################################
### Exercise 1
##################################################

# cancer data
ncog <- read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/ncog.txt",sep=" ",header=T)


# 2. 
library(survival)
library(ggfortify)
autoplot(survfit(Surv(t, d)~arm,data=ncog))
dev.print(device=pdf,"/Users/marco/Library/Mobile Documents/com~apple~CloudDocs/Columbia/teaching/Statistical Infererence and Modelling/Homeworks/HW2/HW2surv1.pdf",onefile=FALSE)


# 3, 4
fit.exp <- survreg(Surv(t,d)~arm, dist="exp",data=ncog)
summary(fit.exp)


# 5
fit<- survfit(Surv(t, d)~arm,data=ncog)
par(mfrow=c(1,1))
plot(fit,col=c("darkred","darkgreen"),xlab="Days",conf.int=TRUE,main="Survival times")

x <- seq(from=0,to=2300,by=1)
lines(x,1-pexp(x,exp(-coef(fit.exp)[1])),col="darkred")
lines(x,1-pexp(x,exp(-sum(coef(fit.exp)))),col="darkgreen")
legend(1600, .9, legend=c("Arm A", "Arm B"), fill=c( "darkred","darkgreen"), cex=1.3)
dev.print(device=pdf,"/Users/marco/Library/Mobile Documents/com~apple~CloudDocs/Columbia/teaching/Statistical Infererence and Modelling/Homeworks/HW2/HW2surv2.pdf",onefile=FALSE)

#6
fit.wei <- survreg(Surv(t,d)~arm,data=ncog)
summary(fit.wei)
2*(fit.wei$loglik[2]-fit.exp$loglik[2])
pchisq(2*(fit.wei$loglik[2]-fit.exp$loglik[2]),df = 1)


par(mfrow=c(1,1))
plot(fit,col=c("darkred","darkgreen"),xlab="Days",conf.int=TRUE,main="Survival times")

x <- seq(from=0,to=2300,by=1)
lines(x,1-pexp(x,exp(-coef(fit.exp)[1])),col="darkred")
lines(x,1-pexp(x,exp(-sum(coef(fit.exp)))),col="darkgreen")
lines(x,1-pweibull(x,fit.wei$scale,exp(coef(fit.wei)[1])),lty=3,col="blue")
lines(x,1-pweibull(x,fit.wei$scale,exp(sum(coef(fit.wei)))),lty=4,col="blue")
legend(1600, .9, legend=c("Arm A", "Arm B","Weibull fit"), fill=c("darkgreen", "darkred","blue"), cex=1.3)
dev.print(device=pdf,"/Users/marco/Library/Mobile Documents/com~apple~CloudDocs/Columbia/teaching/Statistical Infererence and Modelling/Homeworks/HW2/HW2surv3.pdf",onefile=FALSE)



##################################################
### Exercise 2
##################################################
library(SMPracticals)
library(MissMech)

# Student scores data
x1 <- c(NA,53,51,NA,NA,NA,44,49,30,NA,NA,42,NA,NA,NA,17,39,48,46,30,NA,NA)
x2 <- c(63,61,67,69,69,49,61,41,69,59,40,60,63,55,49,53,46,38,40,34,30,26)
x3 <- c(65,72,65,53,61,62,52,61,50,51,56,54,53,59,45,57,46,41,47,43,32,15)
x4 <- c(70,64,65,53,55,63,62,49,52,45,54,49,54,53,48,43,32,44,29,46,35,20)
x5 <- c(63,73,NA,53,45,62,NA,NA,45,51,NA,NA,NA,NA,NA,51,NA,33,NA,18,21,NA)  



StScore <-cbind(x1,x2,x3,x4,x5) 
write.table(StScore,file="scores.txt")

marks=as.matrix(mathmarks)
rand=sample(66)
ind=c(5,6,7,17,18,20,21,22,31,35,36,39,40,41,47,52,60,62,71,74,86,87)
more=(1:88)[-ind] # observations from mathmarks not used in StScore not used 

# add one first observation of mathmarks
StScore1=rbind(StScore,marks[more[1],])
M=as.matrix(StScore1) 
unname(M)

# mean imputation
mean.impute<-function(M)
{
  indNA1=which(is.na(M[,1])==T)
  indNA5=which(is.na(M[,5])==T)
  m1=mean(M[-indNA1,1])
  m5=mean(M[-indNA5,5])
  impM=M
  impM[indNA1,1]=m1
  impM[indNA5,5]=m5
  impM
}
# bootstrap estimates with mean imputation
set.seed(1)
cov.boot <- array(0,dim=c(5,5,1000))
for(b in 1:1000)
{ 
  ind=sample(23,replace=T)
  cov.boot[,,b]=cov(mean.impute(M[ind,]))
}

###### 1
EMS=Mls(M)
pairS=cov(M,use="pairwise")
compS=cov(M,use="complete")
impS=cov(mean.impute(M))



bootimpS<-matrix(0,5,5)
for(i in 1:5)
{
  for(j in 1:5)
  {
    bootimpS[i,j]=mean(cov.boot[i,j,])
  }
}
bootimpS


# 2

eigComp=eigen(compS)
h.lambda=eigComp$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/22
h.lambda+1.96*h.lambda*sqrt(2)/22

eigPair=eigen(pairS)
h.lambda=eigPair$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/22
h.lambda+1.96*h.lambda*sqrt(2)/22

eigImp=eigen(impS)
h.lambda=eigImp$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/22
h.lambda+1.96*h.lambda*sqrt(2)/22



eigEM=eigen(EMS)
h.lambda=eigEM$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/22
h.lambda+1.96*h.lambda*sqrt(2)/22

boot.lambda=NULL
for(b in 1:1000)
{
  boot.lambda[b]=eigen(cov.boot[,,b])$values[1]
}
q.ind=order(boot.lambda)
q.025=boot.lambda[q.ind[25]]
q.975=boot.lambda[q.ind[975]]
q.025
q.975

# 3
eigCov=eigen(cov(mathmarks))
h.lambda=eigCov$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/88
h.lambda+1.96*h.lambda*sqrt(2)/88

# using only the data in  Student Score 
eigCov=eigen(cov(mathmarks[c(1,ind),]))
h.lambda=eigCov$values[1]
h.lambda-1.96*h.lambda*sqrt(2)/88
h.lambda+1.96*h.lambda*sqrt(2)/88



#########################################
### Exercise 4
#########################################

library(timeDate)
CentralPark <- read.csv("CentralPark.csv",header=T)
## Clean data
CP <- subset(CentralPark,select = c('DATE','SNWD'))
CP$time <- as.character(CP$DATE)
CP$time <- timeDate(CP$time,format = '%m/%d/%Y')
keep <- data.frame()
for(i in 1:dim(CP)[1]){
  if(as.numeric(format(CP$time[i],'%Y'))>1912){
    if(!(format(CP$time[i],'%Y') %in% c('1999','2000','2001','2002'))){
      if(format(CP$time[i],'%m')=='12'){
        if((as.numeric(format(CP$time[i],'%d'))>=17)&(as.numeric(format(CP$time[i],'%d'))<=31)){
          keep <- rbind(keep,CP[i,])
        }
      }
    }
  }
}

CP.mc <- list()
for(i in 1:101){
  CP.mc[[i]] <- as.numeric(keep$SNWD[(15*(i-1)+1):(15*i)]>50)
}
## First-order model
N1 <- matrix(0,nrow = 2,ncol = 2)
for(i in 1:101){ 
  for(j in 1:14){
    index1 <- as.numeric(CP.mc[[i]][j] == 1) + 1
    index2 <- as.numeric(CP.mc[[i]][j+1] == 1) + 1
    N1[index1,index2] = N1[index1,index2] + 1
  }
}
P1 <- N1/rowSums(N1)
P1
## Second-order model
N2 <- matrix(0,nrow = 4,ncol = 2)
for(i in 1:101){
  for(j in 1:13){
    index1 <- as.numeric(CP.mc[[i]][j] == 1) + 1 + 2*as.numeric(CP.mc[[i]][j+1] == 1)
    index2 <- as.numeric(CP.mc[[i]][j+2] == 1) + 1
    N2[index1,index2] = N2[index1,index2] + 1
  }
}
P2 <- N2/rowSums(N2)
P2
## Likelihood Ratio Test
loglikihood1 <- sum(log(P1)*N1)
loglikihood2 <- sum(log(P2)*N2)
pchisq(2*(loglikihood2-loglikihood1),df = 2)



#########################################
### Exercise 4
#########################################
library(TSA)
Huron <- read.table('Huron.txt')
## 1
Huron.linear <- data.frame(x = (1:nrow(Huron)),y = Huron$x)
fit.linear <- lm(y ~ x, data = Huron.linear)
summary(fit.linear)
plot(x = (1:98), y = residuals(fit.linear),type = 'l')
## 2
Huron.detrend <- residuals(fit.linear)
acf(Huron.detrend)
pacf(Huron.detrend)
## 3
arima(Huron.detrend, order = c(1,0,0))
arima(Huron.detrend, order = c(2,0,0))
## 4
eacf(Huron.detrend)

res=armasubsets(Huron.detrend,nar=2,nma=2,ar.method = 'ols')
plot(res)

arima(Huron.detrend, order = c(1,0,1))
arima(Huron.detrend, order = c(1,0,2))


