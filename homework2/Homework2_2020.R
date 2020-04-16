rm(list = ls())

##################################################
### Exercise 1
##################################################


tp <- read.table("transplant.txt")
names(tp) <- c("t","type","d")
tp$type<- factor(tp$type)

# 2. 
library(survival)
library(ggfortify)
autoplot(survfit(Surv(t, d)~type,data=tp))


# 3, 4
fit.exp <- survreg(Surv(t,d)~type, dist="exp",data=tp)
summary(fit.exp)


# 5
fit<- survfit(Surv(t, d)~type,data=tp)
par(mfrow=c(1,1))
plot(fit,col=c("orange","blue"),xlab="Months",conf.int=TRUE,main="Survival times")

x <- seq(from=0,to=70,by=1)
lines(x,1-pexp(x,exp(-coef(fit.exp)[1])),col="orange")
lines(x,1-pexp(x,exp(-sum(coef(fit.exp)))),col="blue")
legend(40, 1, legend=c("Type 1", "Type 2"), fill=c( "orange","blue"), cex=0.8)


#6
fit.wei <- survreg(Surv(t,d)~type,data=tp)
summary(fit.wei)
2*(fit.wei$loglik[2]-fit.exp$loglik[2])
pchisq(2*(fit.wei$loglik[2]-fit.exp$loglik[2]),df = 1)


par(mfrow=c(1,1))
plot(fit,col=c("orange","blue"),xlab="Months",conf.int=TRUE,main="Survival times")

x <- seq(from=0,to=70,by=1)
lines(x,1-pexp(x,exp(-coef(fit.exp)[1])),col="orange")
lines(x,1-pexp(x,exp(-sum(coef(fit.exp)))),col="blue")
lines(x,1-pweibull(x,fit.wei$scale,exp(coef(fit.wei)[1])),lty=3,col="red")
lines(x,1-pweibull(x,fit.wei$scale,exp(sum(coef(fit.wei)))),lty=4,col="magenta")
legend(40, 1, legend=c("Type 1", "Type 2","Weibull fit Type 1", "Weibull fit Type 2"),
       fill=c("orange", "blue","red", "magenta"), cex=0.8)


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
h.lambda-1.96*h.lambda*sqrt(2/22)
h.lambda+1.96*h.lambda*sqrt(2/22)

eigPair=eigen(pairS)
h.lambda=eigPair$values[1]
h.lambda-1.96*h.lambda*sqrt(2/22)
h.lambda+1.96*h.lambda*sqrt(2/22)

eigImp=eigen(impS)
h.lambda=eigImp$values[1]
h.lambda-1.96*h.lambda*sqrt(2/22)
h.lambda+1.96*h.lambda*sqrt(2/22)

eigEM=eigen(EMS)
h.lambda=eigEM$values[1]
h.lambda-1.96*h.lambda*sqrt(2/22)
h.lambda+1.96*h.lambda*sqrt(2/22)

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
h.lambda-1.96*h.lambda*sqrt(2/88)
h.lambda+1.96*h.lambda*sqrt(2/88)

# using only the data in  Student Score 
eigCov=eigen(cov(mathmarks[c(1,ind),]))
h.lambda=eigCov$values[1]
h.lambda-1.96*h.lambda*sqrt(2/88)
h.lambda+1.96*h.lambda*sqrt(2/88)



#########################################
### Exercise 3
#########################################

library(timeDate)
CentralPark <- read.csv("CentralPark.csv",header=T)
## Clean data
CP <- subset(CentralPark,select = c('DATE','PRCP'))
CP$DATE <- as.Date(CP$DATE, "%m/%d/%Y")
CP <- subset(CP, format.Date(DATE, "%m")=="07" )
CP.mc<- as.numeric(CP$PRCP>1.5)

##### 3
## First-order model
N1 <- matrix(0,nrow = 2,ncol = 2)
for(i in 1:(length(CP.mc)-1)){ 
  if (i%%31 != 0){
    index1 <- as.numeric(CP.mc[i] == 1) + 1
    index2 <- as.numeric(CP.mc[i+1] == 1) + 1
    N1[index1,index2] = N1[index1,index2] + 1
  }
}

P1 <- N1/rowSums(N1)
P1

#### 4
## significance
1-pnorm((P1[1,1]-P1[2,2])/(P1[1,1]*(1-P1[1,1])/(N1[1,1]+N1[1,2])+P1[2,2]*(1-P1[2,2])/(N1[2,1]+N1[2,2]))^0.5)

#### 5
## Second-order model
N2 <- matrix(0,nrow = 4,ncol = 2)
for(i in 1:(length(CP.mc)-2)){ 
  if ((i%%31 != 0) & (i%%31 != 30)){
    index1 <- as.numeric(CP.mc[i] == 1) + 1 + 2*as.numeric(CP.mc[i+1] == 1)
    index2 <- as.numeric(CP.mc[i+2] == 1) + 1
    N2[index1,index2] = N2[index1,index2] + 1
  }
}
P2 <- N2/rowSums(N2)
P2
## Likelihood Ratio Test
loglikihood1 <- sum(log(P1)*N1)
loglikihood2 <- sum(log(P2)*N2)
1-pchisq(2*(loglikihood2-loglikihood1),df = 2)



#########################################
### Exercise 5
#########################################

#define tau to be the vecotor of times of death/relapse
tau <- tp[which(tp$d==1),]$t 
tau <- sort(tau)

A.df <- tp[which(tp$type==1),] #subset treatment A
B.df <- tp[which(tp$type==2),] #subset treatment B

#initialize
Z.num = 0
Z.den = 0

for(k in 1:length(tau)){
  #calculate relevant quantities
  y_k = sum(A.df$t == tau[k])
  n_A = sum(A.df$t >= tau[k])
  n_B = sum(B.df$t >= tau[k])
  n_d = sum(tp$t == tau[k])
  
  E_k = n_A*n_d/(n_A+n_B)
  if(n_A+n_B == 1){
    V_k = 0
  } else {
    V_k = n_A*n_B*n_d*(n_A+n_B-n_d)/((n_A+n_B)^2*(n_A+n_B-1)) 
  }
  
  #update numerator and denominator
  Z.num = Z.num + y_k - E_k
  Z.den = Z.den + V_k
}

Z.test = Z.num/sqrt(Z.den)
2*(1-pnorm(abs(Z.test))) #p-value








