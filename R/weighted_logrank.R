
# Note, UGH need to be carefule
# these functions are duplicated in "plottting_functions_v0.R"
# But here they are with lower case "w" in the name

N.weighted<-function(x,error,w=rep(1,length(error))){
  sum(w*(error<=x))
}

R.weighted<-function(x,error,w=rep(1,length(error))){
  sum(w*(error>=x))
}


################################################
### Start Kosorok's function surv.Rtest() ######
################################################

# Note: This section is from M. Kosorok's 
# code posted on his webpage

EG.testval<-function(x){
 xs<-max(x)
     xi<-min(x)
  
     if(abs(xs)>abs(xi)){test<-xs}
     else{test<-xi}
     y <- abs(test)
     return(y)
     }



sup.G<-function(x,m=10)
{
     k<-m:0
     (4/pi)*sum(((-1)^k)/(2*k+1)*exp(-(pi^2)*((2*k+1)^2)/(8*x^2)))
}
sup.g<-function(x,m=10)
{
     k<-m:0
     (pi/x^3)*sum(((-1)^k)*(2*k+1)*exp(-(pi^2)*((2*k+1)^2)/(8*x^2)))
}
cnorm<-function(z,thresh=3.6,delta=0.6,kk=4){
check<-F
if(z<0){
     z<-(-1)*z
     check<-T
}
if(z<thresh){
     out<-1-pnorm(z)
}
else{
     term<-1
     tally<-term
     if(kk>1){
          for(k in 1:(kk-1)){
               term<-(-1)*term*(2*k-1)/z^2
               tally<-tally+term
          }
     }
     out<-tally*dnorm(z)/z
     if(z<thresh+delta){
          x<-1-pnorm(z)
          out<-x+(z-thresh)*(out-x)/delta
     }
}
if(check){out<-1-out}
out
}
sup.inverse<-function(alpha,error=1e-8)
{
     x<-qnorm(1-alpha/4)
     temp<-max(1,2/x)
     m<-ceiling((x/pi)*sqrt(2*log(temp/(pi*error)))-0.5)
     if(m<0){m<-0}
     interror<-1
     while(interror>error)
     {
          yx<-sup.G(x,m=m)
          dg<-sup.g(x,m=m)
          delta<-(1-alpha-yx)/dg
          x<-x+delta
          interror<-sup.G(x)-(1-alpha)
     }
     x
}
sup.r<-function(alpha, beta, error=1e-8)
{
     time<-sup.inverse(alpha,error=error)
     y<-1-beta
     ml<-qnorm(1-alpha/2)+qnorm(1-beta)
     x<-ml
     delta<-1
     while(delta>error)
     {
          yx<-cnorm(time-x)+exp(2*x*time)*cnorm(time+x)
          dp<-dnorm(time-x)-exp(2*time*x)*dnorm(time+x)+2*time*exp(2*time*x)*cnorm(time+x)
          delta<-(y-yx)/dp
          x<-x+delta
     }
     (x/ml)^2    
}
surv.Rtest<-function (time, delta, group, rho=0, gamma=0, logrank=F, 
     error=1.0e-8) 
{    
     otime <- order(time)
     time <- time[otime]
     delta <- delta[otime]
     n<-length(time)
     if((rho+gamma)==0){
          weight<-rep(1,n)
     }
     else{
          km.left<-KM.left(time,delta)
          weight<-km.left^rho*(1-km.left)^gamma
     }
     group <- group[otime] - 1
     n2 <- sum(group)
     atrisk2 <- n2 - cumsum(group) + group
     n1 <- n-n2
     atrisk1 <- n1 - cumsum(1 - group) + 1 - group
     delta1 <- delta * (1 - group)
     delta2 <- delta * group
     
          
     y1 <- tapply(atrisk1, time, "max")
     y2 <- tapply(atrisk2, time, "max")
     d1 <- tapply(delta1, time, "sum")
     d2 <- tapply(delta2, time, "sum")
     weight<-tapply(weight,time,"max")
     w <- (y1 * y2)/(y1 + y2)
     
     
     # Modification to output Z(t)
     time.w<-c(0,time[which(w > 0)])
     zrg<-data.frame(cbind(time.w,time.w))
     colnames(zrg)<-c("time","Z.rg")              
     terms <-(d1/y1 - d2/y2)[w > 0]       
     temp<-y1+y2-1
     temp<-ifelse(temp<1,1,temp)
     cc<-1-(d1+d2-1)/temp
     vterms <- (cc*(d1 + d2)/(y1 + y2))[w > 0]
     weight<-weight[w > 0]
     w<-w[w > 0]
     terms <- (weight * w * terms)/sqrt(sum(weight^2 * w * vterms))
     temp<-c(0,cumsum(terms))       
        
    zrg$Z.rg<-c(0,cumsum(terms))
    
     xs<-max(temp)
     xi<-min(temp)
  
     if(abs(xs)>abs(xi)){test<-xs}
     else{test<-xi}
     x <- abs(test)
     m<-ceiling(max(c(1,(x*sqrt(2)/pi)*sqrt(max(c(1,log(1/(pi*error)))))-0.5)))
     p<-1-sup.G(x,m=m)
     out <- NULL
     out$test <- test
     out$p <- p
     out$zrg<-zrg
     if(logrank){
          x<-temp[length(temp)]
          out$test.logrank<-x
          out$p.logrank<-2*cnorm(abs(x))     
     }
     out
}
KM.left<-function(time,delta){
     n<-length(time)
     dtime<-tapply(time,time,"max")
     ddelta<-tapply(delta,time,"sum")
     dy<-tapply(rep(1,n),time,"sum")
     m<-length(dy)
     y<-rep(n,m)-c(0,cumsum(dy)[1:(m-1)])
     km<-1
     km.left<-rep(0,m)
     for(i in 1:m){
          km.left[i]<-km
          km<-km*(1-ddelta[i]/y[i])
     }
     out<-rep(0,n)
     for(i in 1:n){
          out[i]<-min(km.left[dtime==time[i]])
     }
     out
}

######################################################
################## End Kosorok's #####################
######################################################


wt.rg.logrank<-function(S,rho,gamma,tpoints=NULL,t.tau=NULL,w0.tau=0,w1.tau=1){
temp<-S[1:(length(S)-1)]
S.new<-c(1,temp) # S.KM(t-)
term1<-S.new^rho
term2<-(1-S.new)^gamma
wt<-term1*term2
if(!is.null(t.tau)){
# For optimal wt with known changepoint at t.tau
wt<-ifelse(tpoints<=t.tau,w0.tau,w1.tau)  
}
return(wt)
}



weight.Logrank<-function(time,delta,z,rho=0,gamma=0){
is.sorted<-!is.unsorted(time)
if(!is.sorted){
id<-order(time); time<-time[id]; delta<-delta[id]; z<-z[id]
}
at.points<-sort(unique(c(time[delta==1])))

U0<-time[which(z==0)]
D0<-delta[which(z==0)]

# Control group
# Risk and Counting processes
risk.z0<-unlist(lapply(as.list(at.points),R.weighted,error=U0))
counting<-unlist(lapply(as.list(at.points),N.weighted,error=U0[D0==1]))
N.z0<-counting
dN.z0<- diff(c(0, counting))

U1<-time[which(z==1)]
D1<-delta[which(z==1)]

# Control group
# Risk and Counting processes
risk.z1<-unlist(lapply(as.list(at.points),R.weighted,error=U1))
counting<-unlist(lapply(as.list(at.points),N.weighted,error=U1[D1==1]))
N.z1<-counting
dN.z1<- diff(c(0, counting))

dN.pooled<-dN.z0+dN.z1
risk.pooled<-risk.z0+risk.z1

##########################
# Pooled K-M estimator
dN.Risk<-ifelse(risk.pooled>0,dN.pooled/risk.pooled,0)
S.pool<-cumprod(1-dN.Risk)
temp<-S.pool[1:(length(S.pool)-1)]
S.pool<-c(1,temp) # S.pool(t-)
w<-(S.pool^rho)*((1-S.pool)^gamma)

# Checking S.pool
#km.pool<-survfit(Surv(time,delta)~1)
#temp<-summary(km.pool,c(at.points))$surv
#temp<-temp[1:(length(temp)-1)]
#S.pool.check<-c(1,temp) # S.pool(t-)
#plot(at.points,S.pool,type="s",lty=1,col="grey",lwd=3)
#lines(at.points,S.pool.check,lty=2,type="s",col="red",lwd=2)

K<-ifelse(risk.pooled>0,w*(risk.z0*risk.z1)/risk.pooled,0.0)

term0<-sum(ifelse(risk.z0>0,(K/risk.z0)*dN.z0,0.0))
term1<-sum(ifelse(risk.z1>0,(K/risk.z1)*dN.z1,0.0))
# lr>0 --> increased risk for "Control"
lr<-term1-term0

# variance
h0<-ifelse(risk.z0==0,0,(K^2/risk.z0))
h1<-ifelse(risk.z1==0,0,(K^2/risk.z1))
dJ<-ifelse(risk.pooled==1,0,(dN.pooled-1)/(risk.pooled-1))
dL<-ifelse(risk.pooled==0,0,dN.pooled/risk.pooled)
sig2s<-(h0+h1)*(1-dJ)*dL
sig2<-sum(sig2s)

Z.lr<-lr/sqrt(sig2)
result<-list(Z.lr=Z.lr,pval=1-pchisq(Z.lr^2,1),score=lr,sig2=sig2)
return(result)
}

# Evaluate tests 
Z.stat<-function(x,stat="overall"){
if(stat=="overall") return(x[length(x)])
if(stat=="renyi") return(max(abs(x)))
}

# For calculating individual p-values
# for Z statistics based on resampling 
# (SM = synthetic martingale)
# x[1] will be the observed test
# x[-1] the resamples
pvalZ.SM<-function(x){
return(mean(x[-1]>=x[1]))
}


# Gill's K-class of weighted log-rank
# tests. Includes: 
# FH rho-gamma; Gehan; Tarone-ware
ZK.logrank<-function(r,K,Y1,D1,Y0,D0){
#points at which r=r.tau>0
Ybar<-Y1+Y0
Dbar<-D1+D0
# D1/Y1 = dN1(.)/Y1 representing
# differential of Nelson-Aalen
# estimator of CHF
term0<-ifelse(Y0>0,K*(D0/Y0),0.0)
term1<-ifelse(Y1>0,K*(D1/Y1),0.0)
temp<-term0-term1
wK<-cumsum(temp[r>0])

# Variance term
H1<-ifelse(Y1>0,K/Y1,0.0)
H0<-ifelse(Y0>0,K/Y0,0.0)

h1<-ifelse(Y1>0,(K^2/Y1),0.0)
h0<-ifelse(Y0>0,(K^2/Y0),0.0)
temp<-Ybar-1
ybar.mod<-ifelse(temp<1,1,temp)

Qbar<-ifelse(ybar.mod>0,(Dbar-1)/ybar.mod,0.0)
dHbar<-ifelse(Ybar>0,Dbar/Ybar,0.0)

sig2s<-(h1+h0)*(1-Qbar)*dHbar

sig2.wK<-cumsum(sig2s[r>0])
ZK<-ifelse(sig2.wK>0,wK/sqrt(sig2.wK),0.0)
return(list(wK=wK,sig2K=sig2.wK,H1=H1,H0=H0,ZK=ZK))

}

rg.count<-function(x,rg0){
  ifelse(x[1]==rg0[1] & x[2]==rg0[2],1,0)
}



combo2.logrank<-function(time,delta,z,draws=0,rho.gamma=c(0,0),include.gehan=FALSE,L.type="observed",seed.value=NULL,
                         details=FALSE,checking=TRUE,checking.KM=FALSE,draws.status=TRUE,draws.out=FALSE){
# Test H0: F1=F0 vs H1:F1>F0
# Default is combination of 1-sided wLR tests
# For sup|Z(t)| the tests are 2-sided

  if (missing(z))
        stop("Treatment missing")
    n<-length(time)
    if (n!=length(z)) stop("lengths of 'time' and 'z' differ")
    if (n!=length(delta)) stop("lengths of 'time' and 'delta' differ")
        if ( !identical(sort(as.integer(unique(z))),0:1) )
        stop("Treatment (z) values must be 0,1")
        if (is.vector(rho.gamma)&&(length(rho.gamma)==2)) {
        rho.gamma<-list(rho.gamma)
    } else {
        if (!is.list(rho.gamma)) stop("'rho.gamma' must be a vector of length 2 or a list of vectors of length 2")
    }
    rho.gamma.mat<-matrix(unlist(rho.gamma),ncol=2,byrow=TRUE) # each row one pair (rho,gamma)
    m<-nrow(rho.gamma.mat)
# If m>1 then will combine over the 
# specified (rho,gamma) family of weighting schemes

    # Count (r,g) combinations for "max Cox"
    rgs.mat<-cbind(rho.gamma.mat,rep(0,nrow(rho.gamma.mat))) # Initiate counts
    colnames(rgs.mat)<-c("rho","gamma","counts")
    
    
# Control
U0<-time[z==0]
delta0<-delta[z==0]
get0<-order(U0)
U0<-U0[get0]
delta0<-delta0[get0]
# Treatment
U1<-time[z==1]
delta1<-delta[z==1]
get1<-order(U1)
U1<-U1[get1]
delta1<-delta1[get1]

jump.points<-c(sort(unique(time)))

Y0<-c(unlist(lapply(as.list(jump.points),R.weighted,error=U0)))
counting<-unlist(lapply(as.list(jump.points),N.weighted,error=U0[delta0==1]))
counting<-c(0, counting)
N0<-counting
DN0<-diff(counting)

Y1<-c(unlist(lapply(as.list(jump.points),R.weighted,error=U1)))
counting<-unlist(lapply(as.list(jump.points),N.weighted,error=U1[delta1==1]))
counting<-c(0, counting)
N1<-counting
DN1<-diff(counting)

############################################################
# Check KM at jump.points
if(checking.KM){
#par(mfrow=c(1,2))
S0.KM<-cumprod(1-(DN0/Y0))
temp0<-survfit(Surv(U0,delta0)~1)
# match to points in S0.KM(.)
common.points<-intersect(jump.points,temp0$time)
at0<-jump.points %in% common.points
S0.match<-S0.KM[at0]
at0.check<-temp0$time %in% common.points
S0.check<-summary(temp0,c(temp0$time[at0.check]))
diff<-round(c(S0.check$surv-S0.match),2)
if(max(abs(diff))>0){
plot(temp0,col="grey",lwd=5,conf.int=FALSE)
lines(jump.points,S0.KM,type="s",lty=2,col="black",lwd=2)
stop("Problem in KM Estimator: S0.match vs surv()")
}
S1.KM<-cumprod(1-(DN1/Y1))
temp1<-survfit(Surv(U1,delta1)~1)
# match to points in S1.KM(.)
common.points<-intersect(jump.points,temp1$time)
at1<-jump.points %in% common.points
S1.match<-S1.KM[at1]
at1.check<-temp1$time %in% common.points
S1.check<-summary(temp1,c(temp1$time[at1.check]))
diff<-round(c(S1.check$surv-S1.match),2)
if(max(abs(diff))>0){
print(cbind(jump.points[at1],S1.match,temp1$time[at1.check],S1.check$surv,S1.check$surv-S1.match))
plot(temp1,col="grey",lwd=5,conf.int=FALSE)
lines(jump.points,S1.KM,type="s",lty=2,col="black",lwd=2)
stop("Problem in KM Estimator: S1.match vs surv()")
   }
}  # End checking KM

Ybar<-Y0+Y1
Nbar<-N0+N1
DNbar<-DN0+DN1
DN.Risk<-DNbar/Ybar

rbar<-(Y1*Y0/Ybar)

# we will restrict to time points at which
# rbar>0
r.points<-c(0,jump.points[rbar>0])

S.pool<-cumprod(1-(DNbar/Ybar))

n0<-length(U0)
n1<-length(U1)

# Placeholder for Gehan (or other K-class)
if(include.gehan){
K<-Y1*Y0  # Gehan Statistic
temp<-ZK.logrank(r=rbar,K=K,Y1=Y1,D1=DN1,Y0=Y0,D0=DN0)
w.gehan<-temp$wK
sig2.gehan<-temp$sig2K
Z.gehan<-temp$ZK
H1.gehan<-temp$H1
H2.gehan<-temp$H2
rm("temp")
}

# Each m=1,...,M rows of rho.gamma.mat corresponds
# to a (rho,gamma) member
# rho=-1 and gamma=-1 represent RMST

if(draws==0) G0.draws<-NULL; G1.draws<-NULL

# Sampling distribution under the null!
# Calculate individual Z.rg and combination test
# Initialize g*'s for each treatment group
if(draws>0){
if(!is.null(seed.value)) set.seed(seed.value)
G0.draws<-matrix(rnorm(draws*n0),ncol=draws)
G1.draws<-matrix(rnorm(draws*n1),ncol=draws)
}
# Correlation is accounted for by using the same g's
# accross the weighted-log rank statistics

# The "overall" final statistic
Z.rg.obs.vec<-rep(NA,m)
Z.rg.star.mat<-matrix(NA,nrow=m,ncol=draws)

# The Renyi-type statistic
Z.rg2.obs.vec<-rep(NA,m)
Z.rg2.star.mat<-matrix(NA,nrow=m,ncol=draws)

# Each column will store Z.rg^* for d=1,...,draws
# These represent distribution under the null!


# Calculate RMST
rhos<-rho.gamma.mat[,1]
gammas<-rho.gamma.mat[,2]
fit.RMST<-NULL
if(any(rhos==-1 & gammas==-1)){
fit.RMST<-KM.MeanTrunc(time=time,delta=delta,z=z,L=NULL,L.type=L.type,draws=0,G1.draws=G1.draws,G0.draws=G0.draws,get.band=FALSE,details=details)
Z.rmst.star<-fit.RMST$Zdstar
Z.rmst<-fit.RMST$Z.asy
}

##########################################################
# NOTE: rho,gamma calculation here!
# Any additional (rho,gamma) calculations can
# be implemented here:
# Currently, this implements max(Z1,Z2,...,Zk)
# where the Z's indicate k 1-sided weighted log-rank
# statistics (k=(rho,gamma) pairs).
# Also implements max(sup_t|Z1(t)|, ..., sup_t|Zk(t)|).
# The latter is referred to as "Renyi".
############################################################

for(rg in 1:m){

rhogamma<-rho.gamma.mat[rg,]

if(rhogamma[1]>=0){
w.rg<-wt.rg.logrank(S=S.pool,rho=rhogamma[1],gamma=rhogamma[2])

# rho-gamma LR process
K<-ifelse(Ybar>0,w.rg*(Y1*Y0)/Ybar,0.0)
temp<-ZK.logrank(r=rbar,K=K,Y1=Y1,D1=DN1,Y0=Y0,D0=DN0)
wK.rg<-temp$wK
sig2.rg<-temp$sig2K
Hrg.1<-temp$H1
Hrg.0<-temp$H0
ZK.rg<-c(0,temp$ZK)
rm("temp")

# ZK.rg(.) is the wLR Z statistic over time
# Next, calculate the summary statistic
# E.g., if Zstat.type="overall" then overall (final) Z;
# If Zstat.type="renyi" then max{|Z(.)}}

# NOTE: wE CAN CALCULATE Both "overall" and "Renyi" without much additional
# computations ... but focus on 1 version first.

Z.rg.obs.vec[rg]<-Z.stat(ZK.rg,stat="overall")
Z.rg2.obs.vec[rg]<-Z.stat(ZK.rg,stat="renyi")
}

# Note: Can also calculate RMST as function of truncation
# and develop Renyi- sup over time (truncation) version
# But for now, just take at specified L
if(rhogamma[1]==-1 & rhogamma[2]==-1){
# RMST   
Z.rg.obs.vec[rg]<-Z.rmst
Z.rg2.obs.vec[rg]<--99999 
}

if(draws>0){
for(dd in 1:draws){

if(rhogamma[1]>=0){  
g0<-G0.draws[,dd]
counting<-unlist(lapply(as.list(jump.points),N.weighted,error=U0[delta0==1],w=g0[delta0==1]))
counting<-c(0, counting)
DN0.star<- diff(counting)

g1<-G1.draws[,dd]
counting<-unlist(lapply(as.list(jump.points),N.weighted,error=U1[delta1==1],w=g1[delta1==1]))
counting<-c(0, counting)
DN1.star<- diff(counting)
d.perturb<-Hrg.0*DN0.star-Hrg.1*DN1.star
Z.rg.star<-ifelse(sig2.rg>0,c(cumsum(d.perturb[rbar>0]))/sqrt(sig2.rg),0.0)
ZK.rg.star<-c(0,Z.rg.star)
Z.rg.star.mat[rg,dd]<-Z.stat(ZK.rg.star,stat="overall")
Z.rg2.star.mat[rg,dd]<-Z.stat(ZK.rg.star,stat="renyi")
}
  
if(rhogamma[1]==-1 & rhogamma[2]==-1){
Z.rg.star.mat[rg,dd]<-c(Z.rmst.star[dd])
Z.rg2.star.mat[rg,dd]<--9999
}

# Show status for last rg=m weighting
if(rg==m & draws.status){
if(round(dd/draws,digits=3)==0.10) cat("10% resampling done","\n")
if(round(dd/draws,digits=3)==0.25) cat("25% resampling done","\n")
if(round(dd/draws,digits=3)==0.75) cat("75% resampling done","\n")
if(round(dd/draws,digits=3)==0.90) cat("90% resampling done","\n")
}
}  # Finish resampling for rg combination m
   }  # Finish m
#rg.max.star<-rho.gamma.mat[which.max(Z.rg.star.mat[,dd]),]
#rgs.mat[,3]<-rgs.mat[,3]+apply(rgs.mat,1,rg.count.star,rg0=rg.combo)
} 
# we have (Z(rho1,gamma1),Z(rho2,gamma2),...,Z(rhom,gammam)) stored in Z.rg.obs.vec
# Max statistic is then
Z.combo<-max(Z.rg.obs.vec)
Z.combo.renyi<-max(Z.rg2.obs.vec)

# The (rho,gamma) combination that corresponds to the maximum 
# This will be used to calculate the corresponding weighted Cox estimator

rg.locs<-which(rho.gamma.mat[,1]>=0)
rho.gamma.mat2<-rho.gamma.mat[rg.locs,]

rg.max<-rho.gamma.mat2[which.max(Z.rg.obs.vec[rg.locs]),]
rg.max.renyi<-rho.gamma.mat2[which.max(Z.rg2.obs.vec[rg.locs]),]

# Maximizer across all (r,g)'s and RMST 
max.all<-rho.gamma.mat[which.max(Z.rg.obs.vec),]
max.renyi.all<-rho.gamma.mat[which.max(Z.rg2.obs.vec),]

# Each column of Z.rg.star.mat contains a resample for Z.rg
# Calculate max across columns to get resample realization of Z.combo
Z.combo.star<-apply(Z.rg.star.mat,2,max)

Z.combo.renyi.star<-apply(Z.rg2.star.mat,2,max)

if(length(Z.combo.star)!=draws) stop("Z* not of right dimension")

pval.combo<-mean(c(Z.combo.star>=Z.combo),na.rm=TRUE)
pval.renyi<-mean(c(Z.combo.renyi.star>=Z.combo.renyi),na.rm=TRUE)


# Individual 
# Standard 1-sided (superiority) Pvalue
# Note: Aymptotics only available [here] for Zstat.type="overall"
pval.Z.asy<-c(1-pnorm(Z.rg.obs.vec))

# Un-adjusted p-value (will adjust when calculating rejections)
pz.combo.naive<-c(1-pnorm(Z.combo))

if(draws==0) pval.combo<-pz.combo.naive

# SM denotes 'synthetic-martingale' resampling
# ... to differentiate with asymptotic for overall Zstat
Ztemp<-cbind(Z.rg.obs.vec,Z.rg.star.mat)
# The first row is the observed test Ztemp[1,]
# Ztemp[-1,] will be the resamples
# pvalZ.SM calculates mean(Ztemp[1]>=Ztemp[-1])
# = prop(Z*>=Zobs)
pval.Z.SM<-apply(Ztemp,1,pvalZ.SM)

if(checking){ 
# "overall" stat
for(rg in 1:m){
if(details){
cat("rho,gamma=",rho.gamma.mat[rg,],"\n")
cat("Standard weighted Log-Rank: Z,Z^2",c(Z.rg.obs.vec[rg],Z.rg.obs.vec[rg]^2),"\n")
cat("Standard wLR 2-sided Pvalue",c(1-pchisq(Z.rg.obs.vec[rg]^2,1)),"\n")
  }
if(rho.gamma.mat[rg,1]==0 & rho.gamma.mat[rg,2]==0){
checkit<-survdiff(Surv(time,delta)~z)
if(details) print(checkit)
#print(cbind(checkit$chisq,Z.rg.obs.vec[rg]^2))
if(!identical(round(checkit$chisq,6),round(Z.rg.obs.vec[rg]^2,6))) stop("Problem with Z(rho,gamma) calculation")
}
if(rho.gamma.mat[rg,1]>0 & rho.gamma.mat[rg,2]==0){
checkit<-survdiff(Surv(time,delta)~z,rho=rho.gamma.mat[rg,1])
if(details) print(checkit)
if(!identical(round(checkit$chisq,6),round(Z.rg.obs.vec[rg]^2,6))) stop("Problem with Z(rho,gamma) calculation")
}
 }
}


if(m>1){
stats.Z<-cbind(rho.gamma.mat,Z.rg.obs.vec,pval.Z.SM,pval.Z.asy)
out<-list(stats.Z=stats.Z,rg.max=rg.max,Z.combo=Z.combo,pval.combo=pval.combo,rg.max.renyi=rg.max.renyi,pval.renyi=pval.renyi,Z.combo.renyi=Z.combo.renyi,
          G1.draws=G1.draws,G0.draws=G0.draws,pz.combo.naive=pz.combo.naive,max.all=max.all,max.renyi.all=max.renyi.all,fit.RMST=fit.RMST)
colnames(out$stats.Z) = c("rho","gamma","stat","pval","pval.asy")
if(draws.out) out$Z.combo.star<-Z.combo.star
if(draws>=20){# Output first 20 re-samples of Z.rg.star.mat
out$Z.rg.SM<-cbind(rho.gamma.mat,Z.rg.star.mat[,c(1:20)])
colnames(out$Z.rg.SM)=c("rho","gamma",paste("Z.SM",c(1:20),sep="."))
     }
  }

if(m==1){
stats.Z<-c(rho.gamma.mat,Z.rg.obs.vec,pval.Z.SM,pval.Z.asy)
out<-list(stats.Z=stats.Z,rg.max=rg.max,Z.combo=Z.combo,pval.combo=pval.combo,rg.max.renyi=rg.max.renyi,pval.renyi=pval.renyi,Z.combo.renyi=Z.combo.renyi,
          G1.draws=G1.draws,G0.draws=G0.draws,pz.combo.naive=pz.combo.naive)
names(out$stats.Z) = c("rho","gamma","stat","pval","pval.asy")
if(draws.out) out$Z.combo.star<-Z.combo.star
if(draws>=20){# Output first 20 re-samples of Z.rg.star.mat
out$Z.rg.SM<-c(rho.gamma.mat,Z.rg.star.mat[c(1:20)])
names(out$Z.rg.SM)=c("rho","gamma",paste("Z.SM",c(1:20),sep="."))
     }
  }
out
}
