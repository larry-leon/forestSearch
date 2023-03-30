#requires(data.table)

# minp = Minimum prevalence rate 
# Only combinations where all factors have prevalance at least minp are evaluated
# Minimum difference in subgroup sample size
# Max subgroup size is n: After first factor x1 is added, the sample size for 
# the subgroup is then n(x1), say.   If combining with next factor x2 
# does not reduce the sample size by rmin, then consider this combination "redundant"

# Adding d.min to require min(Events) per treatment arm (d0.min for control, d1.min for treatment)
subgroup.search<-function(Y,Event,Treat,ID,Z,n.min=30,d0.min=15,d1.min=15,hr.threshold=1.0,
                          max.minutes=30,time_index=3,
                          minp=0.025,rmin=5,details=FALSE,maxlimit=2^30.75,maxk=10){

 if (!("data.table" %in% utils::installed.packages())) stop("Package data.table required")  
  
  non.redundant<-0.0
  ngroups.fit<-0.0
  ngroups.meet<-0.0
  
  # output max of sg estimates
  # initialize
  max_sg_est<--Inf
  
  id.group<-0.0 # Identify subgroup (Will output subgroups directly)
  subs.out<-NULL
  
  temp<-cbind(Y,Event,Treat,ID,Z)
  temp<-na.exclude(temp)
  yy<-temp[,1]
  dd<-temp[,2]
  tt<-temp[,3]
  pid<-temp[,4]
  zz<-temp[,-c(1,2,3,4)]
  
  n<-length(yy)
  k<-ncol(zz)
  covs.in<-as.matrix(rep(0,k)) # Initially no members
  
  # Store HRs Covs.in
  HR.model.k<-NULL
  
  limit<-(2**k)-1 # Possible configurations including redundancies 
  
  if(details){
    cat("Number of possible subgroups=",c(limit),"\n")
    if(limit>10^6) cat("Number of possible subgroups (in millions)=",c(limit/(10^6)),"\n")
  }
  
  
  # Reset limit to at most maxlimit
  limit<-min(limit,maxlimit)
  
  # Track criteria not met for subgroup evaluation
  crit.fail0<-rep(0,limit)
  # crit.fail0 is first criteria to be met
  # The following criteria are conditional on meeting crit.fail0
  crit.failure<-rep(0,limit)
  crit.faild0<-rep(0,limit)
  crit.faild1<-rep(0,limit)
  
  t.start<-proc.time()[time_index]
  
  t.sofar<-0
  
  for(kk in 1:limit){
    
    t.now<-proc.time()[time_index]
    t.sofar<-c((t.now-t.start)/60)
    
    # Stop search algorithm if time exceeded: Output accumulated results
    if(t.sofar>max.minutes){ 
      out<-NULL
      if(!is.null(HR.model.k)){
        hr.out<-data.table(HR.model.k)
        names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
        rownames(hr.out)<-NULL
        hr.out<-setorder(hr.out,-HR,K)
        out<-list(ngroups.fit=ngroups.fit,hr.subgroups=hr.out,
                  time_search=t.sofar,
                  max_sg_est=c(hr.out[1,"HR"]))
      }
      if(details){ 
        cat("# of subgroups=",c(non.redundant),"\n")
        cat("# of subgroups fitted=",c(ngroups.fit),"\n")
        cat("Minutes=",c(t.sofar),"\n")
      }    
      return(out)
      
      stop("Maximum timing exceeded")
    }
    
    isin<-ztrail(kk)
    covs.in[isin,]<-one.zero(covs.in[isin,])
    
    ## We now have subgroup indices
    ## So extract subgroup membership
    k.in<-sum(covs.in)  
    
    crit0<-(k.in<=maxk)
    
    crit.fail0[kk]<-(k.in>maxk)
    
    if(crit0){
      
      x<-matrix(NA,n,k.in)
      index<-1
      
      for(jj in 1:k){
        if(covs.in[jj]==1){
          x[,index]<-zz[,jj]
          index<-index+1
        } 
      }
      # If any two factors are mutually exclusive then set id.x=NULL (-> next iteration)  
      xpx<-t(x)%*%x
      if(!all(xpx>0)) id.x<-NULL
      
      if(all(xpx>0)){
        id.x<-rep(1,nrow(temp))
        # Elements will be 1 if in subgroup defined by x  
        flag.redundant<-FALSE
        nx.last<-n # To flag redundancy
        for(m in 1:ncol(x)){
          if(!flag.redundant){
            id.x<-c(id.x)*c(x[,m])
            nx.this<-sum(id.x)
            if(nx.last-nx.this<=rmin) flag.redundant<-TRUE
            nx.last<-nx.this
          }
        }
      }  
      
      # Prevalance of each factor
      #prop.x<-apply(x,2,mean)    
      # Use flag.redundant 
      
      d1<-sum(c(dd*tt)[which(id.x==1)]) # number of events for treatment
      d0<-sum(c(dd*(1-tt))[which(id.x==1)]) # number of events for control
      
      # non-redundant flags additional criteria (in addition to flag.redundant above)  
      
      crit1<-(sum(id.x)>n.min)
      crit2<-(!flag.redundant)
      crit3<-(all(apply(x,2,mean)>=minp))
      crit4<-(d1>=d1.min)
      crit5<-(d0>=d0.min)
      
      crit.faild1[kk]<-c(d1<d1.min)
      crit.faild0[kk]<-c(d0<d0.min)
      
      crit.sum<-sum(c(crit1,crit2,crit3,crit4,crit5))
      if(crit.sum<5) crit.failure[kk]<-crit.sum
      
      if(crit.sum==5){
        data.x<-data.frame(Y=yy,E=dd,Treat=tt,ID=pid,id.x=id.x)
        df.x<-subset(data.x,id.x==1)
        #cat("Subgroup size=",c(dim(df.x)),"\n")
        
        non.redundant<-non.redundant+1
        
        # HR and median for subgroup
        # Try to speedup by only calculating HR's withOut CIs 
     
        #bhat.cox<-try(get.Cox(time=df.x$Y,delta=df.x$E,z=df.x$Treat)$bhat,TRUE)
       
        hr.cox<-try(summary(coxph(Surv(Y,E)~Treat,data=df.x,robust=FALSE))$conf.int,TRUE)
       
        if(!inherits(hr.cox,"try-error")){
          
          max_sg_est<-max(c(max_sg_est,hr.cox[1]))
          
          ngroups.fit<-ngroups.fit+1
          
          if(hr.cox[1]>hr.threshold){
   
            ngroups.meet<-ngroups.meet+1
            # Experimental median
            #fit1<-survfit(Surv(Y,E)~1,data=subset(df.x,Treat==1))
            #m1<-summary(fit1)$table["median"]
            
            df1.x<-df.x[which(df.x$Treat==1),]
            fit1<-NA.CHR.Weighted(time=df1.x$Y,Delta=(df1.x$E==1),at.points=sort(df1.x$Y[which(df1.x$E==1)]))
            S1.KM<-fit1$S.KM
            at.points1<-fit1$at.points
            m1<-suppressWarnings(min(at.points1[S1.KM<=0.50]))
            # SAS definition
            #m1a<-suppressWarnings(min(at.points1[S1.KM<=0.50]))
            #m1b<-suppressWarnings(max(at.points1[S1.KM>=0.50]))
            #m1<-(m1a+m1b)/2
            
            # Add m0
            df0.x<-df.x[which(df.x$Treat==0),]
            fit0<-NA.CHR.Weighted(time=df0.x$Y,Delta=(df0.x$E==1),at.points=sort(df0.x$Y[which(df0.x$E==1)]))
            S0.KM<-fit0$S.KM
            at.points0<-fit0$at.points
            m0<-suppressWarnings(min(at.points0[S0.KM<=0.50]))
            rm("df0.x","fit0")
            
            id.group<-id.group+1
            
            hr.this<-c(id.group,sum(covs.in),nrow(df.x),sum(df.x$E),sum(df1.x$E),m1,m0,hr.cox[c(1,3,4)],covs.in)
            
            rm("df1.x","fit1")
            df.group<-df.x
            df.group$id.group<-id.group
            
            HR.model.k<-rbind(HR.model.k,hr.this)
          } # Estimation loop
          
          if(details){
          if(t.sofar>=0.249 & t.sofar<=0.25){
            cat("Minutes so far=",c(t.sofar),"\n")
          }
          if(t.sofar==0.50){
            cat("Minutes so far=",c(t.sofar),"\n")
          }
          if(t.sofar==0.75){
            cat("Minutes so far=",c(t.sofar),"\n")
          }  
          if(t.sofar==1.00){
            cat("Minutes so far=",c(t.sofar),"\n")
          }
              }
        }
      }
    }
  } # maxk 
  
  t.end<-proc.time()[time_index]
  t.min<-(t.end-t.start)/60
  
  if(details){
    cat("# of subgroups based on # variables > k.max and excluded",c(sum(crit.fail0)),"\n")
    cat("k.max=",c(maxk),"\n")
    cat("Events criteria for control,exp=",c(d0.min,d1.min),"\n")
    cat("# of subgroups with events less than criteria: control, experimental",c(sum(crit.faild0),sum(crit.faild1)),"\n")
    cat("# of subgroups meeting all criteria =",c(non.redundant),"\n")
    cat("# of subgroups fitted (Cox model estimable) =",c(ngroups.fit),"\n")
    cat("Minutes=",c(t.min),"\n")
    cat("Number of criteria not met for subgroup evaluation","\n")
    print(table(crit.failure))
    cat("Number of subgroups meeting HR threshold",c(ngroups.meet),"\n")
    }
  
  out.found<-NULL
  if(!is.null(HR.model.k)){
    hr.out<-data.table(HR.model.k)
    names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
    rownames(hr.out)<-NULL
    hr.out<-setorder(hr.out,-HR,K)
    
    out.found<-list(ngroups.fit=ngroups.fit,hr.subgroups=hr.out)
  }
  if(is.null(out.found) & details) cat("No subgroups found","\n")
  
  return(list(out.found=out.found,
              max_sg_est=max_sg_est,
              time_search=t.sofar))
}
