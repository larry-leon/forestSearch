# minp = Minimum prevalence rate 
# Only combinations where all factors have prevalance at least minp are evaluated
# Minimum difference in subgroup sample size
# Max subgroup size is n: After first factor x1 is added, the sample size for 
# the subgroup is then n(x1), say.   If combining with next factor x2 
# does not reduce the sample size by rmin, then consider this combination "redundant"
# Adding d.min to require min(Events) per treatment arm (d0.min for control, d1.min for treatment)
# covs.in[isin,]<-unlist(lapply(covs.in[isin,],function(x){ifelse(x==1,0,1)})) # slower

# Removing ID storage since not needed
subgroup.search<-function(Y,Event,Treat,ID=NULL,Z,n.min=30,d0.min=15,d1.min=15,
hr.threshold=1.0,max.minutes=30,minp=0.05,rmin=5,details=FALSE,maxk=3){

non.redundant<-0.0
ngroups.fit<-0.0
ngroups.meet<-0.0
  
  # output max of sg estimates
  # initialize
  max_sg_est<--Inf
  
  id.group<-0.0 # Identify subgroup (Will output subgroups directly)
  subs.out<-NULL
  
  temp<-cbind(Y,Event,Treat,Z)
  temp<-na.exclude(temp)
  yy<-temp[,1]
  dd<-temp[,2]
  tt<-temp[,3]
  zz<-temp[,-c(1,2,3)]
  
  n<-length(yy)
  L<-ncol(zz)
  covs.in<-as.matrix(rep(0,L)) # Initially no members
  
  # Store results appending to HR.model.k
  HR.model.k<-NULL
  
  limit<-(2**L)-1 # Possible configurations including redundancies 
  
  if(details){
    cat("Number of unique levels (L) and possible subgroups=",c(L,limit),"\n")
    if(limit>10^6) cat("Number of possible subgroups (in millions)=",c(limit/(10^6)),"\n")
  }
  
  # Track criteria not met for subgroup evaluation
  crit.fail0<-rep(0,limit)
  # crit.fail0 is first criteria to be met
  # The following criteria are conditional on meeting crit.fail0
  crit.failure<-rep(0,limit)
  crit.faild0<-rep(0,limit)
  crit.faild1<-rep(0,limit)
  crit.failn<-rep(0,limit)
  
  t.start<-proc.time()[3]
  t.sofar<-0
  
  for(kk in 1:limit){
    
    t.now<-proc.time()[3]
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
      
      for(jj in 1:L){
        if(covs.in[jj]==1){
          x[,index]<-zz[,jj]
          index<-index+1
        } 
      }
      # If any two factors are mutually exclusive then set id.x=NULL (-> next iteration)  
      xpx<-t(x)%*%x
      if(!all(xpx>0)) id.x<-NULL
      
      if(all(xpx>0)){
        id.x<-rep(1,n)
        # Elements will be 1 if in subgroup defined by x  
        flag.redundant<-FALSE
        nx.last <- n # To flag redundancy
        for(m in 1:ncol(x)){
          if(!flag.redundant){
           id.x <- c(id.x)*c(x[,m])
            nx.this <- sum(id.x)
           if(nx.last-nx.this <= rmin) flag.redundant<-TRUE
            nx.last <- nx.this
         }
        } # flag redundant
      }  
      
      d1 <- sum(c(dd*tt)[which(id.x==1)]) # number of events for treatment
      d0 <- sum(c(dd*(1-tt))[which(id.x==1)]) # number of events for control
      nx<-sum(id.x)
      crit.failn[kk]<-(nx<=n.min)
      crit.faild1[kk]<-c(d1<d1.min)
      crit.faild0[kk]<-c(d0<d0.min)
      crit3<-(all(apply(x,2,mean)>=minp))
      
      if(!flag.redundant & d1>=d1.min & d0>=d0.min & nx>n.min & crit3){
      data.x<-data.table(Y=yy,E=dd,Treat=tt,id.x=id.x)
        
       df.x<-subset(data.x,id.x==1)
       non.redundant<-non.redundant+1
    
        hr.cox<-try(summary(coxph(Surv(Y,E)~Treat,data=df.x,robust=FALSE))$conf.int,TRUE)
       
        if(!inherits(hr.cox,"try-error")){
         max_sg_est<-max(c(max_sg_est,hr.cox[1]))
         ngroups.fit<-ngroups.fit+1
          
          if(hr.cox[1]>hr.threshold){
            ngroups.meet<-ngroups.meet+1
            meds<-summary(survfit(Surv(Y,E)~Treat,data=df.x))$table[,"median"]
            m0<-meds[1]
            m1<-meds[2]
            id.group<-id.group+1
            # ANY additions or modifications needs to update consistency function
            hr.this<-c(id.group,sum(covs.in),nx,d0+d1,d1,m1,m0,hr.cox[c(1,3,4)],covs.in)
            HR.model.k<-rbind(HR.model.k,hr.this)
          } # Estimation loop
        }
      }
  } # maxk 
  } # d1, d0, nmin | maxk criteria
  t.end<-proc.time()[3]
  t.min<-(t.end-t.start)/60
  
  if(details){
    cat("# of subgroups based on # variables > k.max and excluded (per million)",c(sum(crit.fail0)/10^6),"\n")
    cat("k.max=",c(maxk),"\n")
    cat("Events criteria for control,exp=",c(d0.min,d1.min),"\n")
    cat("# of subgroups with events less than criteria: control, experimental",c(sum(crit.faild0),sum(crit.faild1)),"\n")
    cat("# of subgroups with sample size less than criteria",c(sum(crit.failn)),"\n")
    cat("# of subgroups meeting all criteria =",c(non.redundant),"\n")
    cat("# of subgroups fitted (Cox model estimable) =",c(ngroups.fit),"\n")
    cat("*Subgroup Searching Minutes=*",c(t.min),"\n")
    cat("Number of subgroups meeting HR threshold",c(ngroups.meet),"\n")
    }
  
  out.found<-NULL
  if(!is.null(HR.model.k)){
    hr.out<-data.table(HR.model.k)
    names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
    # ANY additions or modifications needs to update consistency function
    rownames(hr.out)<-NULL
    hr.out<-setorder(hr.out,-HR,K)
    
    out.found<-list(ngroups.fit=ngroups.fit,hr.subgroups=hr.out)
  }
  
  if(is.null(out.found) & details){
  cat_rule(col="red")
  cat_line("NO subgroup candidate found (FS)",col="red")
  cat_rule(col="red")
  }    
  if(!is.null(out.found) & details){
    cat_rule(col="green")
    cat_line("Subgroup candidate(s) found (FS)",col="blue")
    cat_rule(col="green")
  }    
  
return(list(out.found=out.found,max_sg_est=max_sg_est,time_search=t.sofar,L=L))
}
