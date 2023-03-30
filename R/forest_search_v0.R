# Note: At some point may be interesting to compare to bestglm() 
# from bestglm package.  E.g., consider poisson regression 
# with counts as events and time/followup as offset

# Necessary packages

requiredPackages <- c("grf","policytree","data.table","randomForest","survival","plyr")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    install.packages(new.pkg, dependencies = TRUE, quietly=TRUE, verbose=FALSE, logical.return=FALSE)
  sapply(pkg, require, character.only = TRUE, quietly=TRUE)
  }
}

ipak(requiredPackages)

# Within a single function
forestsearch<-function(df,confounders.name,outcome.name,event.name,treat.name,id.name,df.predict=NULL,n.min=30,hr.threshold=1.0,hr.consistency=1.5,
                       fs.splits=100,m1.threshold=Inf,stop.threshold=0.975,pconsistency.threshold,d0.min=5,d1.min=5,
                       split_method="Random",sg_focus="hr",pstop_futile=0.5,
                       max.minutes=10,minp=0.025,details=FALSE,maxlimit=2^30.75,maxk=7,
                       by.risk=NULL,plot.sg=FALSE,vi.grf.min=0.2){
  
  # Override for Nsg since will sort candidates by sg size (Top of list may have low consistency)
  if(sg_focus=="Nsg"){
  pstop_futile<-0.0  
  stop.threshold<-pconsistency.threshold
  # If want higher stop.threshold, then set pconsistency.thresholds so
  }
  
  if(!(split_method %in% c("SPlit","Random"))) stop("split_method must be either SPlit or Random")
  if(split_method=="SPlit") require(SPlit)
  
  if(!is.null(vi.grf.min)){
  if (!("grf" %in% utils::installed.packages())) stop("Package grf required for screening and ordering of candidate factors")  
  
  # Use GRF for screening and order
  # Covariates need to be converted to numeric scale
  # Note: even if dimension not reduced, the sorting is beneficial
  X<-as.matrix(df[,confounders.name])
  # Convert to numeric
  X<-apply(X,2,as.numeric)
  Y<-df[,outcome.name]
  W<-df[,treat.name]
  D<-df[,event.name]
  tau.rmst<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))
  
  # For screening we take 0.9*tau.rms
  cs.forest<-try(suppressWarnings(causal_survival_forest(X,Y,W,D,horizon=0.9*tau.rmst)),TRUE)
  
  vi.cs<-variable_importance(cs.forest)
  vi.cs2<-data.frame(confounders.name,vi.cs)
  vi.order<-order(vi.cs,decreasing=TRUE)
  vi.cs2<-vi.cs2[vi.order,]
  
  conf.screen<-vi.cs2[,1]
  selected.vars <- which(vi.cs2[,2] / mean(vi.cs2[,2]) > vi.grf.min)
  conf.screen<-conf.screen[selected.vars]
  
  if(details) cat("Confounders per grf screening",conf.screen,"\n")
  # Returned as confounders.evaluated
  }
  
  # Otherwise, do not screen and use initial factors
  if(is.null(vi.grf.min)) conf.screen<-confounders.name
  
  df.confounders<-df[,conf.screen]
  df.confounders<-dummy(df.confounders)
  
  Event<-D
  Treat<-W
  ID<-df[,c(id.name)]
  
  df.fs<-data.frame(Y,Event,Treat,ID,df.confounders)
  
  Z<-as.matrix(df.confounders)
  colnames(Z)<-names(df.confounders)
  
  
  find.grps<-subgroup.search(Y=Y,Event=Event,Treat=Treat,ID=ID,Z=Z,d0.min=d0.min,d1.min=d1.min,n.min=n.min,
                             hr.threshold=hr.threshold,max.minutes=max.minutes,details=details,maxk=maxk)
  
  sg.harm<-NULL; df.est<-NULL; df.pred.out<-NULL; grp.consistency<-NULL
  
  max_sg_est<-find.grps$max_sg_est
  
  if(!is.null(find.grps$out.found) & any(find.grps$out.found$hr.subgroups$HR>hr.consistency)){# Found something?

    if(plot.sg & is.null(by.risk)) by.risk<-round(max(Y)/12,0)
  
    grp.consistency<-subgroup.consistency(df=df.fs,hr.subgroups=find.grps$out.found$hr.subgroups,
    hr.threshold=hr.threshold,hr.consistency=hr.consistency,m1.threshold=m1.threshold,
    sg_focus=sg_focus,pstop_futile=pstop_futile,
    pconsistency.threshold=pconsistency.threshold,n.splits=fs.splits,maxk=maxk,
    stop.threshold=stop.threshold,details=details,by.risk=by.risk,plot.sg=plot.sg,split_method=split_method)
    
    if(!is.null(grp.consistency$result)){
      sg.harm<-grp.consistency$sg.harm
      
    # Dataset with treatment recommend flag
      df.est<-grp.consistency$data
      
      # Prediction dataset for cross-validation:
      # df.train is training dataset for estimation, df.pred=df.test 
      # Or just original analysis dataset which facilitates summaries based on 
      # source factors and comparison with other methods  
      
      df.pred.out<-NULL
      if(!is.null(df.predict)){
        df.pred<-dummy(df.predict)
        df.pred$treat.recommend<-NA
        
        id.harm<-c(paste(sg.harm, collapse="==1 & "))
        id.harm<-paste(id.harm,"==1")
        df.pred.0<-subset(df.pred,eval(parse(text=id.harm)))
        df.pred.0$treat.recommend<-0
        
        id.noharm<-c(paste(sg.harm, collapse="!=1 | "))
        id.noharm<-paste(id.noharm,"!=1")
        df.pred.1<-subset(df.pred,eval(parse(text=id.noharm)))
        df.pred.1$treat.recommend<-1
        df.pred.out<-data.frame(rbind(df.pred.0,df.pred.1))
      }
  
    }
  } # Found something 
  return(list(sg.harm=sg.harm,df.est=df.est,df.pred=df.pred.out,grp.consistency=grp.consistency,
find.grps=find.grps,confounders.evaluated=conf.screen,max_sg_est=max_sg_est,
pstop_futile=pstop_futile,sg_focus=sg_focus))
}




fs.estimates.out<-function(df,fs.est=NULL,dgm=NULL,cox.formula.sim=NULL,cox.formula.adj.sim=NULL,analysis="FS"){
p.cens<-mean(1-df[,event.name])  

Y<-df[,outcome.name]
W<-df[,treat.name]
D<-df[,event.name]
taumax<-min(c(max(Y[W==1 & D==1]),max(Y[W==0 & D==1])))

if(!is.null(fs.est)){
grp.consistency<-fs.est$grp.consistency  
sg.harm.fs<-grp.consistency$sg.harm
}

if(is.null(fs.est)){
  grp.consistency<-NA  
  sg.harm.fs<-NULL
}

fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
hr.itt<-c(fit[1])
l.itt<-c(fit[3])
u.itt<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.adj.sim,data=df,robust=FALSE))$conf.int
hr.adj.itt<-c(fit[1,1])
l.adj.itt<-c(fit[1,3])
u.adj.itt<-c(fit[1,4])
rm("fit")


if(dgm$model=="null" & !is.null(dgm$fs.harm.true)) stop("For dgm model null fs.harm.true should be null")

# dgm$fs.harm.true=NULL --> H does not exist

# If H exists and found something
if(!is.null(dgm$fs.harm.true) & !is.null(sg.harm.fs)){
any.H<-1.0
dfout<-fs.est$df.pred
# PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
# Pred H and true H
aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
# Pred Hc and true H
bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
# Pred H and true Hc         
cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
# Pred Hc and true Hc         
dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))

size.H<-with(dfout,sum(treat.recommend==0))
size.Hc<-with(dfout,sum(treat.recommend==1))

# hr estimate with true H
fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
hr.H.true<-c(fit[1])
l.H.true<-c(fit[3])
u.H.true<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
hr.Hc.true<-c(fit[1])
l.Hc.true<-c(fit[3])
u.Hc.true<-c(fit[4])
rm("fit")

# hr estimate with estimated H
fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
hr.H.hat<-c(fit[1])
l.H.hat<-c(fit[3])
u.H.hat<-c(fit[4])
rm("fit")

fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
hr.Hc.hat<-c(fit[1])
l.Hc.hat<-c(fit[3])
u.Hc.hat<-c(fit[4])
rm("fit")

# Bias1 w.r.t hr.H.true
b1.H<-c(hr.H.hat-hr.H.true)
b2.H<-c(hr.H.hat-dgm$hr.H.true)

b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)

ppv<-aa/(aa+bb)
npv<-dd/(cc+dd)
specificity<-dd/(bb+dd)
sensitivity<-aa/(aa+cc)

#cat("PPV and NPV=",c(ppv,npv),"\n")
#cat("Sensitivity and Specificity=",c(specificity,sensitivity),"\n")

# Identify truth in dgm
found.1<-sum(grepl(dgm$fs.harm.true[1],sg.harm.fs))
found.2<-sum(grepl(dgm$fs.harm.true[2],sg.harm.fs))
found.both<-as.numeric(c(found.1+found.2)==2)
#found.al2<-as.numeric(length(sg.harm.fs)>=2)
found.al3<-as.numeric(length(sg.harm.fs)>=3)

}

# If H exists and did NOT find something
if(!is.null(dgm$fs.harm.true) & is.null(sg.harm.fs)){
  any.H<-0
  size.H<-0
  size.Hc<-nrow(df)
  
  # hr estimate with true H
#  hr.H.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int[1])
#  hr.Hc.true<-c(summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int[1])
  
  fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==1),robust=FALSE))$conf.int
  hr.H.true<-c(fit[1])
  l.H.true<-c(fit[3])
  u.H.true<-c(fit[4])
  rm("fit")
  
  fit<-summary(coxph(cox.formula.sim,data=subset(df,flag.harm==0),robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
 
  # hr estimate with estimated H, in this case, does not exist
  hr.H.hat<-NA
  l.H.hat<-NA
  u.H.hat<-NA
  
  # Define HC as ITT 
 # hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])

 fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
   # Bias1 w.r.t hr.H.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  # no H found, treat.recommend==1 for all
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # NO H found
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  # All treat.recommend=1
  bb<-with(df,sum(flag.harm==1))
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-0
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  dd<-with(df,sum(flag.harm==0))
  
  ppv<-aa/(aa+bb)
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-ifelse(aa==0,0,aa/(aa+cc))
  
  # Identify truth in dgm
  found.1<-0
  found.2<-0
  found.both<-0
  found.al3<-0
  
}


# If H does NOT exist and found something
if(is.null(dgm$fs.harm.true) & !is.null(sg.harm.fs)){
  any.H<-1.0
  dfout<-fs.est$df.pred
  # PPV for H membership: Hpred --> treat.reccomend=0 & flag.harm=1
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # flag.harm=0 for all
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  bb<-0
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-with(dfout,sum(treat.recommend==0))
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  dd<-with(dfout,sum(treat.recommend==1))
  
  size.H<-with(dfout,sum(treat.recommend==0))
  size.Hc<-with(dfout,sum(treat.recommend==1))
  
  # hr estimate with true H
  hr.H.true<-NA
  l.H.true<-NA
  u.H.true<-NA
  
  # no H --> ITT
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
  
  # hr estimate with estimated H
  fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==0),robust=FALSE))$conf.int
  hr.H.hat<-c(fit[1])
  l.H.hat<-c(fit[3])
  u.H.hat<-c(fit[4])
  rm("fit")
  
  fit<-summary(coxph(cox.formula.sim,data=subset(dfout,treat.recommend==1),robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
  # Bias1 w.r.t hr.H.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  ppv<-NA
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-aa/(aa+cc)
  
  # Identify truth in dgm
  # fs.harm.true will be null
  # Still want to track how often v1.1 and v4.1 are found
  found.1<-sum(grepl(c("v1.1"),sg.harm.fs))
  found.2<-sum(grepl(c("v4.1"),sg.harm.fs))
  found.both<-as.numeric(c(found.1+found.2)==2)
  found.al3<-as.numeric(length(sg.harm.fs)>=3)
  
}

# If H does NOT exist and did NOT find something
# Define hat(Hc)=ITT
if(is.null(dgm$fs.harm.true) & is.null(sg.harm.fs)){
  any.H<-0
  # **flag.harm=0 for all**
  # Pred H and true H
  #aa<-with(dfout,sum(treat.recommend==0 & flag.harm==1))
  # flag.harm=0 for all
  aa<-0
  # Pred Hc and true H
  #bb<-with(dfout,sum(treat.recommend==1 & flag.harm==1))
  bb<-0
  # Pred H and true Hc         
  #cc<-with(dfout,sum(treat.recommend==0 & flag.harm==0))
  cc<-0
  # Pred Hc and true Hc         
  #dd<-with(dfout,sum(treat.recommend==1 & flag.harm==0))
  # Did NOT find anything --> treat.recommend=1 for all
  dd<-nrow(df)
  
  size.H<-0
  size.Hc<-nrow(df)
  
  # hr estimate with true H
  hr.H.true<-NA
  l.H.true<-NA
  u.H.true<-NA
  
  # ITT
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.true<-c(fit[1])
  l.Hc.true<-c(fit[3])
  u.Hc.true<-c(fit[4])
  rm("fit")
  
  
  # hr estimate with estimated H
  hr.H.hat<-NA
  l.H.hat<-NA
  u.H.hat<-NA
  
  #hr.Hc.hat<-c(summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int[1])
  fit<-summary(coxph(cox.formula.sim,data=df,robust=FALSE))$conf.int
  hr.Hc.hat<-c(fit[1])
  l.Hc.hat<-c(fit[3])
  u.Hc.hat<-c(fit[4])
  rm("fit")
  
  # Bias1 w.r.t hr.H.true
  b1.H<-NA
  b2.H<-NA
  
  b1.Hc<-c(hr.Hc.hat-hr.Hc.true)
  b2.Hc<-c(hr.Hc.hat-dgm$hr.Hc.true)
  
  ppv<-NA
  npv<-dd/(cc+dd)
  specificity<-dd/(bb+dd)
  sensitivity<-NA
  
  # Identify truth in dgm
  # fs.harm.true will be null
  # Still want to track how often v1.1 and v4.1 are found
  found.1<-sum(grepl(c("v1.1"),sg.harm.fs))
  found.2<-sum(grepl(c("v4.1"),sg.harm.fs))
  found.both<-as.numeric(c(found.1+found.2)==2)
  found.al3<-as.numeric(length(sg.harm.fs)>=3)
  
}

df.res<-data.table(any.H,size.H,size.Hc,ppv,npv,specificity,sensitivity,found.1,found.2,found.both,found.al3,
hr.H.true,hr.Hc.true,hr.H.hat,hr.Hc.hat,b1.H,b2.H,b1.Hc,b2.Hc,p.cens,analysis,taumax,
hr.itt,l.itt,u.itt,
hr.adj.itt,l.adj.itt,u.adj.itt,
l.H.true,u.H.true,
l.Hc.true,u.Hc.true,
l.H.hat,u.H.hat,
l.Hc.hat,u.Hc.hat)
return(df.res)
}










