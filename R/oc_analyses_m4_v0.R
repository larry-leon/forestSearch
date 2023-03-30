# Default parameters
oc_analyses_m4_Nsg<-function(sim){
  
  ans.analyses<-NULL  
  
  x<-sim_aftm4_gbsg(dgm=dgm,n=N,maxFollow=maxFollow,muC.adj=muC.adj,simid=sim)
  
  # Initiate with # and prop H
  sizeH_true<-with(x,sum(flag.harm))
  propH_true<-with(x,mean(flag.harm))
  sizeHc_true<-with(x,sum(!flag.harm))
  propHc_true<-with(x,mean(!flag.harm))
  
  dfHpop<-data.table(sim,sizeH_true,propH_true,sizeHc_true,propHc_true)
  
  if(get.FS){
    ansFS<-NULL
    analysis<-label.analyses[1]
    
    y<-try(forestsearch(df=x,confounders.name=confounders.name,df.predict=x,details=FALSE,sg_focus="hr",
                        outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,id.name=id.name,
                        n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                        stop.threshold=0.95,pstop_futile=0.5,
                        d0.min=d.min,d1.min=d.min,
                        pconsistency.threshold=0.90,max.minutes=max.minutes,maxk=maxk,
                        plot.sg=FALSE,vi.grf.min=vi.grf.min,split_method=split_method),TRUE)
    
    if(!inherits(y,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.est=y,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
    }
    
    if(inherits(y,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.est=NULL,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
      if(round(c(ansFS$hr.Hc.hat-ansFS$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (FS)")
    }
  }
  
  if(get.FS){
    ansFS<-NULL
    analysis<-"FS-M"
        y<-try(forestsearch(df=x,confounders.name=confounders.name,df.predict=x,details=FALSE,sg_focus="Nsg",
                        outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,id.name=id.name,
                        n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                        stop.threshold=0.95,pstop_futile=0.0,
                        d0.min=d.min,d1.min=d.min,
                        pconsistency.threshold=0.95,max.minutes=max.minutes,maxk=maxk,
                        plot.sg=FALSE,vi.grf.min=vi.grf.min,split_method=split_method),TRUE)
    
    if(!inherits(y,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.est=y,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
    }
    
    if(inherits(y,"try-error")){
      ansFS<-fs.estimates.out(df=x,fs.est=NULL,dgm=dgm,
                              cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,
                              analysis=analysis)
      ansFS<-cbind(dfHpop,ansFS)
      ans.analyses<-rbind(ans.analyses,ansFS)
      if(round(c(ansFS$hr.Hc.hat-ansFS$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (FS)")
    }
  }
  
  # frac.tau<-1.0
  if(get.GRF){
    ansGRF<-NULL
    analysis<-label.analyses[2]
    
    grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
                                        event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
                                        dmin.grf=dmin.grf,frac.tau=1.0,details=FALSE),TRUE)
    
    if(!inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
    }
    
    if(inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
      if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
    }
    rm("grf.est")
    
  } # GRF analysis frac.tau=100%
  
  
  if(get.GRF){
    ansGRF<-NULL
    analysis<-label.analyses[7]
    
    grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
                                        event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
                                        dmin.grf=dmin.grf,frac.tau=frac.tau,details=FALSE),TRUE)
    
    if(!inherits(grf.est,"try-error")){
      
      ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=frac.tau)
      
      ansGRF<-cbind(dfHpop,ansGRF)
      
      ans.analyses<-rbind(ans.analyses,ansGRF)
    }
    
    if(inherits(grf.est,"try-error")){
      ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                                cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
      ansGRF<-cbind(dfHpop,ansGRF)
      ans.analyses<-rbind(ans.analyses,ansGRF)
      if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
    }
    rm("grf.est")
  } # GRF analysis frac.tau specified
  
  
  
  if(get.VT){
    # VT
    x2<-get.FG(x,tte.name="y.sim",event.name="event.sim")
    # VT at tau1, censored version 
    tau<-24
    analysis<-label.analyses[3]
    ansVT<-NULL
    vt.estA<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estA,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estA,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estA,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24)")
    }
    
    
    ansVT<-NULL
    analysis<-label.analyses[4]
    # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
    vt.estB<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estB,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estB,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estB,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24 ideal)")
    }
    
    
    
    
    tau<-36
    ansVT<-NULL
    analysis<-label.analyses[5]
    vt.estC<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    
    if(!inherits(vt.estC,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estC,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estC,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36)")
    }
    
    
    ansVT<-NULL
    analysis<-label.analyses[6]
    # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
    vt.estD<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                       seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
    if(!inherits(vt.estD,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=vt.estD,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
    }
    
    if(inherits(vt.estD,"try-error")){
      ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                              cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
      ansVT<-cbind(dfHpop,ansVT)
      
      ans.analyses<-rbind(ans.analyses,ansVT)
      if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36 ideal)")
    }
    
  }
  return(ans.analyses)
}



oc_analyses_m4<-function(sim){

ans.analyses<-NULL  

x<-sim_aftm4_gbsg(dgm=dgm,n=N,maxFollow=maxFollow,muC.adj=muC.adj,simid=sim)

# Initiate with # and prop H
sizeH_true<-with(x,sum(flag.harm))
propH_true<-with(x,mean(flag.harm))
sizeHc_true<-with(x,sum(!flag.harm))
propHc_true<-with(x,mean(!flag.harm))

dfHpop<-data.table(sim,sizeH_true,propH_true,sizeHc_true,propHc_true)

if(get.FS){
  ansFS<-NULL
  analysis<-label.analyses[1]
  
  y<-try(forestsearch(df=x,confounders.name=confounders.name,df.predict=x,details=FALSE,sg_focus=sg_focus,
                      outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,id.name=id.name,
                      n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                      stop.threshold=stop.threshold,
                      d0.min=d.min,d1.min=d.min,
                      pconsistency.threshold=pconsistency.threshold,max.minutes=max.minutes,maxk=maxk,
                      plot.sg=FALSE,vi.grf.min=vi.grf.min,split_method=split_method),TRUE)
  
  
  if(!inherits(y,"try-error")){
    ansFS<-fs.estimates.out(df=x,fs.est=y,dgm=dgm,
                          cox.formula.sim=cox.formula.sim,
                          cox.formula.adj.sim=cox.formula.adj.sim,
                          analysis=analysis)
    ansFS<-cbind(dfHpop,ansFS)
    ans.analyses<-rbind(ans.analyses,ansFS)
  }

  if(inherits(y,"try-error")){
    ansFS<-fs.estimates.out(df=x,fs.est=NULL,dgm=dgm,
                            cox.formula.sim=cox.formula.sim,
                            cox.formula.adj.sim=cox.formula.adj.sim,
                            analysis=analysis)
    ansFS<-cbind(dfHpop,ansFS)
    ans.analyses<-rbind(ans.analyses,ansFS)
    if(round(c(ansFS$hr.Hc.hat-ansFS$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (FS)")
  }
  
}

# frac.tau<-1.0
if(get.GRF){
ansGRF<-NULL
analysis<-label.analyses[2]

grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
dmin.grf=dmin.grf,frac.tau=1.0,details=FALSE),TRUE)

if(!inherits(grf.est,"try-error")){
  ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                         cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
  ansGRF<-cbind(dfHpop,ansGRF)
  ans.analyses<-rbind(ans.analyses,ansGRF)
}

if(inherits(grf.est,"try-error")){
  ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
  cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
  ansGRF<-cbind(dfHpop,ansGRF)
  ans.analyses<-rbind(ans.analyses,ansGRF)
  if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
}
rm("grf.est")

} # GRF analysis frac.tau=100%


if(get.GRF){
  ansGRF<-NULL
  analysis<-label.analyses[7]
  
  grf.est<-try(grf.subg.harm.survival(data=x,confounders.name=confounders.name,outcome.name=outcome.name,
                                      event.name=event.name,id.name=id.name,treat.name=treat.name,n.min=n.min,
                                      dmin.grf=dmin.grf,frac.tau=frac.tau,details=FALSE),TRUE)

  if(!inherits(grf.est,"try-error")){

    ansGRF<-grf.estimates.out(df=x,grf.est=grf.est,dgm=dgm,cox.formula.sim=cox.formula.sim,
                           cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=frac.tau)
  
      ansGRF<-cbind(dfHpop,ansGRF)
    
      ans.analyses<-rbind(ans.analyses,ansGRF)
  }
  
  if(inherits(grf.est,"try-error")){
    ansGRF<-grf.estimates.out(df=x,grf.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
    cox.formula.adj.sim=cox.formula.adj.sim,analysis=analysis,frac.tau=1.0)
    ansGRF<-cbind(dfHpop,ansGRF)
    ans.analyses<-rbind(ans.analyses,ansGRF)
    if(round(c(ansGRF$hr.Hc.hat-ansGRF$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (GRF)")
    }
rm("grf.est")
} # GRF analysis frac.tau specified



if(get.VT){
  # VT
  x2<-get.FG(x,tte.name="y.sim",event.name="event.sim")
  # VT at tau1, censored version 
  tau<-24
  analysis<-label.analyses[3]
  ansVT<-NULL
  vt.estA<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                     seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
  
  if(!inherits(vt.estA,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=vt.estA,dgm=dgm,cox.formula.sim=cox.formula.sim,
                          cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
  }

  if(inherits(vt.estA,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                            cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
    if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24)")
  }
  
  
    ansVT<-NULL
  analysis<-label.analyses[4]
  # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
  vt.estB<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                     seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
  
  if(!inherits(vt.estB,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=vt.estB,dgm=dgm,cox.formula.sim=cox.formula.sim,
                          cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
  }
  
  if(inherits(vt.estB,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                            cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
    if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 24 ideal)")
  }
  
  
  
  
  tau<-36
  ansVT<-NULL
  analysis<-label.analyses[5]
  vt.estC<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="Y.FG",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                     seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
  
  if(!inherits(vt.estC,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=vt.estC,dgm=dgm,cox.formula.sim=cox.formula.sim,
                          cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    ans.analyses<-rbind(ans.analyses,ansVT)
  }
  
  if(inherits(vt.estC,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                            cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
    if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36)")
  }
  
  
  ansVT<-NULL
  analysis<-label.analyses[6]
  # VT at tau1, NON-censored version t.sim (Note: event.YFG are 'dummy' event indicators == 1 for all)
  vt.estD<-try(vt.subg.harm.survival(data=x2,cf.names=confounders.name,tte.name="t.sim",event.name="event.YFG",treat.name="treat",ntree=ntree,
                                     seed=8316951,data.pred=NULL,details=FALSE,n.min=n.min,treat.threshold=treat.threshold,maxdepth=maxdepth,tau=tau),TRUE)
  if(!inherits(vt.estD,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=vt.estD,dgm=dgm,cox.formula.sim=cox.formula.sim,
                          cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
  }
  
  if(inherits(vt.estD,"try-error")){
    ansVT<-vt.estimates.out(df=x,vt.est=NULL,dgm=dgm,cox.formula.sim=cox.formula.sim,
                            cox.formula.adj.sim=cox.formula.adj.sim,vt.threshold=vt.threshold,analysis=analysis)
    ansVT<-cbind(dfHpop,ansVT)
    
    ans.analyses<-rbind(ans.analyses,ansVT)
    if(round(c(ansVT$hr.Hc.hat-ansVT$hr.itt),6)!=0) stop("Cox hr for H-complement NOT equal to ITT (VT 36 ideal)")
  }
  
}
return(ans.analyses)
}

tab_tests<-function(res,analyses=label.analyses){
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[1]))
  names(test)<-analyses[1]
  test.all<-test
  # initiate test.all 
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[2]))
  names(test)<-analyses[2]
  test.all<-cbind(test.all,test)
  
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[7]))
  names(test)<-analyses[7]
  test.all<-cbind(test.all,test)
  
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[3]))
  names(test)<-analyses[3]
  test.all<-cbind(test.all,test)
  
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[4]))
  names(test)<-analyses[4]
  test.all<-cbind(test.all,test)
  
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[5]))
  names(test)<-analyses[5]
  test.all<-cbind(test.all,test)
  
  test<-get.summary.VTFS(result=subset(res,analysis==analyses[6]))
  names(test)<-analyses[6]
  test.all<-cbind(test.all,test)
  return(test.all)
}



tab_tests_Nsg<-function(res,analyses=unique(res$analysis)){
test<-get.summary.VTFS(result=subset(res,analysis==analyses[1]))
names(test)<-analyses[1]
test.all<-test
  for(aa in 2:length(analyses)){
   test<-get.summary.VTFS(result=subset(res,analysis==analyses[aa]))
    names(test)<-analyses[aa]
    test.all<-cbind(test.all,test)
  }  
return(test.all)  
}


out.results<-function(res,dgm,output.file=NULL,t.min,details=TRUE){
if(details){
print(output.file)
print(head(res))
cat("Subgroup HRs: H, H^c, Causal=",c(dgm$hr.H.true,dgm$hr.Hc.true,dgm$hr.causal),"\n")
out.analysis<-subset(res,analysis=="FS")
cat("Simulations=",c(nrow(out.analysis)),"\n")
pC<-mean(out.analysis$p.cens)
cat("Avg censoring=",c(pC),"\n")
cat("Min,Max,Avg tau.max=",c(min(out.analysis$taumax),max(out.analysis$taumax),mean(out.analysis$taumax)),"\n")
}
t.hours<-t.min/60

if(!is.na(dgm$hr.H.true)){
hrH.plim<-mean(res$hr.H.true)
hrH.true<-dgm$hr.H.true

pAnyH.approx <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.plim,pC=pC,n.sg=60,k1=log(1.5),k2=log(1.25))$integral
if(details) cat("P(H) approximation at plim(H), approx=",c(hrH.plim,pAnyH.approx),"\n")

pAnyH.approx <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrH.true,pC=pC,n.sg=60,k1=log(1.5),k2=log(1.25))$integral
if(details) cat("P(H) approximation at hr(H) true, approx=",c(hrH.true,pAnyH.approx),"\n")
}

if(is.na(dgm$hr.H.true)){
hrHc.plim<-mean(res$hr.Hc.true)
hrHc.true<-dgm$hr.Hc.true

pAnyH.approx <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrHc.plim,pC=pC,n.sg=60,k1=log(1.5),k2=log(1.25))$integral
if(details) cat("P(H) approximation at plim(Hrc), approx=",c(hrHc.plim,pAnyH.approx),"\n")

pAnyH.approx <- adaptIntegrate(dens_threshold_both,lowerLimit=c(-Inf,-Inf),upperLimit=c(Inf,Inf),theta=hrHc.true,pC=pC,n.sg=60,k1=log(1.5),k2=log(1.25))$integral
if(details) cat("P(H) approximation at hr(Hc) true, approx=",c(hrHc.true,pAnyH.approx),"\n")

}


if(details) cat("Minutes,hours",c(t.min,t.hours),"\n")
# Output file
oc.tab<-tab_tests(res=res)
if(details) print(oc.tab)
if(!is.null(output.file)){
save(pAnyH.approx,t.hours,dgm,res,oc.tab,file=output.file)
}
}

get.summary.VTFS<-function(result,analysis="VT",dig1=2,dig2=3){
  # summarize finding H
  # remove "found.al3" from table
  dfH1a<-result[,c("any.H","found.1","found.2","found.both","ppv","npv","sensitivity","specificity")]  
  temp1a<-as.numeric(apply(dfH1a,2,mean,na.rm=TRUE))
  temp1a<-round(temp1a,digits=dig1)
  
  # Where hat(H) exists --> hat(Hc) exists
  resH<-subset(result,any.H==1)
  # resH includes all sequences where a subgroup is found
  
  #dfH1b<-resH[,c("sensitivity","specificity")]  
  dfH1b<-as.matrix(resH[,c("npv")])  
  temp1b<-as.numeric(apply(dfH1b,2,mean,na.rm=TRUE))
  temp1b<-round(temp1b,digits=dig1)
  
  temp1<-c(temp1a,temp1b)
  # Size of H
  avg.H<-mean(resH$size.H,na.rm=TRUE)
  min.H<-min(resH$size.H,na.rm=TRUE)
  max.H<-max(resH$size.H,na.rm=TRUE)
  # Hc
  avg.Hc<-mean(result$size.Hc,na.rm=TRUE)
  min.Hc<-min(result$size.Hc,na.rm=TRUE)
  max.Hc<-max(result$size.Hc,na.rm=TRUE)
  
  temp2<-round(c(avg.H,min.H,max.H,avg.Hc,min.Hc,max.Hc),digits=0)
  
  #dfH2<-result[,c("hr.H.true","hr.H.hat")]
  dfH2<-resH[,c("hr.H.true","hr.H.hat")]
  #dfH2<-subset(dfH2,any.H==1)
  temp3<-as.numeric(apply(dfH2,2,mean,na.rm=TRUE))
  temp3<-round(temp3,digits=dig2)
  
  dfH3<-resH[,c("hr.Hc.true","hr.Hc.hat")]  
  temp4<-as.numeric(apply(dfH3,2,mean,na.rm=TRUE))
  temp4<-round(temp4,digits=dig2)
  
  # Include hr.H.true and hr.Hc.true for all samples
  dfH4<-result[,c("hr.H.true","hr.Hc.true","hr.itt","hr.adj.itt")]  
  temp5<-as.numeric(apply(dfH4,2,mean,na.rm=TRUE))
  temp5<-round(temp5,digits=dig2)
  
  res.out<-c(temp1,temp2,temp3,temp4,temp5)
    # remove "found.al3"
  res.out2<-data.frame(res.out)
  
  # Rename ppv=ppH, npv=ppHc, sens=sensH, spec=sensHc
  row.names(res.out2)<-c("any.H","found.1","found.2","found.both","ppH","ppHc","sensH","sensHc","hatH_ppHc",
                         "Avg(#H)","minH","maxH","Avg(#Hc)","minHc","maxHc",
                         "hat(H*)","hat(hat[H])","hat(Hc*)","hat(hat[Hc])","hat(H*)all","hat(Hc*)all","hat(ITT)all","hat(ITTadj)all")
  return(res.out2)
}


hr_threshold_FS<-function(n.sg,pC=0.3,theta,hr.threshold,details=FALSE,alpha=0.025){
  d.sg<-n.sg*(1-pC)
  result<-1-pnorm(sqrt(d.sg/4)*(log(hr.threshold)-log(theta)-(2/sqrt(d.sg))*qnorm(1-alpha)))
 
  if(details) cat("Approximation=",c(result),"\n")
  return(result)
}



dens_threshold_both<-function(x,theta,pC=0.3,n.sg,k1,k2){
  d.sg<-n.sg*(1-pC)
  mu_split<-log(theta)
  sig2_split<-8/d.sg
  sig_split<-sqrt(sig2_split)
  dens_1<-dnorm(x[1],mean=mu_split,sd=sig_split)
  dens_2<-dnorm(x[2],mean=mu_split,sd=sig_split)
  ans<-c((x[1]+x[2]) >= 2*k1)*c(x[1]>=k2)*c(x[2]>=k2)*dens_1*dens_2
    return(ans)
}


# Functions used in bootstrap operating characteristics


# Replicates coxph
getci_Cox<-function(df,est,se,target,alpha=0.025,digits=3){
  a<-df[est]
  b<-df[se]
  if(!is.numeric(target)){
    c<-df[target]
  }
  if(is.numeric(target)){
    c<-target
  }
  log_lb<-log(a)-qnorm(0.975)*(b/a)
  log_ub<-log(a)+qnorm(0.975)*(b/a)
  lb<-exp(log_lb)
  ub<-exp(log_ub)
  cov<-ifelse(c>=lb & c<=ub,1,0)
  return(list(lb=round(lb,3),ub=round(ub,3),cover=cov))
}

# Include Avg lengths of CIs
# Include ppv sens
dfB_out<-function(df,dgm){
  df<-na.omit(df)
  target<-unlist(dgm["hr.H.true"])
  Hstats_a<-with(df,c(mean(H_obs),100*(mean(H_obs)-target)/target,sqrt(var(H_obs)),mad(H_obs),mean(seH_obs),min(H_obs),max(H_obs),mean(H0_length),c(mean(H0_cover1),mean(H0_cover2),mean(H0_cover3))))
  Hstats_b<-with(df,c(mean(H1.bc),100*(mean(H1.bc)-target)/target,sqrt(var(H1.bc)),mad(H1.bc),mean(seH1.bc),min(H1.bc),max(H1.bc),mean(H1_length),c(mean(H1_cover1),mean(H1_cover2),mean(H1_cover3))))
  # Remove H2.bc
  # Keep in case added in future
  #Hstats_c<-with(df,c(mean(H2.bc),median(H2.bc),sqrt(var(H2.bc)),mean(seH2.bc),min(H2.bc),max(H2.bc),100*c(mean(H2_cover1),mean(H2_cover2),mean(H2_cover3))))
  Hstats_c<-NULL
  # H_true
  temp<-getci_Cox(df=df,est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]))
  est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))
  Hstats_d<-with(df,c(mean(H_true),100*(mean(H_true)-target)/target,sqrt(var(H_true)),mad(H_true),mean(seH_true),min(H_true),max(H_true),mean(est_length),NA,NA,mean(temp$cover,na.rm=TRUE)))
  Hstats_e<-with(df,c(mean(hatH_causal),100*(mean(hatH_causal)-target)/target,sqrt(var(hatH_causal)),mad(hatH_causal),NA,min(hatH_causal),max(hatH_causal),NA,NA,NA,NA))
  Hstats<-rbind(Hstats_d,Hstats_e,Hstats_a,Hstats_b,Hstats_c)
  # Hc
  target<-unlist(dgm["hr.Hc.true"])
  Hcstats_a<-with(df,c(mean(Hc_obs),100*(mean(Hc_obs)-target)/target,sqrt(var(Hc_obs)),mad(Hc_obs),mean(seHc_obs),min(Hc_obs),max(Hc_obs),mean(Hc0_length),c(mean(Hc0_cover1),mean(Hc0_cover2),mean(Hc0_cover3))))
  Hcstats_b<-with(df,c(mean(Hc1.bc),100*(mean(Hc1.bc)-target)/target,sqrt(var(Hc1.bc)),mad(Hc1.bc),mean(seHc1.bc),min(Hc1.bc),max(Hc1.bc),mean(Hc1_length),c(mean(Hc1_cover1),mean(Hc1_cover2),mean(Hc1_cover3))))
  Hcstats_c<-NULL
  #Hcstats_c<-with(df,c(mean(Hc2.bc),median(Hc2.bc),sqrt(var(Hc2.bc)),mean(seHc2.bc),min(Hc2.bc),max(Hc2.bc),100*c(mean(Hc2_cover1),mean(Hc2_cover2),mean(Hc2_cover3))))
  # Hc_true
  temp<-getci_Cox(df=df,est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]))
  est_length<-c(unlist(temp$ub))-c(unlist(temp$lb))
  Hcstats_d<-with(df,c(mean(Hc_true),100*(mean(Hc_true)-target)/target,sqrt(var(Hc_true)),mad(Hc_true),mean(seHc_true),min(Hc_true),max(Hc_true),mean(est_length),NA,NA,mean(temp$cover,na.rm=TRUE)))
  Hcstats_e<-with(df,c(mean(hatHc_causal),100*(mean(hatHc_causal)-target)/target,sqrt(var(hatHc_causal)),mad(hatHc_causal),NA,min(hatHc_causal),max(hatHc_causal),NA,NA,NA,NA))
  Hcstats<-rbind(Hcstats_d,Hcstats_e,Hcstats_a,Hcstats_b,Hcstats_c)
  
  Tdf<-as.data.frame(rbind(Hstats,Hcstats))
  
  # Note: includes the causal potential-outcome "placeholder" (Set to null for now)
  rnH<-c("$\\hat\\theta(H)$","$\\theta^{\\dagger}(\\hat{H})$","$\\hat\\theta(\\hat{H})$","$\\hat\\theta^{*}(\\hat{H})$")
  rnHc<-c("$\\hat\\theta(H^{c})$","$\\theta^{\\dagger}(\\hat{H}^{c})$","$\\hat\\theta(\\hat{H}^{c})$","$\\hat\\theta^{*}(\\hat{H}^{c})$")
  rownames(Tdf)<-c(rnH,rnHc)
  # Remove rows 2,6
  # Removes the PO version for H and H^c
  Tdf_res<-Tdf[-c(2,6),]
  # Remove column 9
  # Removes Coverage of PO for H and H^{c}
  Tdf_res<-Tdf_res[,-c(9)]
  # New rownames
  colnames(Tdf_res)<-c("Avg","rBias","SD","MAD","est(SD)","min","max","Avg(L)","C(true)","C(fixed)")
  return(Tdf_res)
}


# Note: this can be simplified (see below)
get_df_bias_Boxplot<-function(res,dgm,target=c("H_known","H_causal","H_fixed")){
  # Boxplots H biases H_true, hatH_causal, dgm$hrH.true
  hrH_true<-dgm$hr.H.true
  hrHc_true<-dgm$hr.Hc.true
  
  res_new<-within(res,{
    b1H_1<-100*(H1.bc-H_true)/H_true
    b1H_2<-100*(H1.bc-hatH_causal)/hatH_causal
    b1H_3<-100*(H1.bc-hrH_true)/hrH_true
    
    b2H_1<-100*(H2.bc-H_true)/H_true
    b2H_2<-100*(H2.bc-hatH_causal)/hatH_causal
    b2H_3<-100*(H2.bc-hrH_true)/hrH_true
    
    b0H_1<-100*(H_obs-H_true)/H_true
    b0H_2<-100*(H_obs-hatH_causal)/hatH_causal
    b0H_3<-100*(H_obs-hrH_true)/hrH_true
  }
  )
  
  df_bc<-NULL
  # hr(H_true) target
  res_new<-data.table(res_new)
  df_res<-res_new[,c("b1H_1")]
  df_res$est<-"BC_1"
  df_res$target<-target[1]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b2H_1")]
  df_res$est<-"BC_2"
  df_res$target<-target[1]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b0H_1")]
  df_res$est<-"Un-adj"
  df_res$target<-target[1]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  # hatH_causal target
  df_res<-res_new[,c("b1H_2")]
  df_res$est<-"BC_1"
  df_res$target<-target[2]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b2H_2")]
  df_res$est<-"BC_2"
  df_res$target<-target[2]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b0H_2")]
  df_res$est<-"Un-adj"
  df_res$target<-target[2]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  # H_causal fixed
  df_res<-res_new[,c("b1H_3")]
  df_res$est<-"BC_1"
  df_res$target<-target[3]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b2H_3")]
  df_res$est<-"BC_2"
  df_res$target<-target[3]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  
  df_res<-res_new[,c("b0H_3")]
  df_res$est<-"Un-adj"
  df_res$target<-target[3]
  names(df_res)[1]<-"bias"
  df_bc<-rbind(df_bc,df_res)
  return(df_bc)
}

df_Hplots<-function(df_res,dgm,return_plot=TRUE,qmax=c(1,1,1,1),alpha_sg=0.025){
dfA<-as.data.frame(df_res)
dfA<-na.omit(dfA)
# Observed
temp<-getci_Cox(df=na.omit(dfA),est="H_obs",se="seH_obs",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
Hobs_l<-unlist(temp$lb)
Hobs_u<-unlist(temp$ub)
Hobs_cover<-unlist(temp$cover)
# Bias-corrected
temp<-getci_Cox(df=na.omit(dfA),est="H1.bc",se="seH1.bc",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
Hbc_l<-unlist(temp$lb)
Hbc_u<-unlist(temp$ub)
Hbc_cover<-unlist(temp$cover)
# True
temp<-getci_Cox(df=na.omit(dfA),est="H_true",se="seH_true",target=unlist(dgm["hr.H.true"]),alpha=alpha_sg)
H_l<-unlist(temp$lb)
H_u<-unlist(temp$ub)
H_cover<-unlist(temp$cover)
nStar<-nrow(dfA)
dH_plot <- data.frame(
  x=c(1:(nStar*3)), 
  value1=c(H_l,Hobs_l,Hbc_l), 
  value2=c(H_u,Hobs_u,Hbc_u),
  value3=c(H_cover,Hobs_cover,Hbc_cover),
  ests=c(rep('H_true',nStar),rep('H_un',nStar), rep('H_bc', nStar))
)
nullH.value<-dgm$hr.H.true
# Hc
# Observed
temp<-getci_Cox(df=na.omit(dfA),est="Hc_obs",se="seHc_obs",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
Hobs_l<-unlist(temp$lb)
Hobs_u<-unlist(temp$ub)
Hobs_cover<-unlist(temp$cover)
# Bias-corrected
temp<-getci_Cox(df=na.omit(dfA),est="Hc1.bc",se="seHc1.bc",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
Hbc_l<-unlist(temp$lb)
Hbc_u<-unlist(temp$ub)
Hbc_cover<-unlist(temp$cover)
# True
temp<-getci_Cox(df=na.omit(dfA),est="Hc_true",se="seHc_true",target=unlist(dgm["hr.Hc.true"]),alpha=alpha_sg)
H_l<-unlist(temp$lb)
H_u<-unlist(temp$ub)
H_cover<-unlist(temp$cover)

dHc_plot <- data.frame(
  x=c(1:(nStar*3)), 
  value1=c(H_l,Hobs_l,Hbc_l), 
  value2=c(H_u,Hobs_u,Hbc_u),
  value3=c(H_cover,Hobs_cover,Hbc_cover),
  ests=c(rep('Hc_true',nStar),rep('Hc_un',nStar), rep('Hc_bc', nStar))
  )
nullHc.value<-dgm$hr.Hc.true

plot_H_Hc<-NULL
if(return_plot){
  
  dfP<-subset(dH_plot,value3==0)
  
  # Modify ylim to 95% quantile
   maxy<-quantile(dfP$value2,qmax[1])
   miny<-min(dfP$value1)
  
    plot1<-ggplot(dfP) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=nullH.value, linetype='dashed', color='black', size=0.5)+
    coord_flip(ylim=c(miny,maxy))+
    theme()+
    labs(title = "", colour="Est") +
    theme(legend.position = "top",
          legend.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
    
  dfP<-subset(dH_plot,value3==1)
  
  maxy<-quantile(dfP$value2,qmax[2])
  miny<-min(dfP$value1)
  
  plot2<-ggplot(dfP) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=nullH.value, linetype='dashed', color='black', size=0.5) +
    coord_flip(ylim=c(miny,maxy))+
    theme()+
    labs(title = "", colour="Est") +
    theme(legend.position = "top",
          legend.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  dfP<-subset(dHc_plot,value3==0)
  
  maxy<-quantile(dfP$value2,qmax[3])
  miny<-min(dfP$value1)
  
  plot3<-ggplot(dfP) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=nullHc.value, linetype='dashed', color='black', size=0.5) +
    coord_flip(ylim=c(miny,maxy))+
    theme()+
    labs(title = "", colour="Est") +
    theme(legend.position = "top",
          legend.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  
  dfP<-subset(dHc_plot,value3==1)
  
  maxy<-quantile(dfP$value2,qmax[4])
  miny<-min(dfP$value1)
  
  plot4<-ggplot(dfP) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey")+ 
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=nullHc.value, linetype='dashed', color='black', size=0.5) +
    coord_flip(ylim=c(miny,maxy))+
    theme()+
    labs(title = "", colour="Est") +
    theme(legend.position = "top",
          legend.text = element_text(size = 8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
plot_H_Hc<-plot_grid(plot1,plot2,plot3,plot4,labels=c("(a)","(b)","(c)","(d)"))
}

return(list(dH_plot=dH_plot,dHc_plot=dHc_plot,nullH.value=nullH.value,nullHc.value=nullHc.value,plot_H_Hc=plot_H_Hc,nStar=nStar))
}


#order_analysis<-c("H true","FS","GRF.60","VT(24)","VT(36)")
#plot_breaks<-c("H true","FS","GRF.60","VT(24)","VT(36)")
#col_values<-c("brown","blue","orange","lightblue","darkgrey")
#details=TRUE
#dig_Hc=2
#dig_H=2
#hrH_this=2.5

# Power plots
df_PowerPlots<-function(flist,order_analysis=c("H true","FS","GRF.60","VT(24)","VT(36)"),
plot_breaks=c("H true","FS","GRF.60","VT(24)","VT(36)"),col_values=c("brown","blue","orange","lightblue","darkgrey"),
order_analysis_b=c("H true","FS","GRF.60","ITT"),
plot_breaks_b=c("H true","FS","GRF.60","ITT"),
col_values_b=c("brown","blue","orange","beige"),
details=TRUE,dig_Hc=2,dig_H=2,hrH_this=2.5){
res_all<-NULL
for(ff in 1:length(flist)){
  load(paste0(resdir,flist[ff]))
  # Add "oracle" analysis
if(details)  cat("dgm: H,H-complement,ITT=",c(dgm$hr.H.true,dgm$hr.Hc.true,dgm$hr.causal),"\n")  
  res_oracle<-subset(res,analysis=="FS")  
if(details)  cat("# of FS simulations=",c(nrow(res_oracle)),"\n")  
  res_GRF<-subset(res,analysis=="GRF.60")  
if(details)  cat("# of GRF.60 simulations=",c(nrow(res_GRF)),"\n")  
  # Set estimates (oracle) to truths
  res_oracle$hr.H.hat<-res_oracle$hr.H.true
  res_oracle$l.H.hat<-res_oracle$l.H.true
  res_oracle$u.H.hat<-res_oracle$u.H.true
  res_oracle$hr.Hc.hat<-res_oracle$hr.Hc.true  
  res_oracle$l.Hc.hat<-res_oracle$l.Hc.true  
  res_oracle$u.Hc.hat<-res_oracle$u.Hc.true  
  res_oracle$analysis<-c("H true")
  # Add ITT
  res_itt<-subset(res,analysis=="FS")
  res_itt$hr.Hc.hat<-res_itt$hr.itt
  res_itt$analysis<-c("ITT")
  res_new<-rbind(res,res_oracle,res_itt)
  res_new$hrH<-round(dgm$hr.H.true,3)
  res_new$hrHc<-round(dgm$hr.Hc.true,3)
  res_all<-rbind(res_all,res_new)
}

if(details){
print(table(res_all$analysis))
print(names(res_all))
}

res_all<-within(res_all,{
  rej1<-ifelse(!is.na(l.H.hat) & l.H.hat>1.0,1,0)
  rej2<-ifelse(!is.na(u.Hc.hat) & u.Hc.hat<1.0,1,0)
  rej12<-rej1*rej2
  anyH<-ifelse(!is.na(hr.H.hat),1,0)
  hrH<-ifelse(is.na(hrH),0.0,hrH)
  ppvHc<-npv
  sensHc<-specificity
  ppvH<-ppv
  sensH<-sensitivity
  ppv_Hc<-factor(cut2.factor(npv,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
  sens_Hc<-factor(cut2.factor(specificity,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
  ppv_H<-factor(cut2.factor(ppv,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
  sens_H<-factor(cut2.factor(sensitivity,cutpoints=c(0,0.9,0.95,1),labels=c("<=90%","90-95%",">95%")))
})

# Restrict to FS, GRF.60, and VT(24)
dfPlots<-subset(res_all,analysis %in% order_analysis)

dfPlots<-within(dfPlots,
                {Hc_rel<-100*(hr.Hc.hat-hrHc)/hrHc
                H_rel<-100*(hr.H.hat-hrH)/hrH
                Est<-factor(analysis,levels=order_analysis)})

# Remove H true analysis as classification metrics not applicable
dfP_ppvHc<-subset(dfPlots,analysis!="H true") %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(ppvHc=as.factor(ppv_Hc)) %>%  
  group_by(Est,ppvHc) %>%
  summarize(Hc_rbias=mean(Hc_rel))

dfP_sensHc<-subset(dfPlots,analysis!="H true") %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(sensHc=as.factor(sens_Hc)) %>%  
  group_by(Est,sensHc) %>%
  summarize(Hc_rbias=mean(Hc_rel))


plot_sensHc<-subset(dfPlots, analysis!="H true") %>% ggplot(aes(sens_Hc,Hc_rel,color=Est))+
  geom_boxplot()+
  scale_color_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
  # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
  #coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="sens(Hc)", y="Hc relative bias")+
  theme(legend.position = "top",
        legend.text = element_text(size = 4))

plot_ppvHc<-subset(dfPlots, analysis!="H true") %>% ggplot(aes(ppv_Hc,Hc_rel,color=Est))+
  geom_boxplot()+
  scale_color_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
  # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
  #coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="ppv(Hc)", y="Hc relative bias")+
  theme(legend.position = "top",
        legend.text = element_text(size = 4))

# For H, restrict to "estimable H's"
# Otherwise, NA's
plot_sensH<-subset(dfPlots, analysis!="H true" & !is.na(hr.H.hat)) %>% ggplot(aes(sens_H,H_rel,color=Est))+
  geom_boxplot()+
  scale_color_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
  # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
  #coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="sens(H)", y="H relative bias")+
  theme(legend.position = "top",
        legend.text = element_text(size = 4))

plot_ppvH<-subset(dfPlots, analysis!="H true" & !is.na(hr.H.hat)) %>% ggplot(aes(ppv_H,H_rel,color=Est))+
  geom_boxplot()+
  scale_color_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=c(0.0,10), linetype='dashed', color=c('black','grey'), linewidth=c(0.5,0.25))+
  # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
  #coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="ppv(H)", y="H relative bias")+
  theme(legend.position = "top",
        legend.text = element_text(size = 4))


dfP_0<-dfPlots %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(hrH=as.factor(hrH)) %>%  
  group_by(Est,hrH) %>%
  summarize(anyH_rate=mean(anyH))

dfP_1<-dfPlots %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(hrH=as.factor(hrH)) %>%  
  group_by(Est,hrH) %>%
  summarize(rej_rate=mean(rej1))

# Conditional on finding H
dfPlots_H<-subset(dfPlots,!is.na(hr.H.hat))

dfP_2<-dfPlots_H %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(hrH=as.factor(hrH)) %>%  
  group_by(Est,hrH) %>%
  summarize(rej_rate=mean(rej2))

dfP_12<-dfPlots %>% 
  mutate(Est=factor(analysis,levels=order_analysis)) %>%
  mutate(hrH=as.factor(hrH)) %>%  
  group_by(Est,hrH) %>%
  summarize(rej_rate=mean(rej12))

plot0<-dfP_0 %>% ggplot(aes(hrH,anyH_rate,fill=Est))+
  geom_col(position="dodge")+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
 # geom_text(aes(label=round(anyH_rate,3)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1)+
  coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="Hazard ratio for H", y="Chance of identifying any H")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6))

plot1<-dfP_1 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
  geom_col(position="dodge")+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
 #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
  coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="Hazard ratio for H", y="Power for H > 1")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6))

# Remove x-axis here
plot2<-dfP_2 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
  geom_col(position="dodge")+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
 #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
  coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="Hazard ratio for H", y="Conditional power for H-complement < 1")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title.y = element_blank())

plot12<-dfP_12 %>% ggplot(aes(hrH,rej_rate,fill=Est))+
  geom_col(position="dodge")+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  geom_hline(yintercept=0.025, linetype='dashed', color='black', linewidth=0.5)+
 #geom_text(aes(label=round(rej_rate,2)),color="white",position=position_dodge(0.9),vjust=0.5,hjust=1,size=2.5)+
  coord_flip()+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="Hazard ratio for H", y="Power for {H > 1 & Hc < 1}")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6),
        axis.title.y = element_blank())


plot_all<-plot_grid(plot1,plot2,plot12,labels=c("(a)","(b)","(c)"), ncol=2, byrow=TRUE, align="hv")


res_all_b<-within(res_all,{
analysis2<-factor(analysis,levels=order_analysis_b)})

if(!any(res_all$analysis=="FS-M")){
df_FS<-subset(res_all_b, analysis=="FS" & !is.na(l.Hc.hat) & hrH==hrH_this)
df_GRF<-subset(res_all_b, analysis=="GRF.60" & !is.na(l.Hc.hat) & hrH==hrH_this)
df_true<-subset(res_all_b, analysis=="H true" & !is.na(l.Hc.hat) & hrH==hrH_this)
df_itt<-subset(res_all_b, analysis=="ITT" & hrH==hrH_this)
df_3<-rbind(df_FS,df_GRF,df_true)
# Include ITT only for complement
df_3b<-rbind(df_FS,df_GRF,df_true,df_itt)
}

if(any(res_all$analysis=="FS-M")){
  df_FS<-subset(res_all_b, analysis=="FS" & !is.na(l.Hc.hat) & hrH==hrH_this)
  df_FS_M<-subset(res_all_b, analysis=="FS-M" & !is.na(l.Hc.hat) & hrH==hrH_this)
  df_GRF<-subset(res_all_b, analysis=="GRF.60" & !is.na(l.Hc.hat) & hrH==hrH_this)
  df_true<-subset(res_all_b, analysis=="H true" & !is.na(l.Hc.hat) & hrH==hrH_this)
  df_itt<-subset(res_all_b, analysis=="ITT" & hrH==hrH_this)
  df_3<-rbind(df_FS,df_FS_M,df_GRF,df_true)
  # Include ITT only for complement
  df_3b<-rbind(df_FS,df_FS_M,df_GRF,df_true,df_itt)
}

#table(df_3$analysis)

# Change analysis name to match above
df_3$Est<-df_3$analysis2
df_3b$Est<-df_3b$analysis2

if(details) print(table(df_3b$analysis))

# Relative bias
plot3<-ggplot(df_3b,aes(Hc_rel,fill=Est))+
  geom_density(alpha=0.7)+
  geom_vline(xintercept=c(0.0,7.5), linetype='dashed', color='black', linewidth=c(0.15,0.35))+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="% Relative bias for H-complement", y="Density")+
  theme(legend.position = "top",
  legend.text = element_text(size = 6), axis.title.y = element_blank())


# Hc Bias
plot4<-ggplot(df_3b,aes(hr.Hc.hat,fill=Est))+
  geom_density(alpha=0.7)+
  geom_vline(xintercept=c(round(dgm$hr.Hc.true,dig_Hc)), linetype='dashed', color='black', linewidth=c(0.35))+
  scale_fill_manual(breaks=plot_breaks_b,values=col_values_b)+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="HR estimates (Hc and ITT)", y="Density")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6), axis.title.y = element_blank())


plot4b<-ggplot(df_3,aes(hr.H.hat,fill=Est))+
  geom_density(alpha=0.7)+
  geom_vline(xintercept=hrH_this, linetype='dashed', color='black', linewidth=c(0.35))+
  scale_fill_manual(breaks=plot_breaks,values=col_values)+
  theme()+
  labs(title = "", colour="Est") +
  labs(x="HR estimates (H)", y="Density")+
  theme(legend.position = "top",
        legend.text = element_text(size = 6), axis.title.y = element_blank())

plot_all2<-plot_grid(plot1,plot2,plot12,plot4,labels=c("(a)","(b)","(c)","(d)"), ncol=2, byrow=TRUE, align="hv")

#if(details) plot(plot_all2)

plot_all3<-plot_grid(plot4b,plot4,plot1,plot12,labels=c("(a)","(b)","(c)","(d)"), ncol=2, byrow=TRUE, align="hv")

if(details){
plot(plot_all3)
#cat("Hc bias: OLS vs size Hc + ppv +npv + specificity +  sensitivity","\n")
#ols_Hcbias<-lm(Hc_rel~ size.Hc+ppv+npv+specificity+sensitivity,data=dfPlots)
#print(summary(ols_Hcbias))
#cat("Hc bias: OLS vs size H + ppv +npv + specificity +  sensitivity","\n")
#ols_Hcbias<-lm(Hc_rel~ size.H+ppv+npv+specificity+sensitivity,data=dfPlots)
#print(summary(ols_Hcbias))
#ols_Hcbias<-lm(Hc_rel~ specificity,data=dfPlots)
#print(summary(ols_Hcbias))

#ols_Hbias<-lm(H_rel~ size.Hc+ppv+npv+specificity+sensitivity,data=dfPlots)
#cat("H bias: OLS vs size Hc + ppv +npv + specificity +  sensitivity","\n")
#print(summary(ols_Hbias))
#cat("H bias: OLS vs size H + ppv +npv + specificity +  sensitivity","\n")
#ols_Hbias<-lm(H_rel~ size.H+ppv+npv+specificity+sensitivity,data=dfPlots)
#print(summary(ols_Hbias))
#ols_Hbias<-lm(H_rel~sensitivity,data=dfPlots)
#print(summary(ols_Hbias))
}


return(list(dgm=dgm,res_all=res_all,dfPlots=dfPlots,plot0=plot0,plot_all2=plot_all2,
            plot4b=plot4b,plot4=plot4,plot_sensHc=plot_sensHc,plot_ppvHc=plot_ppvHc,
            plot_sensH=plot_sensH,plot_ppvH=plot_ppvH,
            df_3=df_3,df_3b=df_3b,plot_all3=plot_all3,plot1=plot1,plot2=plot2,plot12=plot12,
            order_analysis=order_analysis,order_analysis_b=order_analysis_b))
}

# Return avg, SD rounded to 3
SummaryStat<-function(df,name,sigdig=2,includeSD=FALSE,showSD=TRUE){
if(is.data.table(df)){
df<-as.data.frame(df)
}
# include "SD=" in output
if(!includeSD){  
temp<-na.omit(df[,c(name)])
m<-round(mean(temp),sigdig)
sig<-round(sqrt(var(temp)),sigdig)
# Check if binary
if(all(as.numeric(temp) %in% c(0,1))){
sig<-round(sqrt(m*(1-m)/length(temp)),sigdig)
}
if(sigdig<=4 & m<=0.001) m<-"<0.001"
if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
out<-paste0(m," (")
out<-paste0(out,sig)
out<-paste0(out,")")
}
  if(includeSD){  
    temp<-na.omit(df[,c(name)])
    m<-round(mean(temp),sigdig)
    sig<-round(sqrt(var(temp)),sigdig)
    # Check if binary
    if(all(as.numeric(temp) %in% c(0,1))){
      sig<-round(sqrt(m*(1-m)/length(m)),sigdig)
    }
    if(sigdig<=4 & m<=0.001) m<-"<0.001"
    if(sigdig<=4 & sig<=0.001) sig<-"<0.001"
    out<-paste0(m," (SD=")
    out<-paste0(out,sig)
    out<-paste0(out,")")
  }
if(!showSD) out<-c(m)  
return(out)
}

DiffRate<-function(df1,df2,name,sigdig=2){
  if(is.data.table(df1)){
    df1<-as.data.frame(df1)
    df2<-as.data.frame(df2)
      }
temp<-na.omit(df1[,c(name)])
m1<-mean(temp)
temp<-na.omit(df2[,c(name)])
m2<-mean(temp)
out<-round(100*(m2-m1),1)
return(out)
}
  
pow_size<-function(df,minsize=0,sigdig=3){
rej<-with(df,mean(rej12 & size.Hc>=minsize))
return(round(rej,sigdig))
}



