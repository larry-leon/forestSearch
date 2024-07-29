
# Some plotting functions

# For plotting subgroups and the complement

plot.subgroup<-function(tte.name,event.name,treat.name,wgt.name=NULL,sub1,sub1C,xloc1=NULL,xloc2=NULL,details=FALSE,show.logrank=FALSE,
                        ymin=0,cox.cex=0.7,prob.points=c(0,0.25,0.5,0.75,1.0),
                        cex_Yaxis=0.8,title1=NULL,title2=NULL,choose_ylim=TRUE,
                        exp.lab="Treat",con.lab="Control",legend.cex=0.70,risk.cex=0.70,yloc1=0.6,yloc2=0.6,subid=NULL,byrisk=2,fix.rows=TRUE,show.med=FALSE,
                        xlab="Months", ylab="Survival"){
  
  
  if(!is.null(subid) & (!is.null(title1) | !is.null(title2))) stop("subid OR titles need to be specified (not both)")

  # Set tpoints.add to span range of both subgroups
  aa<-max(sub1[,c(tte.name)])
  bb<-max(sub1C[,c(tte.name)])
  tpoints.add<-c(-1,max(aa,bb))
  
 if(fix.rows){
 old_par = par(no.readonly = TRUE)
 layout(matrix(1:2, ncol=2))
}

 par(mar = old_par$mar + c(0, 0, 0, -2))
 
  
  df<-sub1
  
  km.fit <- KM.plot.2sample.weighted(df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
                                    tpoints.add=tpoints.add,
                                    risk.set=TRUE, by.risk=byrisk, risk.cex=risk.cex, censor.cex=risk.cex,
                                    risk_offset=0.15, risk_delta=0.075,
                                    Xlab=xlab,Ylab=ylab, cex_Yaxis=cex_Yaxis,
                                    show.ticks=TRUE,
                                    col.1="black", col.2="blue",
                                    arms = c(exp.lab,con.lab), arm.cex=risk.cex,
                                    ltys=c(1,1),lwds=c(2,2),
                                    show.logrank=FALSE, lr.digits=2, put.legend.lr="top",
                                    qlabel="m =",
                                    show.med=TRUE, med.cex=risk.cex-0.15,
                                    show.cox=TRUE, cox.digits=3, cox.cex=risk.cex-0.10)
  
 
 

  if(!is.null(subid)) title(main=subid,cex.main=0.80)
  if(!is.null(title1)) title(main=title1,cex.main=0.80)

  par(mar = old_par$mar + c(0, -2, 0, 0))
  
  df<-sub1C
  
  km.fit <- KM.plot.2sample.weighted(df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
                                     tpoints.add=tpoints.add,
                                     risk.set=TRUE, by.risk=byrisk, risk.cex=risk.cex, censor.cex=risk.cex,
                                     risk_offset=0.15, risk_delta=0.075,
                                     Xlab=xlab,Ylab="", show.Y.axis=FALSE, cex_Yaxis=cex_Yaxis,
                                     show.ticks=TRUE,
                                     col.1="black", col.2="blue",
                                    show_arm_legend=FALSE, arms = c(exp.lab,con.lab), arm.cex=risk.cex,
                                     ltys=c(1,1),lwds=c(2,2),
                                     show.logrank=FALSE, lr.digits=2, put.legend.lr="top",
                                     qlabel="m =",
                                     show.med=TRUE, med.cex=risk.cex-0.15,
                                     show.cox=TRUE, cox.digits=3, cox.cex=risk.cex-0.10)

  if(!is.null(subid)) title(main="Complement",cex.main=0.80)
  if(!is.null(title2)) title(main=title2,cex.main=0.80)
  
  par(old_par)
  # Restore layout.
  layout(1:1)
  
}


#plot.band<-function(x,mean.value,lower,upper,show.axes=F,band=TRUE,ltype="l",lty=1,xlabel=NULL,ylabel=NULL,color="grey",ylim=c(min(lower,na.rm=TRUE),max(upper,na.rm=TRUE))){
#  plot(x[order(x)],mean.value[order(x)],type="n",axes=show.axes,xlab=xlabel,lty=lty,
#       ylab=ylabel,ylim=ylim)
#  if(band){
#    polygon(c(x[order(x)],rev(x[order(x)])),c(lower[order(x)],rev(upper[order(x)])),col=color,border=FALSE)
#    lines(x[order(x)],mean.value[order(x)],lty=lty,lwd=2.5,type=ltype)
#  }
#}


# KM plotting functions (core)


N.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error<=x))
}

R.Weighted<-function(x,error,W=rep(1,length(error))){
  sum(W*(error>=x))
}

NA.CHR.Weighted<-function(time,Delta,W.n=rep(1,length(time)),W.d=rep(1,length(time)),
                          at.points=sort(time),se.type="greenwood",get.Stute=FALSE,tpoints.add=NULL){
  
  if(!is.null(tpoints.add)) at.points<-sort(c(unique(c(at.points,tpoints.add))))
  
  if(se.type!="greenwood" & se.type!="tsiatis") stop("Invalid se type -- greenwood or tsiatis allowed")
  #is.sorted<-(all(time==sort(time)))
  is.sorted<-!is.unsorted(time)
  if(!is.sorted){
    id<-order(time); time<-time[id]; Delta<-Delta[id]; W.n<-W.n[id]; W.d<-W.d[id]
  }
  risk<-unlist(lapply(as.list(at.points),R.Weighted,error=time,W=W.d))
  ###########################################################################
  ### Adaptive H process correspoinding to N-A rep. via integral wrt M.G ####
  Hmart.chf<-ifelse(risk>0,1/risk,0)
  ############################################################################
  counting<-unlist(lapply(as.list(at.points),N.Weighted,error=time,W=W.n*ifelse(Delta==1,1,0)))
  counting <- c(0, counting)
  dN<-diff(counting)
  dN.risk<-ifelse(risk>0,dN/risk,0.0)
  chf <- cumsum(dN.risk)
  var.chf<-cumsum(ifelse(risk>0,dN/(risk^2),0.0))
  S.KM <- cumprod(1-dN.risk)
  
  S.KM[which(S.KM<0)]<-0.0
  
  S.NA <- exp(-chf)
  var.NA<-(S.NA^2)*var.chf
  # Greenwood variance estimate
  if(se.type=="greenwood"){
    aa<-dN
    bb<-risk*(risk-dN)
    var.KM<-(S.KM^2)*cumsum(ifelse(risk>0,aa/bb,0.0))
    se.KM<-sqrt(var.KM)
  }
  if(se.type=="tsiatis"){
    var.KM<-(S.KM^2)*var.chf
    se.KM<-sqrt(var.KM)
  }
  result<-list(time=time,at.points=at.points,S.NA=S.NA,S.KM=S.KM,chf=chf,se.chf=sqrt(var.chf),
               se.NA=sqrt(var.NA),dN.risk=dN.risk,
               n.risk=risk,dN=dN,
               Hmart.chf=Hmart.chf,se.KM=se.KM)
  return(result)
}

# KM plot functions
###################################################################################################################################################
# K-M is weighted via Weight such as propensity scores (default is Weight==1 for standard KM)
# If censor.mark.all=TRUE, then all censorings are marked which replicates survminer
# If censor.mark.all=FALSE, then curves are marked at each censoring time which is not also a death time

# Default ymin=0, ymax=1 probability range 
# "risk_offset" sets space for risk-set by offsetting ymin to ymin-risk_offset
# Eg, if ymin=0 and risk_offset=0.075 then a space of 0.075 is allowed for risk-set
# The ymin-risk_offset sets risk-set y-coordinate for the treatment arm ("experimental")
# risk_delta is distance between treatment and control arm risk-sets
# Eg, if risk_delta=0.025 then space of risk_offset+risk_delta is provided
# with distance between y-corrdinates of 0.025
# "risk_offset" sets space for risk-set by offsetting ymin to ymin-risk_offset
# Eg, if ymin=0 and risk_offset=0.075 then a space of 0.075 is allowed for risk-set
# The ymin-risk_offset sets risk-set y-coordinate for the treatment arm ("experimental")
# risk_delta is distance between treatment and control arm risk-sets
# Eg, if risk_delta=0.025 then space of risk_offset+risk_delta = 0.1 is provided
# with distance between y-corrdinates of 0.025
###################################################################################################################################################


KM.plot.2sample.weighted<-function(df,tte.name,event.name,treat.name,weight.name=NULL,
                                   strata.name=NULL,show.cox=FALSE,cox.cex=0.725,show.logrank=FALSE,
                                   logrank.cex=0.725, cox.eps=0.001,
                                   lr.eps=0.001,show_arm_legend=TRUE,arms=c("treat1","treat2"),
                                   stop.onerror=TRUE,check.KM=TRUE,details=FALSE,
                                   put.legend.cox='topright',put.legend.lr='top',lr.digits=2,cox.digits=3, 
                                   tpoints.add=c(0),by.risk=NULL,Xlab="time",Ylab="proportion surviving",col.1="black",col.2="blue",show.med=FALSE,
                                   choose_ylim =FALSE, arm.cex=0.725,
                                   quant=0.5, qlabel="median =", med.cex=0.725, ymed.offset=0.25, del.med=0.075, xmed.offset=4,
                                   risk.cex=1, plotit=TRUE,ltys=c(1,1),lwds=c(2,2),
                                   censor.mark.all=TRUE,censor.cex=0.725, show.ticks=TRUE,risk.set=TRUE,
                                   ymin=0,ymax=1, ymin.del=0.035,
                                   ymin2=NULL, risk_offset=0.075, risk_delta=0.025, y.risk2=NULL,show.Y.axis=TRUE,cex_Yaxis=1,
                                   y.risk1=NULL,add.segment=FALSE,risk.add=NULL,xmin=0,xmax=NULL,x.truncate=NULL,
                                   time.zero=0.0, time.zero.label=0.0, prob.points=NULL){
  
  Y <- df[,c(tte.name)]
  E <- df[,c(event.name)]
  Treat <- df[,c(treat.name)]
  if(!is.null(strata.name)) Strata <- df[,c(strata.name)]
  if(is.null(strata.name)) Strata <- rep(1,length(Y))
  
  if(!is.null(weight.name)) Weight <- df[,c(weight.name)]
  
  if(is.null(weight.name)) Weight <- rep(1,length(Y))
  
  tpoints.add.plot <- tpoints.add  
  # Add all events to estimation time points
  # append to tpoints.add with "-1" added to allow for risk-set space
  evs<-sort(unique(Y[which(E==1)]))
  # append to tpoints.add
  # and add max observed event
  tpoints.add<-sort(unique(c(time.zero,tpoints.add,evs,max(Y))))
  
  Weight_orig <- Weight
  x1<-Y[which(Treat==1)]
  eta1<-E[which(Treat==1)]
  w1<-Weight[which(Treat==1)]
  
  x2<-Y[which(Treat==0)]
  eta2<-E[which(Treat==0)]
  w2<-Weight[which(Treat==0)]
  
  coxfit<-coxph(Surv(Y,E)~Treat+strata(Strata))  
  hr.ci<-as.matrix(round(summary(coxfit)$conf.int,3))
  pval.cox<-summary(coxfit)$coefficients[5]
  
  # ipw weighting
  if(!all(Weight_orig==1)){ 
    coxfit <- coxph(Surv(Y,E)~Treat+strata(Strata), weights=Weight_orig, robust=TRUE)
    hr.ci<-as.matrix(round(summary(coxfit)$conf.int,3))
    pval.cox<-summary(coxfit)$coefficients[6]
  }
  
  # Get HR for Treat
  hr.ci<-hr.ci[1,c(1,3,4)]
  
  if(details){
    if(all(Weight_orig==1)){
      cat("Cox un-adjusted HR=",c(paste(hr.ci[1])),"\n")
      cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")
    }
    
    if(!all(Weight_orig==1)){
      cat("Cox ipw-adjusted HR=",c(paste(hr.ci[1])),"\n")
      cat("Cox CIs=",c(paste(hr.ci[2]),paste(hr.ci[3])),"\n")
    }
  }
  
  if(is.null(by.risk)){
    tt<-sort(unique(Y))
    by.risk<-quantile(tt-c(0,tt[-length(tt)]),c(0.75),na.rm=TRUE)
  }
  
  # For log-rank and survival estimates at common points    
  at.points<-sort(c(unique(c(x1,x2,tpoints.add))))
  at.points1<-sort(c(unique(c(x1,tpoints.add))))
  at.points2<-sort(c(unique(c(x2,tpoints.add))))
  
  # Experimental
  fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1)
  S1.KM<-fit1$S.KM
  SE1.KM<-fit1$se.KM
  
  # Control
  fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2)
  S2.KM<-fit2$S.KM
  SE2.KM<-fit2$se.KM
  
  if(choose_ylim){
    ymin <- min(c(S1.KM,S2.KM))
  }
  
  if(is.null(ymin2)) ymin2 <- ymin-risk_offset  
  if(is.null(y.risk1)) y.risk1 <- ymin2
  if(is.null(y.risk2)) y.risk2 <- ymin2+risk_delta
  
  # Points for plotting
  #tpoints.add.plot <- tpoints.add
  # Same at at.points is default
  at.points1.plot<-sort(c(unique(c(x1,tpoints.add.plot))))
  at.points2.plot<-sort(c(unique(c(x2,tpoints.add.plot))))
  
  # LR() denotes the log-rank scale
  # Integral of "LR differences" = log-rank test (non-standardized)
  ZLR.dpoints<-LR.dpoints<-S1.dpoints<-SE1.dpoints<-S2.dpoints<-SE2.dpoints<-NULL
  
  tp1.match<-match(tpoints.add,at.points1)
  tp2.match<-match(tpoints.add,at.points2)
  
  S1.dpoints<-S1.KM[tp1.match]
  SE1.dpoints<-SE1.KM[tp1.match]
  Y1.dpoints<-fit1$n.risk[tp1.match]
  dN1.dpoints<-fit1$dN[tp1.match]
  
  S2.dpoints<-S2.KM[tp2.match]
  SE2.dpoints<-SE2.KM[tp2.match]
  Y2.dpoints<-fit2$n.risk[tp2.match]
  dN2.dpoints<-fit2$dN[tp2.match]
  
  K<-(Y1.dpoints*Y2.dpoints)/(Y1.dpoints+Y2.dpoints)
  
  Y.dpoints<-Y1.dpoints+Y2.dpoints
  dN.dpoints<-dN1.dpoints+dN2.dpoints
  
  term1<-cumsum(ifelse(Y1.dpoints>0,(K/Y1.dpoints)*dN1.dpoints,0.0))
  term2<-cumsum(ifelse(Y2.dpoints>0,(K/Y2.dpoints)*dN2.dpoints,0.0))
  
  # variance
  h0<-ifelse(Y1.dpoints==0,0,(K^2/Y1.dpoints))
  h1<-ifelse(Y2.dpoints==0,0,(K^2/Y2.dpoints))
  dJ<-ifelse(Y.dpoints==1,0,(dN.dpoints-1)/(Y.dpoints-1))
  dL<-ifelse(Y.dpoints==0,0,dN.dpoints/Y.dpoints)
  sig2s<-(h0+h1)*(1-dJ)*dL
  sig2<-cumsum(sig2s)
  
  LR.dpoints<-term2-term1
  ZLR.dpoints<-LR.dpoints/sqrt(sig2)
  # Zlr() is weighted log-rank process
  # Zlr(tau) for tau end of followup is 
  # the "final" z statistic (1-sided)
  z.lr<-ZLR.dpoints[length(ZLR.dpoints)]
  p.lr<-1-pnorm(z.lr)
  
  if(check.KM & all(Weight==1)){
    p.2side <- 1-pchisq(z.lr^2,1)
    logrnk<-survdiff(Surv(Y,E)~Treat)
    p.val<-1-pchisq(logrnk$chisq,length(logrnk$n) - 1)
    if(round(p.2side-p.val,8)>0){
      warning("Discrepancy with log-rank process and survdiff")
    }
    if(details){
      cat("My p-value and survdiff=",c(p.2side,p.val),"\n")  
      cat("My z^2 and survdiff=",c(z.lr^2,logrnk$chisq),"\n")
    }    
  }  
  
  #########################################
  # Note: want to confirm with R function
  #########################################
  #           Check KM fits               #
  if(check.KM){
    # only at observed unique observations to align with survfit
    # Experimental
    S1.check<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=sort(unique(x1)))$S.KM
    S2.check<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=sort(unique(x2)))$S.KM
    KM.fit1<-survfit(Surv(x1,eta1)~1,weights=w1)
    max.error<-round(max(abs(KM.fit1$surv-S1.check)),8)
    if(details) cat("Treat=1: Max error max|KM(survfit)[t]-KM(mine)[t]| =",c(max.error),"\n")
    if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM1 fit")
    KM.fit2<-survfit(Surv(x2,eta2)~1,weights=w2)
    max.error<-round(max(abs(KM.fit2$surv-S2.check)),8)
    if(details) cat("Treat=0: Max error max|KM(survfit)[t]-KM(mine)[t]| =",c(max.error),"\n")
    if(stop.onerror & max.error>=0.00001) stop("Discrepancy in KM2 fit")
    #        End Check KM fits               #
  }
  
  events1<-sort(unique(x1[eta1==1]))
  events2<-sort(unique(x2[eta2==1]))
  
  if(!censor.mark.all){
    # Censored points
    # Curves are marked at each censoring time which is not also a death time
    x1.cens<-x1[eta1==0] 
    x2.cens<-x2[eta2==0]
    x1.cens<-setdiff(x1.cens,events1)  # Censorings that are NOT an event (this is the R version which follows estimation convention)
    x2.cens<-setdiff(x2.cens,events2)
    x1.match<-match(x1.cens,at.points1)
    x2.match<-match(x2.cens,at.points2)
  }
  
  if(censor.mark.all){
    # Censored points
    x1.cens<-x1[eta1==0] 
    x2.cens<-x2[eta2==0]
    x1.match<-match(x1.cens,at.points1)
    x2.match<-match(x2.cens,at.points2)
  }
  
  risk.points<-round(c(seq(time.zero.label,max(at.points),by=by.risk)))
  risk.points<-sort(unique(c(risk.points,risk.add)))
  
  risk.2<-unlist(lapply(as.list(risk.points),R.Weighted,error=x2,W=w2))
  risk.1<-unlist(lapply(as.list(risk.points),R.Weighted,error=x1,W=w1))
  
  risk.2<-round(risk.2)
  risk.1<-round(risk.1)
  
  # Note: "med" denotes general quantile, but median is default
  med.1<-suppressWarnings(min(at.points1[S1.KM<=quant]))
  # log transform
  KM1.lower<-exp(log(S1.KM)-1.96*SE1.KM/S1.KM)
  KM1.upper<-exp(log(S1.KM)+1.96*SE1.KM/S1.KM)
  med.1.lower<-suppressWarnings(min(at.points1[KM1.lower<=quant]))
  med.1.upper<-suppressWarnings(min(at.points1[KM1.upper<=quant]))
  
  if(details){
    cat("Treat=1: n,events=",c(length(x1),sum(eta1)),"\n")
    cat("Median, Lower, Upper=",c(med.1,med.1.lower,med.1.upper),"\n")
    cat("survfit medians","\n")
    print(round(summary(KM.fit1)$table,2))
  }
  
  med.2<-suppressWarnings(min(at.points2[S2.KM<=quant]))
  KM2.lower<-exp(log(S2.KM)-1.96*SE2.KM/S2.KM)
  KM2.upper<-exp(log(S2.KM)+1.96*SE2.KM/S2.KM)
  med.2.lower<-suppressWarnings(min(at.points2[KM2.lower<=quant]))
  med.2.upper<-suppressWarnings(min(at.points2[KM2.upper<=quant]))
  
  if(details){
    cat("Treat=0: n,events=",c(length(x2),sum(eta2)),"\n")
    cat("Median, Lower, Upper=",c(med.2,med.2.lower,med.2.upper),"\n")
    cat("survfit medians","\n")
    print(round(summary(KM.fit2)$table,2))
    
  }
  # Experimental
  fit1<-NA.CHR.Weighted(time=x1,Delta=(eta1==1),W.n=w1,W.d=w1,at.points=at.points1.plot)
  S1.KM.plot<-fit1$S.KM
  # Control
  fit2<-NA.CHR.Weighted(time=x2,Delta=(eta2==1),W.n=w2,W.d=w2,at.points=at.points2.plot)
  S2.KM.plot<-fit2$S.KM
  Y1.plot<-S1.KM.plot 
  Y2.plot<-S2.KM.plot
  
  Y1<-S1.KM
  Y2<-S2.KM
  
  if(is.null(x.truncate)) xmax<-max(c(at.points,xmax))
  if(!is.null(x.truncate)) xmax<-x.truncate
  
  if(plotit){
    plot(at.points1.plot,Y1.plot,type="s",ylim=c(ymin2,ymax),xlim=c(xmin,xmax),lty=ltys[1],col=col.1,lwd=lwds[1],xlab=Xlab,ylab=Ylab,axes=FALSE, cex.lab=cex_Yaxis)
    lines(at.points2.plot,Y2.plot,lty=ltys[2],type="s",col=col.2,lwd=lwds[2])
    
    if(show.ticks==TRUE){
      points(x1.cens,Y1[x1.match],pch=3,col=col.1,cex=censor.cex)
      points(x2.cens,Y2[x2.match],pch=3,col=col.2,cex=censor.cex)
    }
    # horizontal line just below ymin 
    abline(h=ymin-ymin.del,lty=1,col=1)
    
    # Include at least 2-points in-between ymin and ymax
    if(choose_ylim){
      if((ymax-ymin2) <= 0.5) prob.points <- round(seq(ymin,ymax,length=4),2)
    }
    if(!choose_ylim & is.null(prob.points)) prob.points <- seq(0,1,by=0.25)
    
    if(show.Y.axis) axis(2,at=c(prob.points),cex.axis=cex_Yaxis, las=1)
    
    axis(1,at=c(risk.points),labels=c(time.zero.label,risk.points[-1]),
         cex.axis=cex_Yaxis)
    
    box()
    if(risk.set){
      text(c(risk.points),c(y.risk2),c(risk.2),col=col.2,cex=risk.cex)
      text(c(risk.points),c(y.risk1),c(risk.1),col=col.1,cex=risk.cex)
    }
    if(add.segment) segments(-10,1,0,1)  
    
    if(show.med){
      if(med.1!=Inf & med.2!=Inf){
        text(xmax-xmed.offset,ymax-ymed.offset,paste(qlabel,round(med.1,1)),col=col.1,cex=med.cex)
        text(xmax-xmed.offset,ymax-(ymed.offset+del.med),paste(qlabel,round(med.2,1)),col=col.2,cex=med.cex)
      }
      
    }
    if(show.cox){
      rp<-vector('expression',4)
      rp[1]=substitute(expression(italic(hr)== MYVALUE1), 
                       list(MYVALUE1 = format(hr.ci[1],digits = cox.digits)))[2]
      rp[2]=substitute(expression(italic(lb) == MYVALUE2), 
                       list(MYVALUE2 = format(hr.ci[2], digits = cox.digits)))[2]
      rp[3]=substitute(expression(italic(ub) == MYVALUE3), 
                       list(MYVALUE3= format(hr.ci[3], digits = cox.digits)))[2]               
      rp[4]=substitute(expression(italic(p) == MYVALUE3), 
                       list(MYVALUE3= format.pval(pval.cox, eps=cox.eps, digits = cox.digits)))[2]               
      legend(put.legend.cox, legend = rp, bty = 'n',cex=cox.cex)
    }
    
    if(show.logrank){
      # Note, currently I do not have the stratified log-rank
      # So, if stratified, report 2-sided from survdiff
      if(!is.null(strata.name)){
        logrnk<-survdiff(Surv(Y,E)~Treat+strata(Strata))
        p.lr<-1-pchisq(logrnk$chisq,length(logrnk$n) - 1)
        
        rp<-vector('expression',1)
        rp[1]=substitute(expression(italic(p["2-sided"])== MYVALUE), 
                         list(MYVALUE =format.pval(p.lr,eps=lr.eps,digits=lr.digits)))[2]
        legend(put.legend.lr, legend = rp, bty = 'n',cex=logrank.cex)
      }
      
      if(is.null(strata.name)){
        rp<-vector('expression',1)
        rp[1]=substitute(expression(italic(p["1-sided"])== MYVALUE), 
                         list(MYVALUE =format.pval(p.lr,eps=lr.eps,digits=lr.digits)))[2]
        legend(put.legend.lr, legend = rp, bty = 'n',cex=logrank.cex)
      }
      
    }
  }
  
  if(show_arm_legend){
    
    legend("left",c(arms),col=c(col.1,col.2),lty=c(ltys),lwd=c(lwds),
           cex=arm.cex,bty="n")  
    
  }
  
  return(list(S1.KM=S1.KM,S2.KM=S2.KM,KM1.lower=KM1.lower,KM2.lower=KM2.lower,KM1.upper=KM1.upper,KM2.upper=KM2.upper,
              pval.cox=pval.cox,pval.lr=p.lr,cpoints=tpoints.add,
              S1.dpoints=S1.dpoints,SE1.dpoints=SE1.dpoints,
              S2.dpoints=S2.dpoints,SE2.dpoints=SE2.dpoints,
              LR.dpoints=LR.dpoints,ZLR.dpoints=ZLR.dpoints,
              at.points1=at.points1,at.points2=at.points2,
              Y1.plot=Y1.plot,Y2.plot=Y2.plot,
              at.points1.plot=at.points1.plot,at.points2.plot=at.points2.plot,
              med.1=med.1,med.2=med.2,events.1=sum(eta1),events.2=sum(eta2),
              y.risk1=y.risk1,y.risk2=y.risk2,
              risk.1=risk.1,risk.2=risk.2,risk.points=risk.points))
}




