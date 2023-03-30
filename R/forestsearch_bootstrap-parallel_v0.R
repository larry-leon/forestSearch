
# Get operator (parallel or serial)
get_operator <- function() {
  is_par <- foreach::getDoParWorkers() > 1
  if (is_par) {
    res <- foreach::`%dopar%`
  } else {
    res <- foreach::`%do%`
  }
  res
}

#' Parallelize arbitrary user-supplied expression
#'
#' @param expr R expression.
#' @param n Number of iterations.
#' @param seed Initial random seed.
#' @param counter Character string. Name of the loop counter variable. Default is `"i"`.
#' @param export Exported variables to be used in `expr`.
bootPar <- function(expr, boots, seed, counter = "i", export) {
  `%op%` <- get_operator()
  expr <- substitute(expr)
  res <- foreach::foreach(
    i = seq_len(boots),
    .export = export,
    .packages = c("survival","grf","policytree","aVirtualTwins","randomForest","SPlit","data.table","cubature","plyr"),
    .combine = "rbind.fill",
    .errorhandling = "pass"
  ) %op% {
    assign(counter, value = i)
    set.seed(seed + i - 1)
    eval(expr)
  }
  res
}

# use rbind here
bootYstar <- function(expr, boots, seed, counter = "i", export) {
  `%op%` <- get_operator()
  expr <- substitute(expr)
  res <- foreach::foreach(
    i = seq_len(boots),
    .export = export,
    .packages = c("survival","grf","policytree","aVirtualTwins","randomForest","SPlit","data.table","cubature"),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %op% {
    assign(counter, value = i)
    set.seed(seed + i - 1)
    eval(expr)
  }
  res
}


fun_arg_list_boot2<-c("dgm","maxFollow","sim_aftm1_gbsg","forestsearch","fs.estimates.out","vi.grf.min","dummy","subgroup.search",
                "subgroup.consistency","acm.disjctif","max.minutes","ztrail","one.zero","maxk","nmin.fs","d.min","hr.threshold",
                "NA.CHR.Weighted","R.Weighted","N.Weighted","hr.consistency","stop.threshold","fs.splits","m1.threshold","pconsistency.threshold",
                "confounders.name","outcome.name","id.name","treat.name",
                "cox.formula.sim","cox.formula.adj.sim",
                "event.name",
                "get.FG","FG.transform","LS.int.FG","vt.subg.harm.survival","vt.threshold","treat.threshold","maxdepth","n.min","ntree","vt.data2",
                "formatRCTDataset2","vt.estimates.out","label.analyses","grf.subg.harm.survival","grf.estimates.out","dmin.grf",
                "frac.tau","sg_focus",
                "get.FS","get.VT","get.GRF","grf.fs.subg","muC.adj","oc_analyses","N","split_method","quiet",
                "cox.formula.boot","forestsearch_bootstrap","oc_analyses_FSboot","fsBoot.estimates.out","boots","est.loghr",
                "get_Ystar","count.id","fsboot_forparallel","get_Cox_sg")


fun_arg_list_boot<-c("forestsearch","fs.estimates.out","vi.grf.min","dummy","subgroup.search",
                      "subgroup.consistency","acm.disjctif","max.minutes","ztrail","one.zero","maxk","nmin.fs","d.min","hr.threshold",
                      "NA.CHR.Weighted","R.Weighted","N.Weighted","hr.consistency","stop.threshold","fs.splits","m1.threshold","pconsistency.threshold",
                      "confounders.name","outcome.name","id.name","treat.name",
                      "event.name","sg_focus","pstop_futile",
                      "get.FG","FG.transform","LS.int.FG","vt.subg.harm.survival",
                      "formatRCTDataset2","grf.subg.harm.survival","grf.estimates.out","dmin.grf","n.min","frac.tau",
                      "grf.fs.subg","split_method","quiet",
                      "cox.formula.boot","forestsearch_bootstrap","oc_analyses_FSboot","fsBoot.estimates.out","boots","est.loghr",
                      "get_Ystar","count.id","fsboot_forparallel","get_Cox_sg","df_boot_analysis","fs.est","H_obs","Hc_obs")


get_Ystar<-function(boot){
NN<-nrow(df_boot_analysis)
id0<-c(1:NN)
set.seed(8316951+boot*1000)
id_boot<-sample(id0,size=NN,replace=TRUE)
df_boot<-df_boot_analysis[id_boot,]
# re-name id.boot for bootstrap data
df_boot$id_boot<-c(1:nrow(df_boot))
Ystar<-c(unlist(lapply(df_boot_analysis$id,count.id,dfb=df_boot)))
return(c(Ystar))
}



fsboot_forparallel<-function(boot,cox2.formula.boot=as.formula(paste("Surv(Y,Event)~Treat"))){
  #Check names of df and df_boot
  # For simulations remove these checks
  #names_tocheck<-all.vars(cox.formula.boot)
  #check<-unlist(lapply(names_tocheck,grep,names(df_boot_analysis),value=TRUE))
  #check2<-match(names_tocheck,check)
  #if(sum(!is.na(check2))!=length(names_tocheck)) stop("Dataset df does NOT contain cox.formula.boot variables")
  # Create "id" variable if not included
  #if(!any(names(df_boot_analysis)=="id")) df_boot_analysis$id<-c(1:nrow(df_boot_analysis))
  # In forestsearch algorithm, the Cox model variables are generically as
  #cox2.formula.boot<-as.formula(paste("Surv(Y,Event)~Treat"))
  
  NN<-nrow(df_boot_analysis)
  
  # if(is.null(fs.est$grp.consistency$sg.harm)) stop("Forest search subgroup not found")
  
  grp.consistency<-fs.est$grp.consistency  
  sg.harm.fs<-grp.consistency$sg.harm
  
  pstop_futile<-fs.est$pstop_futile
  sg_focus<-fs.est$sg_focus
  
  # Observed data evaluated at H:   hrH_obs
  # Observed data evaluated at H*   hrHstar_obs
  # Bootstrap data evaluated at H:  hrH_star 
  # Bootstrap data evaluated at H*: hrHstar_star
  
  id0<-c(1:NN)
  set.seed(8316951+boot*1000)
  id_boot<-sample(id0,size=NN,replace=TRUE)
  df_boot<-df_boot_analysis[id_boot,]
  # re-name id.boot for bootstrap data
  df_boot$id_boot<-c(1:nrow(df_boot))
 
  ############################################
  # Bootstrap data evaluated at H:  H_star 
  ############################################
  
  fitH_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  H_star<-fitH_star$est_obs
  
  fitHc_star<-get_Cox_sg(df_sg=subset(df_boot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  Hc_star<-fitHc_star$est_obs
  
  # Bias corrections
  H_biasadj_1<-H_biasadj_2<-NA
  Hc_biasadj_1<-Hc_biasadj_2<-NA
  tmins_search<-NA
  max_sg_est<-NA
  
  # Reset treat.recommend
  # df and df_boot have observed treat.recommend estimates
  # Bootstrap versions will create bootstrap version of treat.recommend
  # Can also just drop
  drop.vars<-c("treat.recommend")
  dfnew<-df_boot_analysis[ ,!(names(df_boot_analysis) %in% drop.vars)]
  dfnew_boot<-df_boot[ ,!(names(df_boot) %in% drop.vars)]
  
  
  tempB<-try(forestsearch(df=dfnew_boot,confounders.name=confounders.name,
                          df.predict=dfnew,
                          split_method=split_method,
                          details=FALSE,
                          outcome.name=outcome.name,treat.name=treat.name,event.name=event.name,
                          id.name=id.name,
                          n.min=nmin.fs,hr.threshold=hr.threshold,hr.consistency=hr.consistency,fs.splits=fs.splits,
                          stop.threshold=stop.threshold,d0.min=d.min,d1.min=d.min,
                          pconsistency.threshold=pconsistency.threshold,
                          max.minutes=max.minutes,
                          sg_focus=sg_focus,pstop_futile=pstop_futile,
                          maxk=maxk,plot.sg=FALSE,vi.grf.min=vi.grf.min),TRUE)
  
  
  # FS with error output NA's
  # FS done (w/o error) but NO sg found
  if(!inherits(tempB,"try-error") & is.null(tempB$sg.harm)){
    
    tmins_search<-tempB$find.grps$time_search
    
    if(!is.null(tempB$find.grps$max_sg_est)){
      max_sg_est<-tempB$find.grps$max_sg_est
    }
  }
  # FS done AND sg found
  if(!inherits(tempB,"try-error") & !is.null(tempB$sg.harm)){# Boot estimable & sg found 
    
    df_PredBoot<-tempB$df.pred # Observed data at bootstrap estimate of H
    
    dfboot_PredBoot<-tempB$df.est # Bootstrap data at H* (bootstrap estimate of H)
    
    # names_tocheck<-all.vars(cox2.formula.boot)
    #  check<-unlist(lapply(names_tocheck,grep,names(dfboot_PredBoot),value=TRUE))
    # check2<-match(names_tocheck,check)
    #if(sum(!is.na(check2))!=length(names_tocheck)) stop("Bootstrap dataset (dfboot_PredBoot) df does NOT contain cox2.formula.boot variables")
    
    sg.harm.boot<-tempB$sg.harm
    
    max_sg_est<-tempB$find.grps$max_sg_est
    tmins_search<-tempB$find.grps$time_search
    
    rm("tempB")
    
    ##############################################
    # Observed data evaluated at H*   hrHstar_obs
    ##############################################
    fitHstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
    Hstar_obs<-fitHstar_obs$est_obs
    rm("fitHstar_obs")
    fitHstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==0),cox.formula=cox2.formula.boot,est.loghr=est.loghr)
    Hstar_star<-fitHstar_star$est_obs
    rm("fitHstar_star")
    
    H_biasadj_1<-H_obs-(Hstar_star-Hstar_obs)
    H_biasadj_2<-2*H_obs-(H_star+Hstar_star-Hstar_obs)
    
    # Hc
    fitHcstar_obs<-get_Cox_sg(df_sg=subset(df_PredBoot,treat.recommend==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
    Hcstar_obs<-fitHcstar_obs$est_obs
    rm("fitHcstar_obs")
    
    fitHcstar_star<-get_Cox_sg(df_sg=subset(dfboot_PredBoot,treat.recommend==1),cox.formula=cox2.formula.boot,est.loghr=est.loghr)
    Hcstar_star<-fitHcstar_star$est_obs
    rm("fitHcstar_star")
    
    Hc_biasadj_1<-Hc_obs-(Hcstar_star-Hcstar_obs)
    Hc_biasadj_2<-2*Hc_obs-(Hc_star+Hcstar_star-Hcstar_obs)
    
  } # Bootstrap estimable
  
  dfres<-data.table(H_biasadj_1,H_biasadj_2,
                    Hc_biasadj_1,Hc_biasadj_2,
                    tmins_search,max_sg_est)
  return(dfres)
}




# On original scale of est.loghr
# Does NOT convert 
ci_est<-function(x,sd,alpha=0.025,scale="hr"){
  c_alpha<-qnorm(1-alpha)
  c_low<-x-c_alpha*sd
  c_up<-x+c_alpha*sd
  est<-x
  if(scale=="hr" & est.loghr){
    # Only convert if log(hr) scale
    est<-exp(x)
    sd<-exp(x)*sd
    c_low<-exp(c_low)
    c_up<-exp(c_up)
  }
  length<-c_up-c_low
  return(list(length=length,lower=c_low,upper=c_up,sd=sd,est=est))
}


# For code checking
#H1_adj<-resB$H_biasadj_1
#ystar<-Ystar_mat

get_dfRes<-function(Hobs,seHobs,H1_adj,H2_adj=NULL,ystar,cov_method="standard",cov_trim=0.0,est.scale="hr"){
# Un-adjusted
cest<-ci_est(x=Hobs,sd=seHobs,scale=est.scale)
H0_lower<-cest$lower
H0_upper<-cest$upper
sdH0<-cest$sd
H0<-cest$est
rm("cest")

est<-get_targetEst(x=H1_adj,ystar=ystar,cov_method=cov_method,cov_trim=cov_trim)
q<-est$target_est
se_new<-est$sehat_new
rm("est")

cest<-ci_est(x=q,sd=se_new,scale=est.scale)
H1_lower<-cest$lower
H1_upper<-cest$upper
sdH1<-cest$sd
H1<-cest$est
rm("cest","q","se_new")

outres<-data.table(H0,sdH0,H0_lower,H0_upper,
                   H1,sdH1,H1_lower,H1_upper)

if(!is.null(H2_adj)){
est<-get_targetEst(x=H2_adj,ystar=ystar,cov_method=cov_method,cov_trim=cov_trim)
q<-est$target_est
se_new<-est$sehat_new
rm("est")

cest<-ci_est(x=q,sd=se_new,scale=est.scale)
H2_lower<-cest$lower
H2_upper<-cest$upper
sdH2<-cest$sd
H2<-cest$est
rm("cest","q","se_new")

outres<-data.table(H0,sdH0,H0_lower,H0_upper,
                H1,sdH1,H1_lower,H1_upper,
                H2,sdH2,H2_lower,H2_upper)
}

return(outres)
}



fsBoot.parallel.out<-function(df,cox.formula.boot,FSboots,est.loghr=TRUE,sim_number,est.scale="hr"){
  
  fit<-get_Cox_sg(df_sg=df,cox.formula=cox.formula.sim,est.loghr=est.loghr)
  H_itt<-fit$est_obs
  seH_itt<-fit$se_obs
  
  cest<-ci_est(x=H_itt,sd=seH_itt,scale=est.scale)
  H_itt<-cest$est
  seH_itt<-cest$sd
  H_itt_lower<-cest$lower
  H_itt_upper<-cest$upper
  rm("cest")
  rm("fit")
  
  fit<-get_Cox_sg(df_sg=df,cox.formula=cox.formula.adj.sim,est.loghr=est.loghr)
  H_adj_itt<-fit$est_obs["treat"]
  seH_adj_itt<-fit$se_obs["treat"]
  
  cest<-ci_est(x=H_adj_itt,sd=seH_adj_itt,scale=est.scale)
  H_adj_itt<-cest$est
  seH_adj_itt<-cest$sd
  H_adj_itt_lower<-cest$lower
  H_adj_itt_upper<-cest$upper
  rm("cest")
  rm("fit")
  
  # hatH(H true)
  fitH<-get_Cox_sg(df_sg=subset(df,flag.harm==1),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  H_true<-fitH$est_obs
  seH_true<-fitH$se_obs
  
  if(est.scale=="hr" & est.loghr){
  H_true<-exp(H_true)
  seH_true<-H_true*seH_true
  }
  
  sizeH_true<-with(df,sum(flag.harm))
  propH_true<-with(df,mean(flag.harm))
  
  fitHc<-get_Cox_sg(df_sg=subset(df,flag.harm==0),cox.formula=cox.formula.boot,est.loghr=est.loghr)
  Hc_true<-fitHc$est_obs
  seHc_true<-fitHc$se_obs
  
  if(est.scale=="hr" & est.loghr){
  Hc_true<-exp(Hc_true)
  seHc_true<-Hc_true*seHc_true
  }  
  rm("fitH","fitHc")
  
  anyH_found<-any(names(df)=="treat.recommend")
  
  if(!anyH_found){
    # Causal hat(H)
    hatH_causal<-NA
    # Causal hat(Hc)
    hatHc_causal<-NA
  }
  
  if(is.null(FSboots) | !anyH_found){
    # initiate as NA
    H_obs<-seH_obs<-H1.bc<-seH1.bc<-NA
    H2.bc<-seH2.bc<-NA
    Hc_obs<-seHc_obs<-Hc1.bc<-seHc1.bc<-NA
    Hc2.bc<-seHc2.bc<-NA
    H0_cover1<-H0_cover2<-H0_cover3<-H0_length<-NA
    H1_cover1<-H1_cover2<-H1_cover3<-H1_length<-NA
    H2_cover1<-H2_cover2<-H2_cover3<-H2_length<-NA
    Hc0_cover1<-Hc0_cover2<-Hc0_cover3<-Hc0_length<-NA
    Hc1_cover1<-Hc1_cover2<-Hc1_cover3<-Hc1_length<-NA
    Hc2_cover1<-Hc2_cover2<-Hc2_cover3<-Hc2_length<-NA
  }
  
  
  if(anyH_found){
    # Causal hat(H)
    hr_1_hatH<-with(subset(df,treat.recommend==0),mean(h1.potential))
    hr_0_hatH<-with(subset(df,treat.recommend==0),mean(h0.potential))
    hatH_causal<-hr_1_hatH/hr_0_hatH
    
    # Causal hat(Hc)
    hr_1_hatHc<-with(subset(df,treat.recommend==1),mean(h1.potential))
    hr_0_hatHc<-with(subset(df,treat.recommend==1),mean(h0.potential))
    hatHc_causal<-hr_1_hatHc/hr_0_hatHc
  
    H_causal<-dgm$hr.H.true
    Hc_causal<-dgm$hr.Hc.true
    
    if(est.loghr & est.scale=="loghr"){
    hatH_causal<-log(hatH_causal)
    hatHc_causal<-log(hatHc_causal)
    H_causal<-log(dgm$hr.H.true)
    Hc_causal<-log(dgm$hr.Hc.true)
    }
  
    }
  
  
  if(!is.null(FSboots) & anyH_found){
    
    # Estimates are on est.scale (converted to HR if est.loghr=TRUE & est.scale="hr")
    
    dfH<-FSboots$H_estimates
   
    H_obs<-dfH$H0
    seH_obs<-dfH$sdH0
    
    H1.bc<-dfH$H1
    seH1.bc<-dfH$sdH1
  
    H2.bc<-dfH$H2
    seH2.bc<-dfH$sdH2
    
    # Un-adjusted CI
    cest<-ci_cover(lower=dfH$H0_lower,upper=dfH$H0_upper,target=c(H_true,hatH_causal,H_causal))
    H0_cover1<-cest$cover[1]
    H0_cover2<-cest$cover[2]
    H0_cover3<-cest$cover[3]
    H0_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfH$H1_lower,upper=dfH$H1_upper,target=c(H_true,hatH_causal,H_causal))
    H1_cover1<-cest$cover[1]
    H1_cover2<-cest$cover[2]
    H1_cover3<-cest$cover[3]
    H1_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfH$H2_lower,upper=dfH$H2_upper,target=c(H_true,hatH_causal,H_causal))
    H2_cover1<-cest$cover[1]
    H2_cover2<-cest$cover[2]
    H2_cover3<-cest$cover[3]
    H2_length<-cest$LC
    rm("cest")
    rm("dfH")
    
    # Hc
    dfHc<-FSboots$Hc_estimates
    
    Hc_obs<-dfHc$H0
    seHc_obs<-dfHc$sdH0
    
    Hc1.bc<-dfHc$H1
    seHc1.bc<-dfHc$sdH1
    
    Hc2.bc<-dfHc$H2
    seHc2.bc<-dfHc$sdH2
    
    # Un-adjusted CI
    cest<-ci_cover(lower=dfHc$H0_lower,upper=dfHc$H0_upper,target=c(Hc_true,hatHc_causal,Hc_causal))
    Hc0_cover1<-cest$cover[1]
    Hc0_cover2<-cest$cover[2]
    Hc0_cover3<-cest$cover[3]
    Hc0_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfHc$H1_lower,upper=dfHc$H1_upper,target=c(Hc_true,hatHc_causal,Hc_causal))
    Hc1_cover1<-cest$cover[1]
    Hc1_cover2<-cest$cover[2]
    Hc1_cover3<-cest$cover[3]
    Hc1_length<-cest$LC
    rm("cest")
    
    cest<-ci_cover(lower=dfHc$H2_lower,upper=dfHc$H2_upper,target=c(Hc_true,hatHc_causal,Hc_causal))
    Hc2_cover1<-cest$cover[1]
    Hc2_cover2<-cest$cover[2]
    Hc2_cover3<-cest$cover[3]
    Hc2_length<-cest$LC
    rm("cest")
    rm("dfHc")
      }
  
  df.res<-data.table(sim_number,sizeH_true,propH_true,
                     H_true,Hc_true,
                     hatH_causal,hatHc_causal,
                     H_obs,H1.bc,H2.bc,
                     Hc_obs,Hc1.bc,Hc2.bc,
                     seH_true,seHc_true,
                     seH_obs,seHc_obs,
                     seH1.bc,seH2.bc,
                     seHc1.bc,seHc2.bc,
                     H0_cover1,H0_cover2,H0_cover3,H0_length,
                     H1_cover1,H1_cover2,H1_cover3,H1_length,
                     H2_cover1,H2_cover2,H2_cover3,H2_length,
                     Hc0_cover1,Hc0_cover2,Hc0_cover3,Hc0_length,
                     Hc1_cover1,Hc1_cover2,Hc1_cover3,Hc1_length,
                     Hc2_cover1,Hc2_cover2,Hc2_cover3,Hc2_length,
                     H_itt,seH_itt,H_itt_lower,H_itt_upper,
                     H_adj_itt,seH_adj_itt,H_adj_itt_lower,H_adj_itt_upper)
  return(df.res)
}


