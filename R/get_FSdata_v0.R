get_defaultparms <- FALSE
if(get_defaultparms){
df <- x
grf_cuts<-NULL
use_grf <- FALSE
use_lasso <- FALSE
confounders.name <- Allconfounders.name
conf_force<-NULL
conf.cont_medians<-NULL
conf.categorical_only <- NULL
conf.cont_medians_force <- NULL
cont.cutoff <- 4
}


# Wrapper for setting up dataset for ForestSearch
# which includes options to apply lasso for dimension reduction
# and GRF cuts

# Note: make cutoff an argument in wrapper
is.continuous<-function(x,cutoff=4){ifelse(length(unique(x))>=cutoff,1,2)}

#conf.categorical<-c("hemo","homo","drugs","race","gender","oprior","symptom")
# Force inclusion of categorical factors
#conf.cat_force<-c("gender") 
# Force inclusion of continuous cut at medians
#conf.cont_medians_force<-c("age")


get_FSdata<-function(df,use_lasso=FALSE,use_grf=FALSE,use_grf_only=FALSE,grf_cuts,confounders.name,
cont.cutoff=4,conf_force=NULL, conf.cont_medians=NULL, conf.cont_medians_force=NULL,
outcome.name=NULL, event.name=NULL, details=TRUE){
  
  if(use_grf & is.null(grf_cuts)){
    warning("GRF cuts need to be specified (run prior): Considering no cuts per GRF turning off GRF")
    use_grf<-FALSE
  }  
  if(use_lasso &  (is.null(outcome.name) | is.null(event.name))) stop("Cox variable names needed for lasso")
  
  # Flag continuous and categorical (continuous if unique values exceeds con.cutoff)
  flag_continuous<-rep(NA,length(confounders.name))
  for(ccs in 1:length(confounders.name)){
    aa<-df[,c(confounders.name[ccs])]
    flag_continuous[ccs]<-(is.continuous(aa,cutoff=cont.cutoff)==1)  
  }
  
if(details){
cat("# of continuous/categorical characteristics",c(sum(flag_continuous),sum(!flag_continuous)),"\n")
if(sum(flag_continuous)>0)  cat("Continuous characteristics:",c(confounders.name[flag_continuous]),"\n")
if(sum(!flag_continuous)>0)   cat("Categorical characteristics:",c(confounders.name[!flag_continuous]),"\n")
}

  if(sum(flag_continuous)==0) conf.categorical <- confounders.name
  if(sum(flag_continuous) > 0) conf.categorical <- confounders.name[!flag_continuous]
  
  if(is.null(conf.cont_medians) & is.null(conf.cont_medians_force) & sum(flag_continuous)>0){
  conf.cont_medians<-confounders.name[flag_continuous]
  }  
  
  # Note: If lasso does not select any factors, then defaults to 
  # splitting continuous at medians and all categorical
  # However, if use_grf then median cuts will be replaced by grf cuts
  # and keep all categorical
  
  lassokeep <- NULL
  lassoomit <- NULL
  
  if(use_lasso){
    set.seed(8316951)
    # Reduce dimension via Cox lasso
    xx<-as.matrix(df[,confounders.name])
    yy<-as.matrix(df[,c(outcome.name,event.name)])
    colnames(yy)<-c("time","status")
    
    cvfit <- cv.glmnet(xx, yy, family="cox") #first do 10-fold cross-validation to select lambda
    
  if(details) cat("CV lambda =",c(cvfit$lambda.min),"\n")
    
    m <- glmnet(xx, yy, family="cox", lambda=cvfit$lambda.min) #plugin the optimal lambda
    
    lassokeep<-confounders.name[which(m$beta!=0)]
    lassoomit<-confounders.name[which(m$beta==0)]
    
    if(details){
      print(m$beta)
      cat("Cox-LASSO selected:",c(lassokeep),"\n")
      cat("Cox-LASSO not selected:",c(lassoomit),"\n")
    }
    
    
if(length(lassokeep)>0){
   if(!is.null(conf.cont_medians)){
      to_keep<-which(conf.cont_medians %in% lassokeep)
      if(length(to_keep)>0) conf.cont_medians <- conf.cont_medians[to_keep]
      if(length(to_keep)==0) conf.cont_medians <- NULL
      }
    if(!is.null(conf.categorical)){
      to_keep<-which(conf.categorical %in% lassokeep)
      if(length(to_keep)>0)  conf.categorical <- conf.categorical[to_keep]
      if(length(to_keep)==0) conf.categorical <- NULL
      }
    } # If any selected per Lasso
  } # Done Lasso
  
  
  if(use_grf_only) conf.cont_medians <- NULL
  # If forcing any median cuts, then done below
  
if(details & use_lasso){ 
  cat("Median cuts after Lasso:",c(conf.cont_medians),"\n")
  cat("Categorical after Lasso:",c(conf.categorical),"\n")
}  
    
  if(!is.null(conf.cont_medians_force)){
    conf.cont_medians<-c(conf.cont_medians,conf.cont_medians_force)
  }
  

if(details & use_grf){
cat("Factors per GRF:",c(grf_cuts),"\n")
}    
  
if(details & use_grf & length(conf.cont_medians)>0){
cat("Medians prior to removing if also in GRF:",c(conf.cont_medians),"\n")
}    
  
  
# Remove any factors to cut at median if already in GRF
if(use_grf & length(conf.cont_medians)>0){

if(length(grf_cuts)>0){  
    # Which to omit  
    cMed_flag<-NULL
    for(ccs in 1:length(grf_cuts)){
    thiscut<-grf_cuts[ccs]
    # confs = c(x1,x2,x3)
    # eg "x1 <= 78"
    # "x3 <= 43"
    flag_omit<-c(unlist(lapply(conf.cont_medians,function(x){grepl(x,thiscut)})))
    # first ccs
    # flag_omit = T,F,F
  if(any(flag_omit))  cMed_flag<-c(cMed_flag,c(conf.cont_medians[flag_omit]))
   }  


if(length(cMed_flag)>0){
to_exclude <- (conf.cont_medians %in% cMed_flag)

# If excluding at least 1, the following remain
if(sum(!to_exclude)>0)  conf.cont_medians <- conf.cont_medians[!to_exclude]

# If excluding all then set to null
if(sum(!to_exclude)==0)  conf.cont_medians <- NULL

if(details){
cat_line("***cMed_flag***=",c(cMed_flag),col="red")
cat_line("***to_exclude***=",c(to_exclude),col="red")
cat_line("***conf.cont_medians***=",c(conf.cont_medians),col="orange")
}
      } # Otherwise, NO factors omitted
     }
  } 
  
if(details & use_grf & length(conf.cont_medians)>0)  cat("Factors after removing any duplicates also in GRF:",c(conf.cont_medians),"\n")
  
if(details & use_grf & length(conf.cont_medians)==0){
cat_line("***conf.cont_medians is NULL --> NO MEDIAN CUTS***",col="purple")
}  
  
# Re-introduce conf.cont_force_medians

if(!is.null(conf.cont_medians_force)) conf.cont_medians <- c(conf.cont_medians,conf.cont_medians_force)

  conf.cont_Medcuts<-NULL

    if(!is.null(conf.cont_medians)){
    # Default is to specify median cuts
    # Continuous covariate indexes
    for(aa in 1:length(conf.cont_medians)){
      cc1<-paste0(conf.cont_medians[aa]," <= ",sep="")  
      cc2<-paste0("median","(", sep="")
      cc3<-paste0(cc2,conf.cont_medians[aa],sep=")")
      cc4<-paste0(cc1,cc3,sep="")
      conf.cont_Medcuts<-c(conf.cont_Medcuts,cc4)
    }
  }
  
  confs<-c(conf.categorical,conf.cont_Medcuts)

  if(details & use_grf & length(confs)>0){
    cat_line("***Factors per lasso after omitting GRF dups***=",c(confs),col="orange")
  }  
  
  
  # At this stage, these are confs per Lasso (GRF step is next)
  if(use_lasso) confs_lasso <- confs 

if(use_grf){
  # If element in grf_cuts matches conf.categorical
  # then this is redundant per FS handling of categorical
  # factors (Any level will already be incorporated)
  # So remove grf_cuts already in conf.categorical
  if(details)  cat("Initial GRF cuts included",c(grf_cuts),"\n")  
  if(length(grf_cuts)>0){  
  grf_cuts_keep<-NULL
  for(ccs in 1:length(grf_cuts)){
   thiscut<-grf_cuts[ccs]
  # Covariate corresponding to
  flag_omit<-c(unlist(lapply(confs,function(x){grepl(x,thiscut)})))

  # Only keep if does not already exist in conf.categorical 
  if(all(!flag_omit)) grf_cuts_keep<-c(grf_cuts_keep,thiscut)
   }  
  confs <- unique(c(confs,grf_cuts_keep))
  }
  } # use_grf
  
  # Factors include per GRF not in confs_lasso

if(use_lasso & use_grf){
which_both <- (confs %in% confs_lasso)
if(details & length(which_both) >0){
cat("Factors included per GRF (not in lasso)",c(confs[!which_both]),"\n")
}
  }
  
 
  if(!is.null(conf_force)) confs<-unique(c(confs,conf_force)) 
  
  n_confs<-length(confs)
  if(n_confs==0) stop("Error in FS dataset prior to flag drop")
  
  # Arrange by continuous and categorical
  # Classify as continuous or categorical
  flag_continuous<-rep(NA,length(confs))
  flag_drop<-rep(NA,length(confs)) # Flag factors with only 1 level and drop
  for(ccs in 1:length(confs)){
    thiscut<-confs[ccs]
    # Covariate corresponding to
    cov_index<-c(unlist(lapply(confounders.name,function(x){grepl(x,thiscut)})))
    cut_name<-confounders.name[cov_index]
    aa<-df[,c(cut_name)]
    flag_continuous[ccs]<-(is.continuous(aa,cutoff=cont.cutoff)==1)  
    flag_drop[ccs]<-(length(unique(aa))<=1)

  }
  
  if(details & any(flag_drop)){
    cat("Dropping variables (flag_drop):",c(confs[flag_drop]),"\n")
  }
  
  conf.categorical<-confs[!flag_continuous & !flag_drop]
  
  conf.cont_cuts <- NULL
  if(sum(flag_continuous)>0){
  conf.cont_cuts<-confs[flag_continuous & !flag_drop]
  }
  
  # Note, possibility that a cut point only has 1 level
  # Drop any such factors

flag_drop2<-NULL    
if(length(conf.cont_cuts)>=1){
  # Note, need to copy (otherwise qcheck will not be saved) 
    df_check<-copy(df)
  flag_drop2<-rep(NA,length(conf.cont_cuts))
    for(ccs in 1:length(conf.cont_cuts)){
      thiscut<-conf.cont_cuts[ccs]
      df_check <- within(df_check,{
        qcheck <- as.factor(ifelse(eval(parse(text=thiscut)),1,0))
      })
      flag_drop2[ccs]<-(length(unique(df_check$qcheck))<=1)
    }
  
  if(details & any(flag_drop2)){
    cat("Dropping variables (flag_drop2)",c(conf.cont_cuts[flag_drop2]),"\n")
  }
conf.cont_cuts<-conf.cont_cuts[!flag_drop2]
rm("df_check")
  } # flag_drop2   
  
# n_confs after drop1 and drop2
n_confs<-n_confs-(sum(flag_drop)+sum(flag_drop2))

confs<-c(conf.cont_cuts,conf.categorical)

if(length(confs)!=n_confs) stop("Error in categorical and continuous classifications")
  
  # These will be names q1,..,qK (say)
  names_new<-c(unlist(lapply(c(1:length(confs)),function(x){paste0("q",x,sep="")})))
  
  # For the continuous
  if(length(conf.cont_cuts)>=1){
    for(ccs in 1:length(conf.cont_cuts)){
      thiscut<-conf.cont_cuts[ccs]
      df<-within(df,{
        q <- as.factor(as.numeric(ifelse(eval(parse(text=thiscut)),1,0)))
      })
      loc_names<-which(names(df)=="q")
      names(df)[loc_names]<-names_new[ccs]
    }
  }
  
  if(length(conf.categorical)>=1){
    for(ccs in 1:length(conf.categorical)){
      thiscut<-conf.categorical[ccs]
      df<-within(df,{
        q <- as.factor(as.numeric(eval(parse(text=thiscut))))
      })
      loc_names<-which(names(df)=="q")
      names(df)[loc_names]<-names_new[ccs+length(conf.cont_cuts)] # Add cont index 
    }
  }

# NOTE: include check that all confs are 0,1 (not TRUE, FALSE)
#check_factors<-apply(df[,c(names_new)],2,function(x){c(unique(x))})
#if(!all(c(check_factors) %in% c(0,1))) stop("Error in factor setup")
 
check_factors<-apply(df[,c(names_new)],2,function(x){c(unique(x))})
if(any(c(check_factors) %in% c("TRUE","FALSE"))) stop("Error in factor setup")
  
  if(details){
    cat_rule(col="orange")
    cat("# of candidate subgroup factors=",c(length(confs)),"\n")
    print(confs)
    cat_rule(col="orange")
    if(any(flag_drop) | any(flag_drop2)){
    cat("At least 1 variable dropped","\n")
    cat_line("***A variable was dropped (check)***=",col="red")
    }

  }  
  
  return(list(df=df,confs_names=names_new,confs=confs,lassokeep=lassokeep,lassoomit=lassoomit))
}




FS_labels<-function(Qsg,confs_labels){
  label_index<-c(1:length(confs_labels))
  # Check how many elements
  # More than 4 --> factor is 10th or greater 
  flag_4<-as.numeric(substr(Qsg,5,5))
  if(is.na(flag_4)){
    # extract index from Qsg
    qindx<-as.numeric(substr(Qsg,2,2))
    qaction<-as.numeric(substr(Qsg,4,4))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  if(!is.na(flag_4)){
    # extract index from Qsg
    qindx<-as.numeric(substr(Qsg,2,3))
    qaction<-as.numeric(substr(Qsg,5,5))
    sg_label<-confs_labels[qindx]
    if(qaction==0) sg_label <- paste0("!{",sg_label,sep="}")
    if(qaction==1) sg_label <- paste0("{",sg_label,sep="}")
  }
  return(sg_label)
}
