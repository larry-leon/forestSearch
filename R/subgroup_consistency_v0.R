
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 


#if(all(!is.na(hr.subgroups$m1))) found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold & hr.subgroups$m1<=m1.threshold,]
#if(any(is.na(hr.subgroups$m1))) found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold,]

# sg_focus = hr --> sort by HR's to prioritize largest HR's
# sg_focus= Nsg ---> prioritize SG size

subgroup.consistency<-function(df,hr.subgroups,hr.threshold=1.0,hr.consistency=1.0,pconsistency.threshold=0.5,m1.threshold=Inf,n.splits=100,details=FALSE,
sg_focus="hr",stop.threshold=0.95,by.risk=2,plot.sg=FALSE,split_method="Random",maxk=7,pstop_futile=0.5){

if(!(sg_focus %in% c("hr","Nsg","NsgOnly","hrOnly"))) stop("sg_focus either hr or Nsg in subgroup.consistency function")

 names.Z<-c(names(hr.subgroups)[-c(1:6)])
  
 # Only implement m1.threshold criteria if ALL sg's treatment medians are estimable
  
  if(m1.threshold==Inf) found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold,]  
  # Implement m1.threshold criteria for sg's where treatment medians are estimable
  if(m1.threshold!=Inf){
  hr.subgroups<-hr.subgroups[which(!is.na(hr.subgroups$m1)),]
  found.hrs<-hr.subgroups[hr.subgroups$HR>=hr.threshold & hr.subgroups$m1<=m1.threshold,]
  }

# sort by HR 
if(sg_focus %in% c("hr","hrOnly")) found.hrs<-found.hrs[order(found.hrs$HR,decreasing=TRUE),]
# Sort by SG sample size if sg_focus="Nsg"
if(sg_focus %in% c("Nsg","NsgOnly"))  found.hrs<-found.hrs[order(found.hrs$n,decreasing=TRUE),]
  
if(details){
cat("Subgroups (1st 10) meeting overall screening thresholds (HR, m1) sorted by focus: (m1,sg_focus)=",c(m1.threshold,sg_focus),"\n")
  print(round(found.hrs[1:10,-c(1)],2))
}  
  n.groups<-0 
  index.Z<-found.hrs[,-c(1:6)]
  
  result_new<-NULL
  result<-NULL
  this.1<-NULL
  # Store group id and consistency %
  any.found<-0
  p.consistency<-0
  stop_futile<-FALSE
    # Count failure threshold
  count_failure<-round((1-pconsistency.threshold)*n.splits,1)
  
    for(m in 1:nrow(found.hrs)){
  # Stopping at first sg exceeding stop.threshold
  # Or stopping after first sg below pstop_futile
  # Since HR's are sorted, chance of meeting 90% likely "futile"
    if(p.consistency<stop.threshold & !stop_futile){
      
      indexm<-as.numeric(unlist(index.Z[m,]))
      this.m<-names.Z[indexm==1]
      n.groups<-n.groups+1
      
      id.m<-c(paste(this.m, collapse="==1 & "))
      id.m<-paste(id.m,"==1")
      # data for subgroup "x"
      df.sub<-subset(df,eval(parse(text=id.m)))
      df.x<-data.table(df.sub)
      # Random splitting
      flag.consistency<-rep(0,n.splits)
      split_failures<-0.0
      
      for(bb in 1:n.splits){

      if(split_failures < count_failure){  
      set.seed(8316951+100*bb)
      # Random splitting
      if(split_method=="Random"){  
        n.split<-round(nrow(df.x)/2,0)
        datax.id<-c(1:nrow(df.x))
        idx.split1<-sample(datax.id,size=n.split,replace=FALSE)
        df.x.split1<-df.x[idx.split1,]
        idx.split2<-setdiff(datax.id,idx.split1)
        df.x.split2<-df.x[idx.split2,]
      }
      if(split_method=="SPlit"){
      SPlitIndices<-quiet(SPlit(data=df.x[,c("Y","Event","Treat")], splitRatio=0.5,  tolerance = 1e-4, maxIterations=100, nThreads=1))
      df.x.split1 <- df.x[SPlitIndices, ]
      df.x.split2 <- df.x[-SPlitIndices, ]
      }
        hr.split1<-try(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split1,robust=FALSE))$conf.int,TRUE)
        hr.split2<-try(summary(coxph(Surv(Y,Event)~Treat,data=df.x.split2,robust=FALSE))$conf.int,TRUE)
        
        if(!inherits(hr.split1,"try-error") & !inherits(hr.split2,"try-error")){
          
          flag.consistency[bb]<-c(hr.split1[1,1]>hr.consistency & hr.split2[1,1]> hr.consistency)
  
          split_failures<-sum((1-flag.consistency)[1:bb],na.rm=TRUE)
          
          } # both split estimable
        }
          } # end splits 
      
      p.consistency<-mean(flag.consistency,na.rm=TRUE)
      
      # Trigger stopping for futility
      if(p.consistency<=pstop_futile) stop_futile<-TRUE
      
      if(details) cat("Consistency",p.consistency,"\n")
      
      if(p.consistency>pconsistency.threshold){
        any.found<-any.found+1
        k<-length(this.m)
        covsm<-rep("M",maxk) # store at most kmax names 
        mindex<-c(1:maxk)
        Mnames<-paste(covsm,mindex,sep=".")
        mfound<-matrix(rep("",maxk))
        mfound[c(1:k)]<-this.m
        mfound<-c(mfound)
        
        resultk<-c(p.consistency,found.hrs$n[m],found.hrs$grp[m],m,k,mfound)
        
        resultk<-transpose(as.data.table(resultk))
        
        names(resultk)<-c("p.consistency","Nsg","group.id","m.index","K",Mnames)
        result<-rbind(result,resultk)
 
        if(details) cat("Splitting method, # of splits=",c(split_method,n.splits),"\n")
        if(details)  cat("Model, % Consistency Met=",c(this.m,p.consistency),"\n")
        }
      }
    }
# Choice of subgroup  
if(any.found>0){
if(details) cat("Number of subgroups meeting consistency criteria=",c(any.found),"\n")

if(details) print(result)

result$p.consistency<-as.numeric(result$p.consistency)
result$Nsg<-as.numeric(result$Nsg)
result$K<-as.numeric(result$K)
  
result_new<-copy(result)

## Among SG's meeting (minimal) consistency criteria
# Choose sg with highest consistency with priority on smallest K   
if(sg_focus=="hr") setorder(result_new,-p.consistency,K)
if(sg_focus=="hrOnly") setorder(result_new,-p.consistency,-Nsg)
# Choose sg with largest size 
# Emphasis by Nsg, consistency, and smallest K
if(sg_focus=="Nsg") setorder(result_new,-Nsg,-p.consistency,K)
if(sg_focus=="NsgOnly") setorder(result_new,-Nsg,-p.consistency)
}
  
  this.1<-NULL
  data.found<-NULL
  grp1<-NULL
  if(!is.null(result_new)){
 # show (up to) top 10
 if(details) print(result_new[1:min(10,nrow(result_new)),])
    
    # Plot top 1
    # Subgroup corresponding to the top
    if(nrow(result_new)==1){
    top_result<-result_new  
    }
   
    if(nrow(result)>1){
    top_result<-result_new[1,]
    }
   
    grp1<-as.numeric(top_result$group.id)
    m1<-as.numeric(top_result$m.index)
   
    index1<-as.numeric(unlist(index.Z[m1,]))
    this.1<-names.Z[index1==1]
    
    id.harm<-c(paste(this.1, collapse="==1 & "))
    id.harm<-paste(id.harm,"==1")
    df.sub<-subset(df,eval(parse(text=id.harm)))
    df.sub$treat.recommend<-0
    
    id.noharm<-c(paste(this.1, collapse="!=1 | "))
    id.noharm<-paste(id.noharm,"!=1")
    df.subC<-subset(df,eval(parse(text=id.noharm)))
    df.subC$treat.recommend<-1
    
    data.found<-rbind(df.sub,df.subC)
    
    # Original analysis dataset is output with harm flag
    
    if(details & plot.sg){
      plot.subgroup(sub1=df.sub,sub1C=df.subC,tte.name="Y",event.name="Event",treat.name="Treat",subid=paste(this.1,collapse="/"),byrisk=by.risk)
      #print(head(df[m1,c(1:6)]))
      print(this.1)
      cat("% consistency criteria met=",c(top_result$p.consistency),"\n")
    }
  }
  out<-list(result=result_new,result_start=result,sg.harm=this.1,data=data.found,sg.harm.id=grp1)
  return(out)
}

