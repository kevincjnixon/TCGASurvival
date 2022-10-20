library(survival)
library(survminer)


#survival<-read.delim("rawData/Survival_SupplementalTable_S1_20171025_xena_sp")
#rownames(survival)<-gsub("-",".",survival$sample)
#survival$group<-ifelse(survival$vital_status=="Alive",1,2)
# sample<-read.delim("~/../Downloads/TCGA_phenotype_denseDataOnlyDownload.tsv.gz")
# rownames(sample)<-gsub("-",".", sample$sample)
# sample<-sample[match(rownames(survival), rownames(sample)),]
# all.equal(rownames(sample), rownames(survival))
# survival<-cbind(sample, survival)

#exp<-read.delim("rawData/geneLevel.txt.gz")
#rownames(exp)<-exp$X
#exp<-exp[,-1]

#survival<-survival[which(rownames(survival) %in% colnames(exp)),]
#exp<-exp[,which(colnames(exp) %in% rownames(survival))]

#saveRDS(exp, "glTPM.RDS")
#saveRDS(survival, "survival.RDS")



message("Setting up functions...")
subDat<-function(exp, survival, type=NULL, gender=NULL, race=NULL, rmNorm=T, prim=F){
  if(!is.null(type)){
    survival<-survival[which(survival$cancer.type.abbreviation %in% type),]
  }
  if(!is.null(gender)){
    survival<-survival[which(survival$gender %in% gender),]
  }
  if(!is.null(race)){
    survival<-survival[which(survival$race %in% race),]
  }
  if(isTRUE(rmNorm)){
    survival<-survival[which(!survival$sample_type_id==11),]
  }
  if(isTRUE(prim)){
    survival<-survival[which(survival$sample_type_id==1),]
  }
  #print(dim(survival))
  return(list(exp=exp[,which(colnames(exp) %in% rownames(survival))],
              survival=survival))
}

rowMedians<-function(x){
  return(apply(x, 1, FUN=function(x) as.numeric(median(x, na.rm=T))))
}

plotSurv<-function(stratum, time, status, title="", timeLim=NULL, retPlot=F, show.conf=T){
  require(survminer)
  dat<-data.frame(stratum=stratum, time=time, status=status)
  if(!is.null(timeLim)){
    dat<-subset(dat, time<=timeLim)
  }

  fit<-surv_fit(Surv(time,status)~stratum, data=dat)
  #fit<-coxph.fit(Surv(time,status)~stratum, data=dat)
  #print(dim(dat))
  g<-ggsurvplot(
    fit=fit,
    risk.table=T,
    xlab = "Days",
    ylab="Probability",
    pval=T,
    conf.int=show.conf,
    title=title
  )
  if(isTRUE(retPlot)){
    return(g)
  } else {
    print(g)
  }
  #return(g)
}

ExpStrat<-function(x, y=NULL, method="median", useRatio=F, retStrat=T, pc=0, retDat=F){
  colMedian<-function(x){
    return(apply(x, 2, FUN=function(x) as.numeric(median(x, na.rm=T))))
  }
  res.dat<-data.frame(row.names=colnames(x), expression=as.numeric(x[1,]))
  x.name<-rownames(x[1,])
  threshold<-c()
  if(nrow(x)<2){
    x<-as.numeric(x)
    threshold<-median(x, na.rm=T)
    if(method=="mean"){
      threshold<-mean(x, na.rm=T)
    }
    if(method=="geoMean"){
      threshold<-BinfTools:::geoMean(x)
    }
    #return(ifelse(x>threshold,"High","Low"))
  } else {
    threshold<-NULL
    if(method=="median"){
      threshold<-median(unlist(x), na.rm=T)
      x<-colMedian(as.matrix(x))
    }
    if(method=="mean"){
      threshold<-mean(unlist(x), na.rm=T)
      x<-colMeans(x, na.rm=T)
    }
    # if(method=="geoMean"){
    #   threshold<-BinfTools:::geoMean(unlist(x))
    #   x<-BinfTools:::rowGeoMean(x)
    # }
  }
  res.dat<-cbind(res.dat, ifelse(x>threshold, "high","low"))
  colnames(res.dat)<-c(paste0(x.name,".expression"),paste0(x.name,".category"))
  if(is.null(y)){
    if(isTRUE(retDat)){
      return(res.dat)
    }
    if(isTRUE(retStrat)){
      return(ifelse(x>threshold, "High","Low"))
    } else {
      return(x)
    }
  } else {
    x.ratio<-x
    if(isFALSE(useRatio)){
      x<-ifelse(x>threshold, "x.High","x.Low")
    }
    res.dat<-cbind(res.dat, as.numeric(y[1,]))
    y.name<-rownames(y[1,])
    if(nrow(y)<2){
      #message("convert to numeric")
      y<-as.numeric(y)
      #message("get median")
      threshold<-median(y, na.rm=T)
      if(method=="mean"){
        threshold<-mean(y, na.rm=T)
      }
      if(method=="geoMean"){
        threshold<-BinfTools:::geoMean(y)
      }
      #return(ifelse(y>threshold,"High","Low"))
    } else {
      threshold<-NULL
      if(method=="median"){
        threshold<-median(unlist(y), na.rm=T)
        y<-colMedian(as.matrix(y))
      }
      if(method=="mean"){
        threshold<-mean(unlist(y), na.rm=T)
        y<-colMeans(y, na.rm=T)
      }
      # if(method=="geoMean"){
      #   threshold<-BinfTools:::geoMean(unlist(y))
      #   y<-BinfTools:::rowGeoMean(y)
      # }
    }
    res.dat<-cbind(res.dat, ifelse(y>threshold,"high","low"))
    colnames(res.dat)<-c(paste0(x.name,".expression"),paste0(x.name,".category"),
                         paste0(y.name,".expression"),paste0(y.name,".category"))
    if(isTRUE(retDat)){
      return(res.dat)
    }
    if(isFALSE(useRatio)){
      if(isTRUE(retStrat)){
        #message("set threshold")
        y<-ifelse(y>threshold, "y.High","y.Low")
        return(paste(x,y,sep="_"))
      } else {
        return(data.frame(x=x, y=y))
      }
    } else {
      #message("find ratio")
      ratio<-log2((pc+x.ratio)/(pc+y))
      med.rat<-median(ratio)
      #print(head(ratio))
      #print(med.rat)
      if(isTRUE(retStrat)){
        return(ifelse(ratio>med.rat, "High","Low"))
      } else {
        return(ratio)
      }
    }
  }
}

optStrat<-function(x,y=NULL, data, time, event, type=c("rl","cph"),
                   title="", useRatio=F, retPlot=F, timeLim=NULL, retDat=F, show.conf=T){
  dat<-data.frame(cbind(rowMedians(t(data[which(rownames(data) %in% x),])),time, event))
  if(!is.null(timeLim)){
    dat<-subset(dat, time<=timeLim)
  }
  #print(head(dat))
  res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
  res.cat<-surv_categorize(res.cut)
  res.dat<-data.frame(row.names=rownames(res.cat), expression=dat$V1, category=res.cat$V1)
  if(length(x)<2){
    colnames(res.dat)<-c(paste0(x,".expression"), paste0(x,".category"))
  }
  #print(head(res.cat))
  if(!is.null(y)){
    x.dat<-dat
    dat<-data.frame(cbind(rowMedians(t(data[which(rownames(data) %in% y),])), time, event))
    if(!is.null(timeLim)){
      dat<-subset(dat, time<=timeLim)
    }
    res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
    cat2<-surv_categorize(res.cut)
    res.dat<-cbind(res.dat, dat$V1, cat2$V1)
    if(length(y)<2){
      colnames(res.dat)<-c(paste0(x,".expression"), paste0(x,".category"), paste0(y,".expression"), paste0(y,".category"))
    }
    if(isTRUE(useRatio)){
      #print(head(x.dat))
      #print(head(dat))
      dat$V1<-log2(x.dat$V1/dat$V1)
    }
    res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
    cat2<-surv_categorize(res.cut)
    if(isFALSE(useRatio)){
      #print(head(cat2))
      res.cat<-data.frame(time=res.cat$time, event=res.cat$event, V1=paste(res.cat$V1, cat2$V1, sep="."))
      #print(head(res.cat))
    }
  }

  fit<-surv_fit(Surv(time, event)~V1, data=res.cat)
  p<-TRUE
  if(type[1]=="cph"){
    #print(lengths(split(res.cat$V1, res.cat$V1)))
    res.cox<-coxph(Surv(time, event)~V1, data=res.cat)
    p<-round(summary(res.cox)$logtest[3],digits=3)
    new_df<-with(res.cat,
                data.frame(V1=c("high","low")))
    #print(res.cox)
    fit<-survfit(res.cox, newdata=new_df)
  }

  #print(class(fit))

  g<-ggsurvplot(
    fit=fit,
    data=res.cat,
    risk.table = T,
    xlab = "Days",
    ylab="Probability",
    pval=p,
    conf.int=show.conf,
    title=title
  )
  if(isTRUE(retDat)){
    return(list(data=res.dat, fit=fit))
  }
  if(isTRUE(retPlot)){
    return(g)
  } else {
    print(g)
  }
}


# dat<-optStrat("PCSK4","FAM153A",exp, survival$OS.time, survival$OS, "", retDat=T)
# dat<-ExpStrat(exp[which(rownames(exp) %in% "PCSK4"),], exp[which(rownames(exp) %in% "FAM153A"),],
#               retDat=T)

plotStrat<-function(dat){
  require(dplyr)
  toPlot<-dat[,1,drop=F] %>% tidyr::gather(key="Gene", value="Expression") %>% dplyr::mutate(category=dat[,2])
  if(ncol(dat)>2){
    tmp<-dat[,3,drop=F] %>% tidyr::gather(key="Gene", value="Expression") %>% dplyr::mutate(category=dat[,4])
    toPlot<-rbind(toPlot, tmp)
  }
  toPlot$Gene<-sapply(strsplit(toPlot$Gene,".",T),'[[',1)
  toPlot$Expression<-log2(0.001+toPlot$Expression)
  g<-ggplot2::ggplot(toPlot, ggplot2::aes(x=Gene, y=Expression, fill=category))+
    ggplot2::geom_violin(position=ggplot2::position_dodge(width=0.9))+theme_minimal()+
    ggplot2::geom_point(position=ggplot2::position_jitterdodge(seed=1, dodge.width = 0.9), alpha=0.01) +
    ggplot2::labs(y="Log2(0.001 + TPM)", x="")
  print(g)
}

splitConf<-function(x){
  res<-list()
  for(i in 1:ncol(x)){
    res[[i]]<-lengths(split(rownames(x), x[,i]))
    names(res)[i]<-colnames(x)[i]
  }
  return(res)
}

coxForest<-function(x,y=NULL, data, time, event, confounders=NULL, minConf=2,
                   title="", useRatio=F, retPlot=F, timeLim=NULL, retDat=F, refLow=F){
  dat<-data.frame(cbind(rowMedians(t(data[which(rownames(data) %in% x),])),time, event))
  if(!is.null(timeLim)){
    dat<-subset(dat, time<=timeLim)
  }
  #print(head(dat))
  res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
  res.cat<-surv_categorize(res.cut)
  res.dat<-data.frame(row.names=rownames(res.cat), expression=dat$V1, category=res.cat$V1)
  if(length(x)<2){
    colnames(res.dat)<-c(paste0(x,".expression"), paste0(x,".category"))
  }
  #print(head(res.cat))
  if(!is.null(y)){
    x.dat<-dat
    dat<-data.frame(cbind(rowMedians(t(data[which(rownames(data) %in% y),])), time, event))
    if(!is.null(timeLim)){
      dat<-subset(dat, time<=timeLim)
    }
    res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
    cat2<-surv_categorize(res.cut)
    res.dat<-cbind(res.dat, dat$V1, cat2$V1)
    if(length(y)<2){
      colnames(res.dat)<-c(paste0(x,".expression"), paste0(x,".category"), paste0(y,".expression"), paste0(y,".category"))
    }
    if(isTRUE(useRatio)){
      #print(head(x.dat))
      #print(head(dat))
      dat$V1<-log2(x.dat$V1/dat$V1)
    }
    res.cut<-survminer::surv_cutpoint(dat, "time", "event", variables="V1")
    cat2<-surv_categorize(res.cut)
    if(isFALSE(useRatio)){
      #print(head(cat2))
      res.cat<-data.frame(time=res.cat$time, event=res.cat$event, V1=paste(res.cat$V1, cat2$V1, sep="."))
      #print(head(res.cat))
    }
  }
  #Add in confounders - must be in the form of a data frame (same row order as data), must be factor object with levels specified
  conf<-NULL
  if(!is.null(confounders)){
    #Check to see if any confounders have fewer than minConf entries for a bin- if so, mask them
    if(any(lapply(splitConf(confounders), min)<=minConf)){
      lens<-splitConf(confounders)
      for(i in 1:length(lens)){
        if(class(confounders[,i])=="factor"){
          if(any(lens[[i]])<=minConf){
            confounders[which(confounders[,i] %in% unique(names(lens[[i]][which(lens[[i]]<=minConf)]))),]<-rep(NA, ncol(confounders))
          }
        }
      }
    }
    print(head(confounders))
    #These are factor confounders (like gender, race, etc)
    res.cat<-cbind(res.cat, confounders)
    print(lapply(res.cat, class))
    conf<-paste(colnames(confounders), collapse=" + ")
  }
  if(is.null(y)){
    if(isTRUE(refLow)){
      res.cat$V1<-factor(res.cat$V1, levels=c("low","high"))
    }
    print(head(res.cat))
    colnames(res.cat)[which(colnames(res.cat) %in% "V1")]<-x
  }
  if(!is.null(y)){
    tmp<-paste(x,y,sep=".")
    if(isTRUE(useRatio)){
      tmp<-paste(x,y,sep="/")
    }
    colnames(res.cat)[which(colnames(res.cat) %in% "V1")]<-tmp
    x<-tmp
  }
  #print(head(res.cat))
  form<-as.formula(paste("Surv(time,event)~",x))
  res.cox<-coxph(form, data=res.cat)
  if(!is.null(conf)){
    tmp<-paste("Surv(time,event)~",x,"+")
    form<-as.formula(paste(tmp, conf))
    print(form)
    res.cox<-coxph(form, data=res.cat)
  }
  g<-ggforest(res.cox, data=res.cat, main=title)
  if(isTRUE(retDat)){
    return(res.dat)
  }
  if(isTRUE(retPlot)){
    return(g)
  } else {
    print(g)
  }
}


message("Ready!")


#optStrat(x="PCSK4", data=subDat(exp, survival, type="HNSC")$exp,
#         time=subDat(exp, survival, type="HNSC")$survival$death_days_to,
#         event=subDat(exp, survival, type="HNSC")$survival$group, title="PCSK4")

#plotSurv(ExpStrat(x=subDat(exp, survival, type="HNSC")$exp[which(rownames(exp) %in% "PCSK4"),],
#                  method="median"), time=subDat(exp, survival, type="HNSC")$survival$death_days_to,
#         status=subDat(exp, survival, type="HNSC")$survival$group, title="PCSK4")
