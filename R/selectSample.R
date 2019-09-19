#' @title Select sample and returns the \code{data.frame}
#'
#' @description Sample and outputs a data.frame with desired technique and size.
#' @param data The data.frame containing the population data.
#' @param stratum A vector with strata variable names.
#' @param method Select between \code{srswr}, \code{srswor} or \code{systematic}.
#' @param size The size of sample. Must provide a vector if the sample is stratified.
#' @param seed The seed for the random generator.
#' @param resample.data A data frame used as pilot sample. Its useful to ignore the observations sampled before and to maintain the original fpc.
#' @return A \code{data.frame} additionally with \code{Prob} and \code{fpc} columns. If the systematic
#' method of selection has been chosen, returns the \code{k} and cumulated \code{k} also.
#' @import sampling
#' @example
#' pilot <- selectSample(siasgdw[siasgdw$id_unidade == 170010, ],
#' method = "srswor", size = 30, seed = 150)

selectSample <- function(data,stratum=NULL,method,size,seed,resample.data=NULL){
  if(!is.null(resample.data) && !is.data.frame(resample.data)){
    stop("Resample data must be a data.frame.")
  }
  if(!is.null(resample.data) && length(intersect(tolower(names(resample.data)),
                                                 c("id_unit", "fpc"))) != 2){
    stop("Resample data must have id units and fpc columns.")
  }
  if(!is.null(resample.data)){
    fpc_column=which(tolower(names(resample.data)) == "fpc")
    if(!is.null(stratum) && any(tolower(names(resample.data)) == "stratum")){
      stratum_column=which(names(resample.data) == stratum)
      gfpc=unique(resample.data[,c(fpc_column,stratum_column)])
    }else{
      gfpc=unique(resample.data[,fpc_column])
    }
    data=data[-resample.data$ID_unit,]
  }
  if(!method %in% c("srswor","srswr","systematic")){
    stop("Method argument must be 'srswor','srswr' or 'systematic'.")
  }
  if(!is.null(stratum)){
    if(is.character(stratum)){stratum=which(names(data)==stratum)}
    if(method == "systematic"){
      sample=data.frame()
      if(length(stratum)>1){
        strata=do.call(paste, c(data[stratum], sep = "_"))
        data$stratum=strata
        strata=unique(strata)
      }else{
        strata=as.character(unique(data[,stratum]))
        data$stratum=data[,stratum]
      }
      for(i in seq_along(strata)){
        d=data[which(data$stratum==strata[i]),]
        d$ID_unit=1:nrow(d)
        d$Prob=size[i]/nrow(d)
        d$k=nrow(d)/size[i]
        set.seed(seed+i)
        s=d[as.logical(sampling::UPsystematic(pik=rep(size[i]/nrow(d),
                                                              nrow(d)))),]
        s$k_sum=cumsum(s$k)
        s=transform(s, fpc = round(sum(1/Prob),0))
        sample=rbind.data.frame(sample,s)
      }
      if(!is.null(resample.data)){
        if(!is.null(stratum)){
          sample$fpc=NULL
          id=as.numeric(c(rownames(resample.data), rownames(sample)))
          sample=merge(sample,gfpc,by.x="stratum", by.y = names(gfpc)[2])
          sample=smartbind(resample.data,sample)
          sample$ID_unit=id
          sample=sample[order(sample$ID_unit),]
          return(sample)
        }else{
          sample$fpc=gfpc
          id=as.numeric(c(rownames(resample.data), rownames(sample)))
          sample=smartbind(resample.data,sample)
          sample$ID_unit=id
          sample=sample[order(sample$ID_unit),]
          return(sample)
        }
      }else{
        return(sample)
      }

    }else{
    set.seed(seed)
    st=sampling::strata(data=data,stratanames = names(data)[stratum],
                 size=size,method=method)
    sample=sampling::getdata(data,st)
    nfpc=as.numeric(tapply(sample$Prob, sample$Stratum, function(x) sum(1/x)))
    sample$fpc=round(rep(nfpc, as.numeric(table(sample$Stratum))),0)
    if(!is.null(resample.data)){
      if(!is.null(stratum)){
        sample$fpc=NULL
        id=as.numeric(c(rownames(resample.data), rownames(sample)))
        sample=merge(sample,gfpc,by.x=names(sample)[stratum_column], by.y = names(gfpc)[2])
        sample=smartbind(resample.data,sample)
        sample$ID_unit=round(id,0)
        sample=sample[order(sample$ID_unit),]
        return(sample)
      }
    }else{
      return(sample)
    }
  }
  }else{
    if(method=="srswor"){
      data$ID_unit=1:nrow(data)
      data$Prob=size/nrow(data)
      set.seed(seed)
      sample=data[sample(1:nrow(data),size=size),]
      sample$fpc=sum(1/sample$Prob)
      sample=sample[order(sample$ID_unit),]
      if(!is.null(resample.data)){
        sample$fpc=gfpc
        id=as.numeric(c(rownames(resample.data), rownames(sample)))
        sample=smartbind(resample.data,sample)
        sample$ID_unit=id
        sample=sample[order(sample$ID_unit),]
        return(sample)
      }else{
        return(sample)
      }
    }
    if(method=="srswr"){
      data$ID_unit=1:nrow(data)
      data$Prob=size/nrow(data)
      set.seed(seed)
      sample=data[sample(1:nrow(data),size=size,replace=TRUE),]
      sample$fpc=sum(1/sample$Prob)
      sample=sample[order(sample$ID_unit),]
      if(!is.null(resample.data)){
        sample$fpc=gfpc
        id=as.numeric(c(rownames(resample.data), rownames(sample)))
        sample=smartbind(resample.data,sample)
        sample$ID_unit=round(id,0)
        sample=sample[order(sample$ID_unit),]
        return(sample)
      }else{
        return(sample)
      }
    }
    if(method=="systematic"){
      data$ID_unit=1:nrow(data)
      data$Prob=size/nrow(data)
      data$k=nrow(data)/size
      set.seed(seed)
      sample=data[as.logical(sampling::UPsystematic(pik=rep(size/nrow(data),
                             nrow(data)))),]
      sample$k_sum=cumsum(sample$k)
      sample$fpc=sum(1/sample$Prob)
      sample=sample[order(sample$ID_unit),]
      if(!is.null(resample.data)){
        sample$fpc=gfpc
        id=as.numeric(c(rownames(resample.data), rownames(sample)))
        sample=smartbind(resample.data,sample)
        sample$ID_unit=id
        sample=sample[order(sample$ID_unit),]
        return(sample)
      }else{
        return(sample)
      }
    }
  }
}

smartbind <- function(df1,df2){
  rbind.data.frame(
    data.frame(c(df1, sapply(setdiff(names(df2), names(df1)), function(x) NA))),
    data.frame(c(df2, sapply(setdiff(names(df1), names(df2)), function(x) NA)))
  )
}
