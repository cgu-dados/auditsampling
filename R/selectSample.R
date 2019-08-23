#' @title Select sample and returns the \code{data.frame}
#'
#' @description Sample and outputs a data.frame with desired technique and size.
#' @param data The data.frame containing the population data.
#' @param stratanames A vector with strata variable names.
#' @param method Select between \code{srswr}, \code{srswor} or \code{systematic}.
#' @param size The size of sample. Must provide a vector if the sample is stratified.
#' @param seed The seed for the random generator.
#' @return A \code{data.frame} additionally with \code{Prob} and \code{fpc} columns. If the systematic
#' method of selection has been chosen, returns the \code{k} and cumulated \code{k} also.
#' @import sampling
#' @example
#' pilot <- selectSample(siasgdw[siasgdw$id_unidade == 170010, ], method = "srswor", 
#' size = 30, seed = 150)

selectSample <- function(data,stratum=NULL,method,size,seed){
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
        strata=unique(data[,stratum])
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
      return(sample)
    }else{
    set.seed(seed)
    st=sampling::strata(data=data,stratanames = names(data)[stratum],
                 size=size,method=method)
    sample=sampling::getdata(data,st)
    nfpc=as.numeric(tapply(sample$Prob, sample$Stratum, function(x) sum(1/x)))
    sample$fpc=round(rep(nfpc, as.numeric(table(sample$Stratum))),0)
    return(sample)}
  }else{
    if(method=="srswor"){
      data$ID_unit=1:nrow(data)
      data$Prob=size/nrow(data)
      set.seed(seed)
      sample=data[sample(1:nrow(data),size=size),]
      sample$fpc=sum(1/sample$Prob)
      return(sample[order(sample$ID_unit),])
    }
    if(method=="srswr"){
      data$ID_unit=1:nrow(data)
      data$Prob=size/nrow(data)
      set.seed(seed)
      sample=data[sample(1:nrow(data),size=size,replace=TRUE),]
      sample$fpc=sum(1/sample$Prob)
      return(sample[order(sample$ID_unit),])
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
      return(sample[order(sample$ID_unit),])
    }
  }
}

