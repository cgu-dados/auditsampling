#' @title Outputs a table with SRS based inferences
#'
#' @description This function do Simple Random Sample inferences from a sample data,
#' given a confidence level. It can calculate estimates for means and totals and proportions too.
#' It have option to print a Word *.docx table with the estimate parameters.
#' @param sample.data A data.frame representing a sample with the collected data.
#' @param num.cols The number of the columns where the data is continuous.
#' @param cat.cols The number of the columns where the data is categorical.
#' @param alpha Confidence level.
#' @param type If you want means or totals.
#' @param print.report \code{FALSE} by default. If \code{TRUE}, prints a .docx table.
#' @param labels An optional string vector of labels to be putted into printed table.
#' Must have the same length of the sum of numerical columns and all levels of categorical columns.
#' @import flextable
#' @return A data.frame containing the fields:
#' \describe{
#' \item{parameter} The label of the parameter.
#' \item{n} The effective sample size.
#' \item{cv} Variation coefficient
#' \item{point.estimate} Point estimate.
#' \item{interval} Confidence interval.
#' }
#' @export
#' @import flextable
#' @example
#' data('iris)
#'
#' infereceSRS(iris,num.cols = 1:4,cat.cols = 5,alpha = 0.9, N=56892,
#' type = "total")
#'
#' rotulos <- c("Media da espessura da petala","Media do comprimento da petala",
#' "Media da espessura da sepala","Media do comprimento da sepala",
#' "Proporcao de setosas","Proporcao de versicolor", "Proporcao de virginica")
#' infereceSRS(iris,num.cols = 1:4,cat.cols = 5,alpha = 0.95, N=56892,
#'              type = "mean",labels = rotulos,print.report = TRUE)

infereceSRS <- function(sample.data,num.cols=NULL,cat.cols=NULL,alpha,N=Inf,
                    type=c("mean","total"),print.report=FALSE,labels=NULL){
  numeric.data=sample.data[,num.cols]
  formula=as.formula(paste0("~",names(sample.data)[cat.cols],"-1",collapse="+"))
  categoric.data=model.matrix(formula,sample.data)
  inferences=data.frame(matrix(nrow = 1,ncol=5))
  if(!is.null(labels) && length(labels)!=(ncol(numeric.data)+ncol(categoric.data))){
    stop("Your label vector do not have the same length of the parameters to be estimated.")
  }
  inference.data=cbind.data.frame(numeric.data,categoric.data)
  z_alpha=qnorm(1 - ((1 - alpha) / 2))
  for(i in 1:length(num.cols)){
    if(type=="mean"){
      p.estimate=mean(inference.data[,i],na.rm=TRUE)
      n=length(inference.data[!is.na(inference.data[,i]),i])
      sup.int=p.estimate+z_alpha*(sd(inference.data[,i],na.rm=TRUE)/sqrt(n))
      inf.int=p.estimate-z_alpha*(sd(inference.data[,i],na.rm=TRUE)/sqrt(n))
      if(!is.infinite(N)){
        sup.int=p.estimate+z_alpha*(sd(inference.data[,i],na.rm=TRUE)/sqrt(n))*
        (1 - (n/N))
        inf.int=p.estimate-z_alpha*(sd(inference.data[,i],na.rm=TRUE)/sqrt(n))*
          (1 - (n/N))
      }
      inferences[i,1]=paste0("M?dia de ",names(inference.data)[i])
      inferences[i,2]=n
      inferences[i,3]=round(sd(inference.data[,i],na.rm=TRUE)/
        mean(inference.data[,i],na.rm=TRUE),3)
      inferences[i,4]=round(p.estimate,3)
      inferences[i,5]=paste0("[",round(inf.int,3),";",
                             round(sup.int,3),"]")
    }
    if(type=="total"){
      p.estimate=N*mean(inference.data[,i],na.rm=TRUE)
      n=length(inference.data[!is.na(inference.data[,i]),i])
      sup.int=p.estimate+z_alpha*((N*sd(inference.data[,i],na.rm=TRUE))/sqrt(n))
      inf.int=p.estimate-z_alpha*(N*(sd(inference.data[,i],na.rm=TRUE))/sqrt(n))
      if(!is.infinite(N)){
        sup.int=p.estimate+z_alpha*((N*sd(inference.data[,i],na.rm=TRUE))/sqrt(n))*
          (1 - (n/N))
        inf.int=p.estimate-z_alpha*((N*sd(inference.data[,i],na.rm=TRUE))/sqrt(n))*
          (1 - (n/N))
      }
      inferences[i,1]=paste0("Total de ",names(inference.data)[i])
      inferences[i,2]=n
      inferences[i,3]=round(sd(inference.data[,i],na.rm=TRUE)/
        mean(inference.data[,i],na.rm=TRUE),3)
      inferences[i,4]=round(p.estimate,3)
      inferences[i,5]=paste0("[",round(inf.int,3),";",
                             round(sup.int,3),"]")
    }
  }
  for(i in (length(num.cols)+1):ncol(inference.data)){
    if(type=="mean"){
      p.estimate=mean(inference.data[,i],na.rm=TRUE)
      n=length(inference.data[!is.na(inference.data[,i]),i])
      sup.int=p.estimate+z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)
      inf.int=p.estimate-z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)
      if(!is.infinite(N)){
        sup.int=p.estimate+z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)*sqrt((N-n)/(N-1))
        inf.int=p.estimate-z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)*sqrt((N-n)/(N-1))
      }
      inferences[i,1]=paste0("Propor??o m?dia de ",names(inference.data)[i])
      inferences[i,2]=n
      inferences[i,3]=round(sqrt(p.estimate*(1 - p.estimate))/
        p.estimate,3)
      inferences[i,4]=round(p.estimate,3)
      inferences[i,5]=paste0("[",round(inf.int,3),";",
                             round(sup.int,3),"]")
    }
    if(type=="total"){
      p.estimate=mean(inference.data[,i],na.rm=TRUE)
      p.total=N*p.estimate
      n=length(inference.data[!is.na(inference.data[,i]),i])
      sup.int=p.estimate+z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)
      inf.int=p.estimate-z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)
      if(!is.infinite(N)){
        sup.int=p.estimate+z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)*sqrt((N-n)/(N-1))
        inf.int=p.estimate-z_alpha*sqrt((p.estimate*(1 - p.estimate))/n)*sqrt((N-n)/(N-1))
      }
      inferences[i,1]=paste0("Quantidade total de ",names(inference.data)[i])
      inferences[i,2]=n
      inferences[i,3]=round(sqrt(p.estimate*(1 - p.estimate))/
        p.estimate,3)
      inferences[i,4]=p.total
      inferences[i,5]=paste0("[",round(N*inf.int,0),";",
                             round(N*sup.int,0),"]")
    }

  }
  names(inferences)=c("parameter","n","cv","point.estimate","interval")
  if(!is.null(labels)){
    inferences[,1]=labels
  }
  if(print.report==FALSE){return(inferences)}
  if(print.report==TRUE){
    inf.table=flextable::flextable(inferences)
    inf.table=flextable::set_header_labels(inf.table,
                                           parameter="Par?metro",n="n",cv="CV",
                                           point.estimate="Estimativa pontual",
                                           interval="Intervalo")
    inf.table=flextable::width(inf.table,width = c(2.5,0.5,0.5,1,1.5))
    inf.table=flextable::align(inf.table,j=1,align = "left",part="all")
    inf.table=flextable::align(inf.table,j=2:3,align = "center",part="all")
    print(inf.table,preview="docx")
  }
}
