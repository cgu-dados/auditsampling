#' @title Outputs a table with SRS and Stratified based inferences
#'
#' @description This function do Simple Random Sample and Stratified inferences from a sample data,
#' given an alpha level. It can calculate estimates for means, totals, proportions and ratios too.
#' It have option to print a Word *.docx table with the estimate parameters and two plotting graph types.
#' @param sample.data A data.frame representing a sample with the collected data.
#' @param num.y The index number of the columns (vector) where the data is continuous.
#' @param denom.y The index number of the columns (vector) that will be used as denominator for ratio estimators. Use NULL if you don't want to use ratio estimators or combine regular estimators with ratio ones.
#' @param total.x If you are using ratio estimators and want to estimate totals, you gonna need a total of \eqn{\bar{X}}. Use vectors of the same length as your estimated variables for estimate totals. If you don't provide totals, they will be estimated from data.
#' @param cat.y The index number of the columns where the data is categorical.
#' @param strat The index number of the columns to be used as stratified design.
#' @param post.strat The index number of the columns to used as a post stratified design.
#' @param alpha 1 - Confidence level.
#' @param N Population size.
#' @param type A vector of characters. Allowed means, totals or ratios. The function use \code{mean} as default.
#' @param fpc The index number of the column to be used as the factor of population correction.
#' @param labels An optional string vector of labels to be putted into printed table.
#' @param print.report \code{FALSE} by default. If \code{TRUE}, prints a .docx table.
#' Must have the same length of the sum of numerical columns and all levels of categorical columns.
#' @import flextable
#' @return A data.frame containing the following:
#' \describe{
#' \item{parameter}{The label of the parameter.}
#' \item{n}{The effective sample size.}
#' \item{se}{Standard error.}
#' \item{point.estimate}{Point estimate.}
#' \item{interval}{Confidence interval.}
#' }
#' @export
#' @import flextable
#' @import survey
#' @import extrafont
#' @example
#' library(survey)
#' data('api')
#'
#' inference(sample.data = apisrs, num.y = c(24,25),
#'           alpha = 0.05)
#' inference(sample.data = apisrs, num.y = c(36,12,13,13,25),
#'           cat.y = 16:19, alpha = 0.05, denom.x = c(NA,13,NA,12,NA,NA,NA,NA,NA))
#' inference(sample.data = apistrat, strat = 2, num.y = c(36,12,13,13,25),
#'           cat.y = 16:19, alpha = 0.05, denom.x = c(NA,13,NA,12,NA,NA,NA,NA,NA))
#' inference(sample.data = apistrat, post.strat = 2, num.y = c(36,12,13,13,25), N=nrow(apipop),
#'           cat.y = 16:19, alpha = 0.05, denom.x = c(NA,13,NA,12,NA,NA,NA,NA,NA))

inference <- function(sample.data,num.y=NULL,denom.x=NULL,total.x=NULL,
                      cat.y=NULL,strat=NULL,post.strat=NULL, alpha,N=Inf,
                    type=rep("mean",length(num.y)+length(cat.y)),fpc=NULL,labels=NULL,print.report=FALSE){
  ### Checks
  if(!is.null(post.strat) && is.infinite(N)){
    stop("If you input 'post.strata' you must provide 'N'.")
  }
  if(length(type) != length(num.y) + length(cat.y)){
    stop("The 'type' vector needs to be the same length as the sum of 'num.y' and 'cat.y'.")
  }
  if(is.null(denom.x)){
    denom.x <- rep(NA, length(num.y) + length(cat.y))
  }
  if(length(denom.x) != length(num.y) + length(cat.y)){
    stop("The 'denom.x' vector needs to be the same length as the sum of 'num.y' and 'cat.y'.")
  }
  if(is.null(total.x)){
    total.x <- rep(NA,length(num.y) + length(cat.y))
  }
  if(length(total.x) != length(num.y) + length(cat.y)){
    stop("The 'total.x' vector needs to be the same length as the sum of 'num.y' and 'cat.y'.")
  }
  if(alpha < 0  && alpha > 1){
    stop("The 'alpha' argument must to be between 0 and 1.")
  }
  if(any(num.y %in% cat.y)){
    stop("The vectors 'num.y' and 'cat.y' must be unique and contain no intersections.")
  }
  if(is.null(fpc)){
    fpc <- which(tolower(names(sample.data)) == "fpc")
    if(is.na(fpc)){
      fpc <- NULL
    }
  }
  ### Object svydesing creation
  FPC <- as.matrix(sample.data[,fpc])
  z <- qnorm(1 - (alpha/2))
  if(is.null(strat)){
    if(is.infinite(N)){
      options(warn = -1)
      svy <- svydesign(ids=~1,data=sample.data)
      options(warn = 1)
      N <- NULL
    }else{
      svy <- svydesign(ids=~1,fpc=FPC,data=sample.data)
    }}else{
    svy <- svydesign(ids=~1,strata=sample.data[,strat],
                     data=sample.data,fpc=FPC)
    }
  if(is.null(N) && !is.null(fpc) && !is.null(strat)){
    fp <- cbind.data.frame(svy$strata,svy$allprob)
    N <- sum(tapply(fp$probs, fp$strata, function(x) sum(1/x)))
  }
  if(is.null(N) && !is.null(fpc) && !is.null(post.strat)){
    fp <- cbind.data.frame(svy$strata,svy$allprob)
    N <- sum(tapply(fp$probs, fp$strata, function(x) sum(1/x)))
  }
  if(!is.null(post.strat)){ # Post stratification design
    var.name <- names(sample.data)[post.strat]
    aggr <- aggregate(as.formula(paste0("1:nrow(sample.data)",
                      " ~ ", paste0(var.name, collapse = " + "))),
                      data = sample.data, length)
    names(aggr)[ncol(aggr)] <- "Freq"
    aggr$Freq <- round(aggr$Freq / sum(aggr$Freq) * N, 0)
    svy <- postStratify(svy, as.formula(paste0("~",paste0(var.name, collapse = " + "))),aggr)
  }
  ### Inferences calculation
  inferences <- data.frame()
  j <- 1
  for(i in num.y){
    if(type[j] == "mean"){
      if(!is.na(denom.x[j])){
        if(is.na(total.x[j])){
          svy.cal <- calibrate(svy, formula = as.formula(paste0("~",names(sample.data)[denom.x[j]]," - 1")),
                               population = mean(sample.data[,denom.x[j]], na.rm = TRUE) *
                                 sum(unique(sample.data[,fpc])))
        }else{
          svy.cal <- calibrate(svy, formula = as.formula(paste0("~",names(sample.data)[denom.x[j]]," - 1")),
                               population = total.x[j])
        }
        infs <- svymean(sample.data[,i],svy.cal,na.rm=TRUE)
        point <- abs(as.numeric(infs))
        se <- as.numeric(SE(infs))
        interval <- c(point - (z * se), point + (z * se))
        n <- length(sample.data[,i][!is.na(sample.data[,i])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("M?dia de ",names(sample.data)[i]),
                                                  n=n, SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[1],3),
                                                  upr_int=round(interval[2],3)))
      }else{
        infs <- svymean(sample.data[,i],svy,na.rm=TRUE)
        point <- as.numeric(infs)
        se <- as.numeric(SE(infs))
        interval <- c(point - (z * se), point + (z * se))
        n <- length(sample.data[,i][!is.na(sample.data[,i])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("M?dia de ",names(sample.data)[i]),
                                                  n=n, SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[1],3),
                                                  upr_int=round(interval[2],3)))
      }
    }
    if(type[j] == "total"){
      if(!is.na(denom.x[j])){
        infs <- svyratio(sample.data[,i],sample.data[,denom.x[j]],svy,na.rm=TRUE)
        if(is.na(total.x[j])){
          infs <- predict(infs, total = mean(sample.data[,denom.x[j]], na.rm = TRUE) *
                            sum(unique(sample.data[,fpc])))
        }else{
          infs <- predict(infs, total = total.x[j])
        }
        point <- as.numeric(infs$total)
        se <- as.numeric(infs$se)
        interval <- c(point - (z * se), point + (z * se))
        n <- length(sample.data[,i][!is.na(sample.data[,i])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("Total de ",names(sample.data)[i]),
                                                  n=n, SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[1],3),
                                                  upr_int=round(interval[2],3)))
      }else{
        infs <- svytotal(sample.data[,i],svy,na.rm=TRUE)
        point <- as.numeric(infs)
        se <- as.numeric(SE(infs))
        interval <- c(point - (z * se), point + (z * se))
        n <- length(sample.data[,i][!is.na(sample.data[,i])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("Total de ",names(sample.data)[i]),
                                                  n=n, SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[1],3),
                                                  upr_int=round(interval[2],3)))
      }

    }
    if(type[j] == "ratio"){
      infs <- svyratio(sample.data[,i],sample.data[,denom.x[j]],svy,na.rm=TRUE)
      point <- as.numeric(infs$ratio)
      se <- as.numeric(SE(infs))
      interval <- c(point - (z * se), point + (z * se))
      n <- length(sample.data[,i][!is.na(sample.data[,i])])
      inferences <- rbind.data.frame(inferences,
                                     data.frame(var=paste0("Raz?o entre ",names(sample.data)[i],
                                                           " e ",names(sample.data)[denom.x[j]]),
                                                n=n, SE=round(se,4), p_estimate=round(point,3),
                                                lwr_int=round(interval[1],3),
                                                upr_int=round(interval[2],3)))
    }
    j <- j + 1
  }
  for(i in seq_along(cat.y)){
    formulas <- paste0("~I(", names(sample.data)[cat.y[i]],")=='",
    as.character(unique(sample.data[,cat.y[i]])) ,"'")
    if(type[j] == "mean"){
      for(k in seq_along(formulas)){
        infs <- svyciprop(as.formula(formulas[k]), svy,
                          method = "logit", level = 1 - alpha)
        point <- as.numeric(infs)
        se <- as.numeric(SE(infs))
        interval <- as.numeric(confint(infs))
        n <- length(sample.data[,cat.y[i]][!is.na(sample.data[,cat.y[i]])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("Propor??o de ",names(sample.data)[cat.y[i]],
                                                             " | ", unique(sample.data[,cat.y[i]])[k]),
                                                  n=n,
                                                  SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[1],3),
                                                  upr_int=round(interval[2],3)))
      }
    }
    if(type[j] == "total"){
      for(k in seq_along(formulas)){
        infs <- svytotal(as.formula(formulas[k]), svy)
        point <- round(as.numeric(infs[2]),0)
        se <- as.numeric(SE(infs)[2])
        interval <- as.numeric(confint(infs))
        n <- length(sample.data[,cat.y[i]][!is.na(sample.data[,cat.y[i]])])
        inferences <- rbind.data.frame(inferences,
                                       data.frame(var=paste0("Total de ",names(sample.data)[cat.y[i]],
                                                             " | ", unique(sample.data[,cat.y[i]])[k]),
                                                  n=n,
                                                  SE=round(se,4), p_estimate=round(point,3),
                                                  lwr_int=round(interval[2],3),
                                                  upr_int=round(interval[4],3)))
      }
    }
    j <- j + 1
  }
  if(!is.null(labels) && length(labels) == nrow(inferences)){
    inferences$var <- labels
  }
  return(inferences)
}


print_report <- function(data, alpha){
  names(data)=c("parameter","n","se","point.estimate","interval_lwr","interval_upr")
  inf.table=flextable::flextable(data)
  inf.table=flextable::font(inf.table, i = NULL, j = NULL, "Humanst521 BT", part = "all")
  inf.table=flextable::set_header_labels(inf.table,
                                         parameter="Par?metro",n="n",se="Erro Padr.",
                                         point.estimate="Estimativa pontual",
                                         interval_lwr=paste0("Intervalos*"),
                                         interval_upr="")
  inf.table=flextable::autofit(inf.table, add_w = 0)
  # inf.table=flextable::width(inf.table,width = c(2.5,0.5,0.5,0.5,1,1))
  inf.table=flextable::merge_at(inf.table, i = 1, j = 5:6, part = "header")
  inf.table=flextable::align(inf.table,j=1,align = "left",part="all")
  inf.table=flextable::align(inf.table,j=2:6,align = "right",part="all")
  inf.table=flextable::align(inf.table,j=5,align = "center",part="header")
  inf.table=flextable::add_footer_lines(inf.table, values = c(paste0("\n* N?vel de confian?a de ",
                       (1 - alpha) * 100,"%")))
  inf.table=flextable::font(inf.table, i = NULL, j = NULL, "Humanst521 BT", part = "footer")
  print(inf.table,preview="docx")
  inf.table=flextable::fontsize(inf.table, size = 10, part = "footer")
}

# print_report(table, 0.05)

install_cgu_themes <- function(){
  if(!"Humanst521 BT" %in% extrafont::fonts()){
    dir.create(paste0(getwd(),"/fonts"))
    download.file("https://github.com/cgu-dados/auditsampling/raw/master/fonts/humanst521-bt.zip",
                  destfile = "./fonts/humanst521-bt.zip", mode="wb",
                  quiet = TRUE)
    unzip("./fonts/humanst521-bt.zip", exdir = "./fonts")
    extrafont::font_import("./fonts", prompt = FALSE)
    extrafont::loadfonts(device="win",quiet = TRUE)
    unlink("./fonts", recursive = TRUE)
  }
}
