#' @title Returns the number of elements per strata based on Neyman optimum allocation
#'
#' @description Based on data from sample and population, this function returns the optimum allocation of elements.
#' Users can input strata vectors from sample and from population, but the last is not required.
#' The output is a vector.
#' @param n The size of sample.
#' @param strat.vector The vector of strata from population data.
#' @param sample.vector The vector of strata from sample data. Required.
#' @param sd.vector The data sample vector from which the standards deviations will be calculated.
#' @return Vector with the same length as the unique labels of strata containing the size of each strata for sampling.
#' @example
#' library(survey)
#' data('api')
#'
#' ### CONTINUOUS 
#' neymanAllocation(n = 65, strat.vector = apipop$stype, sample.vector = apistrat$stype,
#'                  sd.vector = apistrat$enroll)
#' ### CATEGORICAL
#' neymanAllocation(n = 65, strat.vector = apipop$stype, sample.vector = apistrat$stype,
#'                  sd.vector = apistrat$sch.wide)

neymanAllocation <- function(n,strat.vector=NULL,sample.vector,sd.vector){
  if(length(sample.vector) != length(sd.vector)){
    stop("The 'sample.vector' must have the same length as the 'sd.vector'")
  }
  if(!is.null(strat.vector)){
    Wh <- as.numeric(table(strat.vector)) / sum(table(strat.vector))
  }else{
    Wh <- as.numeric(table(sample.vector)) / sum(table(sample.vector))
  }
  if(is.numeric(sd.vector)){
    sd <- tapply(sd.vector,sample.vector,sd)
  }
  if(is.factor(sd.vector) | is.character(sd.vector)){
    sd <- tapply(sd.vector,sample.vector,
                 function(x) prop.table(table(x)) * (1-prop.table(table(x))))
    sd <- as.numeric(levels(factor(unlist(sd))))
  }
  nh <- n * (Wh * sd / sum(Wh * sd))
  return(round(nh,0))
}
