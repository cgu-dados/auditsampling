#' @title Calculates the statistics for hypotesis tests
#'
#' @description This function calculates de statistics for hypotesys tests considering the following hypotesis:
#' \deqn{H_{0}:\text{ } \mu \neq \mu_{0}}
#' \deqn{H_{0}:\text{ } \bar{p} \neq p_{0}}
#' \deqn{H_{0}:\text{ } \mu \leq \text{ or } \geq \mu_{0}}
#' \deqn{H_{0}:\text{ } \bar{p} \leq \text{ or } \geq p_{0}}
#' @param sample.data The source \code{data.frame} from which the funtion will search the columns.
#' @param num.col Vector of the numeric column numbers.
#' @param prop.col Vector of the binary column numbers.
#' @param mu0 Hypothesized mean.
#' @param p0 Hypothesized proportion.
#' @param alpha Level of significance.
#' @param alternative \code{1} if is a greater/lesser test and \code{2} for equal/different one.
#' @param labels Label vector for the outputs.
#' @param plot.graph Plot a confidence intervals for proportions only. Not implemented yet.
#' @return \code{data.frame} containing the following:
#' \describe{
#' \item{var}{The label of the variable/parameter.}
#' \item{point}{Point estimate.}
#' \item{tolerance}{The hypothesized mean/proportion.}
#' \item{sdd_err}{Standard error.}
#' \item{z_value}{Statistics z.}
#' \item{pval}{P-value.}
#' \item{power}{Test power, defined as (\eqn{1 - \beta}).}
#' }
#' @export
#' @examples
#' n=168
#' set.seed(101)
#' dados <- data.frame(CONTROLE1=rbinom(n,1,prob = 0.4),
#'                     CONTROLE2=rbinom(n,1,prob = 0.02),
#'                     CONTROLE3=rbinom(n,1,prob = 0.1),
#'                     CONTROLE4=rbinom(n,1,prob = 0.01),
#'                     CONTROLE5=rbinom(n,1,prob = 0.2),
#'                     CONTROLE6=rnorm(n,10,3),
#'                     CONTROLE7=rnorm(n,250,20))
#' test(sample.data = dados, num.col = 6:7, prop.col = 1:5, alpha = 0.05,
#'      mu0 = c(10,100),p0=c(0.05,0.05,0.10,0.10,0.1),
#'      alternative = c(rep(1,5),rep(2,2)))

test <- function(sample.data, num.col=NULL,prop.col=NULL,mu0=NULL,p0=NULL,
                 alpha, alternative=NULL,labels=NULL, plot.graph=FALSE){
  ### Checks
  options(scipen = 999)
  if(!is.data.frame(sample.data)){
    stop("'sample.data' must be a data.frame.")
  }
  if(is.null(num.col) && is.null(prop.col)){
    stop("You need to input 'num.col' and/or 'prop.col'.")
  }
  if(is.null(mu0) && is.null(p0)){
    stop("You need to input 'mu0' and/or 'p0'.")
  }
  if(!is.null(num.col) && is.null(mu0) | !is.null(prop.col) && is.null(p0)){
    stop("You need to input 'mu0' and/or 'p0'.")
  }
  if((length(mu0) + length(p0)) != (length(num.col) + length(prop.col))){
    stop("The vector lengths of 'num.col' + 'prop.col' must be the same as 'mu0' + 'p0'.")
  }
  if(any(p0 > 1 && po < 0)){
    stop("The 'p0' vector must contain values between 0 and 1.")
  }
  if(any(!alternative %in% c(1,2))){
    stop("Alternative must be '1 = greater/lesser' or '2 = equal/different'.")
  }
  if(alpha > 1 && alpha < 0){
    stop("The 'alpha' vector must contain values between 0 and 1.")
  }
  ### Tests
  tests <- data.frame()
  cols <- c(num.col,prop.col)
  if(length(labels) != length(cols)){
    warning("'labels' vector is not the same length as 'cols' vectors and will be ignored.")
  }
  tols <- c(mu0,p0); tols <- as.numeric(formatC(tols))
  for(i in seq_along(cols)){
    n <- length(sample.data[, cols[i]][!is.na(sample.data[, cols[i]])])
    if(cols[i] %in% num.col){
      point <- mean(sample.data[,cols[i]], na.rm = TRUE)
      sigma <- sd(sample.data[,cols[i]], na.rm = TRUE)
      std_err <- sigma / sqrt(n)
      if(alternative[i] == 1){
        other_data <- one_sided_mean(mu=point,mu0=tols[i],sigma=sigma,n=n,alpha=alpha)
      }else if(alternative[i] == 2){
        other_data <- two_sided_mean_equal(mu=point,mu0=tols[i],sigma=sigma,n=n,alpha=alpha)
      }
      tests <- rbind.data.frame(tests,
               cbind.data.frame(var=names(sample.data)[cols[i]],
               point=point,tolerance=tols[i],std_err=std_err,other_data))
    }else if(cols[i] %in% prop.col){
      point <- mean(sample.data[,cols[i]], na.rm = TRUE)
      pq <- point * (1 - point)
      std_err <- sqrt(pq / n)
      if(alternative[i] == 1){
        other_data <- one_sided_prop(p=point,p0=tols[i],n=n,alpha=alpha)
      }else if(alternative[i] == 2){
        other_data <- two_sided_prop_equal(p=point,p0=tols[i],n=n,alpha=alpha)
      }
      tests <- rbind.data.frame(tests,
               cbind.data.frame(var=names(sample.data)[cols[i]],
               point=point,tolerance=tols[i],std_err=std_err,other_data))
    }
  }
  is.num <- sapply(tests, is.numeric)
  tests[is.num] <- lapply(tests[is.num], round, 3)
  if(length(labels) == length(cols)){
    tests$var <- labels
  }
  return(tests)
}

one_sided_mean <- function(mu,mu0,sigma,n,alpha){
  z <- (mu-mu0)/(sigma/sqrt(n))
  pval <- round(pnorm(abs(z),lower.tail = FALSE),3)
  power <- pnorm(abs(z)-qnorm(1-alpha))
  return(data.frame(z_value=z,pval=pval,power=power))
}

one_sided_prop <- function(p,p0,n,alpha){
  z <- (p-p0)/sqrt(p0*(1-p0)/n)
  pval <- round(pnorm(abs(z),lower.tail = FALSE),3)
  power <- pnorm(sqrt(p0*(1-p0)/p/(1-p))*(abs(z)-qnorm(1-alpha)))
  return(data.frame(z_value=z,pval=pval,power=power))
}

two_sided_mean_equal <- function(mu,mu0,sigma,n,alpha){
  z <- (mu-mu0)/(sigma/sqrt(n))
  pval <- round(2*pnorm(-abs(z)),3)
  power <- pnorm(z-qnorm(1-alpha/2))+pnorm(-z- qnorm(1-alpha/2))
  return(data.frame(z_value=z,pval=pval,power=power))
}

two_sided_prop_equal <- function(p,p0,n,alpha){
  z <- (p-p0)/sqrt(p0*(1-p0)/n)
  pval <- round(2*pnorm(-abs(z)),3)
  power <- pnorm(z-qnorm(1-alpha/2))+pnorm(-z-qnorm(1-alpha/2))
return(data.frame(z_value=z,pval=pval,power=power))
}


