#' @title Calculates the sample size for hypotesis tests
#'
#' @description This function calculates de sample size for hypotesys tests considering the following hypotesis:
#' \deqn{H_{0}:\text{ } \mu \neq \mu_{0}}
#' \deqn{H_{0}:\text{ } \bar{p} \neq p_{0}}
#' \deqn{H_{0}:\text{ } \mu \leq \text{ or } \geq \mu_{0}}
#' \deqn{H_{0}:\text{ } \bar{p} \leq \text{ or } \geq p_{0}}
#' @param p Population proportion.
#' @param p0 Hypothesized proportion.
#' @param mu Population mean.
#' @param mu0 Hypothesized mean.
#' @param alpha Level of significance.
#' @param alternative \code{1} if is a greater/lesser test and \code{2} for equal/different one.
#' @return The sample size.
#' @export
#' @examples
#' nPowerSample(p=0.03, p0=0.05)
#' nPowerSample(p=0.02, p0=0.05)
#' nPowerSample(mu=100, mu0=120,sd=24)
#' nPowerSample(mu=100, mu0=110,sd=24,beta=0.05,alternative = 2)


nPowerSample <- function(p=NULL,p0=NULL,mu=NULL,mu0=NULL,sd=NULL,alpha=0.05,beta=0.2,alternative=1){
  if(!alternative %in% c(1, 2)){
    stop("Alternative must be '1 = greater/lesser' or '2 = equal/different'.")
  }
  if(any(!is.null(p),!is.null(p0)) && any(!is.null(mu),!is.null(mu0))){
    stop("Input p or mu, not both.")
  }
  if(alternative == 1){
    if(!is.null(p) && !is.null(p0)){
      (n=p0*(1-p0)*((qnorm(1-alpha)+qnorm(1-beta)*sqrt(p*(1-p)/p0/(1-p0)))/(p-p0))^2)
      return(ceiling(n))
    }
    if(!is.null(mu) && !is.null(mu0)){
      (n=(sd*(qnorm(1-alpha)+qnorm(1-beta))/(mu-mu0))^2)
      return(ceiling(n))
    }
  }
    if(alternative == 2){
      if(!is.null(p) && !is.null(p0)){
        (n=p*(1-p)*((qnorm(1-alpha/2)+qnorm(1-beta))/(p-p0))^2)
        return(ceiling(n))
      }
      if(!is.null(mu) && !is.null(mu0)){
        (n=(sd*(qnorm(1-alpha/2)+qnorm(1-beta))/(mu-mu0))^2)
        return(ceiling(n))
      }
  }
}
