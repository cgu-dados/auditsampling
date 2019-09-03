#' @title Calculates the remaing sample size when using pilot sample
#' @description Based on a pilot sample, this function calculates the remain elements to the final sample
#' given the variance considered.
#' @param s2 The variance, as \eqn{\sigma^2}. Not necessary if \code{pq} is not \code{NULL}.
#' @param pq The variance, as \eqn{P(1-P)}. Not necessary if \code{s2} is not \code{NULL}.
#' @param n1 The size of the pilot sample.
#' @param alpha 1 - confidence level.
#' @param moe Margin of error.
#' @return The final sample is
#' \deqn{n = n_1 + n_2}
#' Where \eqn{n_2} is the value returned by this function.
#' @references COCHRAN, William Gemmell. Sampling techniques-3. 1977.
#' @export
#' @examples
#' nPilot(s2=13.5,n1=30,alpha=0.05,moe=5)
#' nPilot(pq=0.15,n1=30,alpha=0.1,moe=0.05)


nPilot <- function(s2=NULL, pq=NULL, n1, alpha, moe, N){
  V=(moe/qnorm(1 - (alpha / 2)))^2
  if(is.null(s2) && is.null(pq)){
    stop('You must provide s2 or pq argument.')
  }
  if(!is.null(s2)){
    var=s2
    n_linha=(var/V)*(1 + (2 / n1))
    n_linha=n_linha+n1
    n2=n_linha/(1 + (n_linha / N))
    n=n2-n1
  }else{
    var=pq
    n_linha=(pq / V)+((3-(8*pq))/pq)+((1 - (3 * pq))/(V*n1))
    n_linha=n_linha+n1
    n2=n_linha/(1 + (n_linha / N))
    n=n2-n1
    }
  return(ceiling(n))
}
