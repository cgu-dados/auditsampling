#' @title Calculates the remaing sample size when using pilot sample
#' @description Based on a pilot sample, this function calculates the remain elements to the final sample
#' given the variance considered.
#' @param S2 The variance, as \eqn{\sigma^2}. Not necessary if `S2` is not `NULL`.
#' @param pq The variance, as \eqn{P(1-P)}. Not necessary if `S2` is not `NULL`.
#' @param n1 The size of the pilot sample.
#' @param V The target variance.
#' @return The final sample is
#' \deqn{n = n_1 + n_2}
#' Where \eqn{n_2} is the value returned by this function.
#' @references COCHRAN, William Gemmell. Sampling techniques-3. 1977.
#' @export
#' @examples
#' nPilot(s2=13.5,n1=30,V=1)
#' nPilot(pq=0.15,n1=30,V=0.025)


nPilot <- function(s2=NULL,pq=NULL,n1,V){
  if(is.null(s2) && is.null(pq)){
    stop('You must provide s2 or pq argument.')
  }
  if(!is.null(s2)){
    var=s2
    n=(var/V)*(1 + (2 / n1))
  }else{
    var=pq
    n=(pq / V)+((3-(8*pq))/pq)+((1 - (3 * pq))/(V*n1))
    }
  return(ceiling(n))
}