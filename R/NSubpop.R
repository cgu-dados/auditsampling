#' @title Calculates the new subpopulation due elements thar not belong to the populational universe
#' @description This function returns the number of the new subpopulation.
#' @param data The data.frame containing the variable denoted as the cluster.
#' @param N Population size.
#' @param n Sample size.
#' @param p Proportion or the number of elements that belong to the universe intended.
#' @param z The z-value of \eqn{\alpha}.
#' @return The value (integer) of the new population (subpopulation) number.
#' @export
#' @examples
#' NSubpop(45000,100,59)
#' NSubpop(45000,100,5)

NSubpop <- function(N,n,p,z=3){
  if(p > 1){p = p / n}
  fpc = (N - n) / (N - 1)
  lim_sup = p + (z * sqrt((p * (1 - p) / n) * fpc))
  return(ceiling(N * lim_sup))
}
