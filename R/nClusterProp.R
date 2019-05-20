#' @title Calculates the number os cluster to a cluster sample design
#' @description This function returns the number of clusters to be sampled in a proportion based survey.
#' By default, the argument `pq` is 0.25 (maximum variance).
#' @param data The data.frame containing the variable denoted as the cluster.
#' @param clustername The variable name of the considered cluster.
#' @param alpha Confidence level.
#' @param moe Margin of error.
#' @param pq The variance, as \eqn{P(1-P)}. By default, it is 0.25.
#' @return The value (integer) of the number of clusters to sample in cluster sampling.
#' @export
#' @examples
#' library(survey)
#' data("api")
#' nClusterProp(data=apipop, clustername = "dnum", alpha = 0.95, moe = 0.1)
#' nClusterProp(data=apipop, clustername = "dnum", alpha = 0.95, moe = 0.05)

nClusterProp <- function(data,clustername,alpha,moe,pq=0.25){
  m_linha=mean(table(data[,clustername]))
  z_alpha=qnorm(1 - ((1 - alpha) / 2))
  n=(z_alpha^2*pq)/(m_linha*moe^2)
  return(ceiling(n))
}
