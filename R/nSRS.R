#' @title Calculates the sample size for SRS
#'
#' @description This function calculates de sample size in a SRS design based on one of three possible
#' arguments: \eqn{\sigma^2}, \eqn{pq} or \eqn{CV}. Provide just one of then.
#' @param moe Margin of error.
#' @param alpha 1 - (confidence level).
#' @param CV Variaton coefficient.
#' @param S2 \eqn{\sigma^2}, population variance.
#' @param pq \eqn{P(1 - P)}, population variance for proportions.
#' @param N Population size.
#' @return The value of \eqn{n}.
#' @export
#' @examples
#' nSRS(moe=0.05,alpha=0.05,S2=6,N=1500)
#' nSRS(moe=0.05,alpha=0.05,CV=0.2,N=1500)
#' nSRS(moe=0.05,alpha=0.05,pq=0.25,N=1500)

nSRS <- function(moe, alpha,CV=NULL,S2=NULL,pq=NULL,N=Inf){
  n.sam <- NULL
  if (alpha <= 0 || alpha >= 1)
    stop("alpha must be in (0,1).\n")
  z_alpha = qnorm(1 - (alpha / 2))
  if (sum(sapply(list(moe, N), is.null) != 0))
    stop("e and N cannot be NULL.\n")
  if (any(alpha > 1)) stop("alpha must be scalar.\n")

  #    if (any(N <= 0, S2 <= 0, CVpop <= 0))
  #        stop("None of N, S2, or CVpop can be <= 0\n")
  N.chk <- any(N <= 0, S2 <= 0, CV <= 0)
  if (N.chk) stop("None of N, S2, or CV can be <= 0\n")
  if(!is.null(S2)){
    n.sam=z_alpha^2*S2/(moe^2+z_alpha^2*S2/N)
  }

  if (!is.null(CV)){
    n.sam=z_alpha^2*CV^2/(moe^2+z_alpha^2*CV^2/N)
  }
  if(!is.null(pq)){
    n.sam=z_alpha^2*pq/moe^2
    if(!is.infinite(N)){
      n.sam=n.sam/(1 + (n.sam/N))
    }
  }
  if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
  return(ceiling(n.sam))
}
