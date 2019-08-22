#' @title Returns a data.frame with total resources to be spend in sampling process, based on the thresholds defined
#' @description This function returns a data.frame with the calculations of time to be spend on sampling process, listing respectives sampling parameters, like confidence levels and margin of errors.
#' Need a sample size function to work.
#' @param hh The cost, measured in time (HH), for collecting data for each sample unit.
#' @param max.hh The cost threshold from which results and their parameters will be excluded.
#' @param max.n The sample size threshold from which results and their parameters will be excluded.
#' @param max.alpha The maximum z-value of \eqn{\alpha} from which results and their parameters will be excluded.
#' @param max.beta The maximum z-value of \eqn{\beta} from which results and their parameters will be excluded. Applicable only to \code{nPowerSample} function.
#' @param fun Function used to calculate samplesize. Can be \code{nSRS}, \code{nStrata} or \code{nSampleSize}.
#' @param ... Others parameters passed to be used with the sample size functions. See examples.
#' @return A data.frame with
#' \describe{
#' \item{alpha}{The \eqn{\alpha}.}
#' \item{moe}{The margin of error, when using \code{nSRS} or \code{nStrata}.}
#' \item{beta}{The \eqn{\alpha}, when using \code{nSampleSize}.}
#' \item{n}{Respective sample size.}
#' \item{hh}{Cost, in time, as the argument passed to the function.}
#' \item{hh_total}{Total cost, in time.}
#' }
#' @export
#' @examples
#' ################### nSRS
#' nPlan(hh = (10/60), fun = "nSRS", max.alpha = 0.1, max.n = 176, max.hh = 60,
#'       moe = seq(0.01,0.1,by=0.01), pq = 0.25, N = 16548)
#'
#' ################### nPowerSample
#' nPlan(hh = (5/60), fun = "nPowerSample", max.hh = 40, max.alpha = 0.1, max.beta = 0.2,
#'       p=0.02, p0=0.05, alternative = 1)


nPlan <- function(hh, max.hh = NULL, max.n = NULL, max.alpha = NULL,
                  max.beta = NULL, fun = "nSRS",...){
  if(!fun %in% c("nSRS","nStrata","nPowerSample")){
    stop("'fun' argument must be 'nSRS', 'nStrata' or 'nPowerSample'.")
  }
  alpha = seq(0.01,0.4,by=0.01)
  args = list(...)
  if(fun == "nSRS"){
    grid = expand.grid(alpha,args$moe)
    names(grid) <- c("alpha","moe")
    grid$n <- nSRS(alpha = grid$alpha, moe = grid$moe, pq = args$pq, CV = args$CV,
                   S2 = args$S2, N = args$N)
    grid$hh <- hh
    grid$hh_total <- grid$n * grid$hh
    if(!is.null(max.hh)) grid <- grid[grid$hh_total <= max.hh, ]
    if(!is.null(max.n)) grid <- grid[grid$n <= max.n, ]
    if(!is.null(max.alpha)) grid <- grid[grid$alpha <= max.alpha, ]
    grid <- grid[order(grid$alpha,grid$moe,-grid$hh_total),]
  }
  if(fun == "nStrata"){
    grid = expand.grid(alpha,args$moe)
    names(grid) <- c("alpha","moe")
    if(args$method != "none"){
      grid$n <- nStrata(data = args$data, alpha = grid$alpha,
                        moe = grid$moe, pq = args$pq, V = args$V, stratanames = args$stratanames,
                        S2 = args$S2, N = args$N, method = args$method)$`n`
    }else{
      grid$n <- nStrata(data = args$data, alpha = grid$alpha,
                        moe = grid$moe, pq = args$pq, V = args$V, stratanames = args$stratanames,
                        S2 = args$S2, N = args$N, method = args$method)
    }
    grid$hh <- hh
    grid$hh_total <- grid$n * grid$hh
    if(!is.null(max.hh)) grid <- grid[grid$hh_total <= max.hh, ]
    if(!is.null(max.n)) grid <- grid[grid$n <= max.n, ]
    if(!is.null(max.alpha)) grid <- grid[grid$alpha <= max.alpha, ]
    grid <- grid[order(grid$alpha,grid$moe,-grid$hh_total),]
  }
  if(fun == "nPowerSample"){
    beta = seq(0.01,0.4,by=0.01)
    grid = expand.grid(alpha,beta)
    names(grid) <- c("alpha","beta")
    grid$n <- nPowerSample(p = args$p, p0 = args$p0, mu = args$mu, mu0 = args$mu0,
                           sd = args$sd, alpha = grid$alpha, beta = grid$beta,
                           alternative = args$alternative)
    grid$hh <- hh
    grid$hh_total <- grid$n * grid$hh
    if(!is.null(max.hh)) grid <- grid[grid$hh_total <= max.hh, ]
    if(!is.null(max.n)) grid <- grid[grid$n <= max.n, ]
    if(!is.null(max.alpha)) grid <- grid[grid$alpha <= max.alpha, ]
    if(!is.null(max.beta)) grid <- grid[grid$beta <= max.beta, ]
    grid <- grid[order(grid$alpha,grid$beta,-grid$hh_total),]
  }
  return(grid)
}


