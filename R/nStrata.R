#' @title Calculates the sample size in a stratified desgin
#' @description This function returns the number of elements to be sampled in a stratified based survey.
#' If the method choosen is `prop` or `optimum`, the function return the allocation too.
#' @param data The data.frame containing the variable denoted as the cluster.
#' @param stratanames The variable name of the considered strata. Can be a vector of variable names.
#' @param alpha Confidence level.
#' @param moe Margin of error.
#' @param S2 The variance, as \eqn{\sigma^2}. Not necessary if `pq` or `V` is not `NULL`.
#' @param pq The variance, as \eqn{P(1-P)}. Not necessary if `S2` or `V` is not `NULL`.
#' @param V The variance, as \eqn{(d/t)^2}. Not necessary if `S2` or `pq` is not `NULL`.
#' @param N Population size.
#' @param method A string with `none`, `prop` or `optimum`. With `optimum`, it uses Neyman allocation.
#' @return The value (integer) of the number of elements to sample in stratified sampling. If `prop` or `optimum`
#' where used, the function returns a list with the components
#' \describe{
#' \item{n}{Total number of elements in the sample.}
#' \item{nh}{In case of method `prop` or `optimum`, the allocation in each strata.}
#' }
#' @export
#' @examples
#' library(survey)
#' data("api")
#' nStrata(apipop,stratanames = "stype", alpha = 0.95, moe = 0.05,
#'         pq = 0.25, N = Inf, method = "none")
#' nStrata(apipop,stratanames = "stype", alpha = 0.95, moe = 0.05,
#'         pq = 0.25, N = nrow(apipop), method = "prop")
#' nStrata(apipop,stratanames = "stype", alpha = 0.95, moe = 0.05,
#'         pq = 0.25, N = nrow(apipop), method = "optimum")
#' ### Cochran's example (Cochran, W. G. (2007). Sampling techniques. John Wiley & Sons. p. 106-107)
#' data_cochran=data.frame(stratum=rep(c(1:6),c(13,18,26,42,73,24)))
#' nStrata(data_cochran,stratanames = "stratum", N=196, method = "optimum",
#' S2=c(325^2,190^2,189^2,82^2,86^2,190^2),V=7974976)




nStrata <- function(data,stratanames,alpha,moe,S2=NULL,pq=NULL,V=NULL,N=Inf,
                    method=c("none","prop","optimum")){
  if(is.null(S2) && is.null(pq) && is.null(V)){
    stop("You must provide 'V', 'S2' or 'pq' argument.")
  }
  if(missing(method)) {warning("The method is not specified; by default, the method is 'none'.")
    method="none"
  }
  if(!(method %in% c("none","prop","optimum"))){
    stop("The method name is not in the list. You must input 'none', 'prop' or 'optimum'.")}
  if(method %in% c("prop","optimum") && N==Inf){
    stop("You have to input N if method is 'prop' or 'optimum'.")
  }
  dados=as.data.frame(data)
  if(!any(stratanames %in% names(dados))){
    stop("One or more of your stratanames are wrong.")
  }
  m=match(stratanames,colnames(dados))
  dados$stratum=do.call(paste, c(dados[m], sep = "_"))
  n_strata=length(unique(dados$stratum))
  if(!is.null(S2) && length(S2) != 1 && length(S2) != n_strata){
    stop("The S2 vector must be of the same length as the number os strata. \n
         In this case, your strata number is ",n_strata,".")
  }
  if(!is.null(pq) && length(pq) != 1 && length(pq) != n_strata){
    stop("The S2 vector must be of the same length as the number os strata. \n
         In this case, your strata number is ",n_strata,".")
  }
  if(!is.null(S2)){var=S2
  if(length(S2)==1){var=rep(S2,n_strata)}
  }else{var=pq
  if(length(pq)==1){var=rep(pq,n_strata)}}

  if(is.null(V)){
    z_alpha=qnorm(1 - ((1 - alpha) / 2))
    V=(moe/z_alpha)^2}
  Wh=table(dados$stratum)/sum(table(dados$stratum))
  Nh=table(dados$stratum)
  Wh2var2=Wh^2*var
  wh=Wh *sqrt(var)/sum(Wh*sqrt(var))
  if(method=="none"){
    if(!is.null(S2)){
        n_linha=(1/V)*sum(Wh2var2/wh)
        n=ceiling(n_linha/(1+(1/(N*V))*sum(Wh*var)))
        # n=ceiling(sum(WhSh)^2/(moe^2/(z_alpha^2)*(1/N)*sum(Wh*var^2)))
        if(n < (length(n_strata) * 5)){n=n_strata * 5}
        if(n_linha < (length(n_strata) * 5)){n_linha=n_strata * 5}
        if(is.infinite(N)){return(ceiling(n_linha))}
        if(!is.infinite(N)){return(n)}
      }else{
        wh=Wh*sqrt(var)
        n_linha=(z_alpha^2/moe^2)*sum(Wh*var/wh)
        n=ceiling(n_linha/(1+(n_linha/N)))
        if(n < (length(n_strata) * 5)){n=n_strata * 5}
        if(n_linha < (length(n_strata) * 5)){n_linha=n_strata * 5}
        if(is.infinite(N)){return(ceiling(n_linha))}
        if(!is.infinite(N)){return(n)}
      }
    }
  if(method=="prop"){
    if(!is.null(S2)){
      n_linha=sum(Wh * var)/V
      n=ceiling(n_linha/(1+(n_linha/N)))
      if(n < (length(n_strata) * 5)){n=n_strata * 5}
      nh=round(n*prop.table(table(dados$stratum)),0)
      if(any(nh<2)){
        id=which(nh<2)
        nh[id]=2
        n=sum(nh)
      }
      return(structure(list(n=n,nh=nh)))
    }else{
      n_linha=(z_alpha^2 *sum(Wh*var))/moe^2
      n=ceiling(n_linha/(1+(n_linha/N)))
      if(n < (length(n_strata) * 5)){n=n_strata * 5}
      nh=round(n*prop.table(table(dados$stratum)),0)
      if(any(nh<2)){
        id=which(nh<2)
        nh[id]=2
        n=sum(nh)
      }
      return(structure(list(n=n,nh=nh)))
      }
    }
  if(method=="optimum"){
    if(!is.null(S2)){
      NhSh <- Nh * sqrt(var)
      n_linha=sum(NhSh)^2/V
      n=ceiling(n_linha/(1 + (1/V) * sum(Nh * var)))
    }else{
      n_linha=(z_alpha^2*sum(Wh*sqrt(var))^2)/moe^2
      n=n_linha/(1+(z_alpha^2/moe^2)*sum(Wh*var))
    }

    if(n < (length(n_strata) * 5)){n=n_strata * 5}
    neyman=(Wh*var)/sum(Wh*var)
    nh=round(n*neyman,0)
    return(structure(list(n=n,nh=nh)))
  }
}
