#' @name siasgdw
#' @title Data from brazillian government bid processes
#'
#' @description Data from Sistema Integrado de Administracao de Servicos Gerais (SIASG) from goods and services purchases by the Brazilian government
#'
#' A dataset containing all bids made in 2018
#' It has some simulated data, to be used for learning purposes.
#'
#' @docType data
#'
#' @usage data(siasgdw)
#'
#' @format A dataset of bids made by Brazillian government offices, totalizing 164372 observations,
#' on the  following variables
#' \describe{
#' \item{id_compra}{ID of the bidding process}
#' \item{no_processo}{Number of the processo protocol}
#' \item{dt_compra}{Reference date of the purchase}
#' \item{id_unidde}{ID of the department of a public agency responsible for the purchase}
#' \item{no_unidade}{Name of the department a public agency responsible for the purchase}
#' \item{orgao}{Name of the public agency to which the one responsible for the purchase is subordinate}
#' \item{orgao_sup}{Name of the supervisor public agency of the respective responsible for the purchase}
#' \item{modalidae}{The category of juridic rules and procedures involved in the bidding process}
#' \item{valor_compra}{The purchase amount involved in the bidding process}
#' \item{qtd_forn}{The number of bidders involved in the bidding process}
#' \item{qtd_itens}{The number of different goods and/or services demanded in the bidding process}
#' \item{controle_01 / 12}{Internal controls failure occurrences (binary)}
#' \item{risco_01 / 05}{Risk occurrences (binary)}
#' \item{risco_06 / 08}{Risk occurrences (financial losses)}
#' \item{tempo_total}{Total time spent, in days, for the conclusion of the bidding process}
#' \item{tempo_fase_interna}{Total time spent, in days, for the planning of the bidding process}
#' }
#'
#' @keywords datasets
#'
#' @source Sistema Integrado de Administração de Serviços Gerais(SIASG)
#'
#' @examples
#' library(auditsampling)
#' data('siasgdw')
#' ############# Subscripting
#' uasg_cgu <- siasgdw[siasgdw$id_unidade == 370003,]
#' summary(uasg_cgu)
"siasgdw"
