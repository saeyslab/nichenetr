#' @title Construct ligand-target matrix
#'
#' @description \code{construct_ligand_target_matrix} construct a ligand-target matrix from input (weighted) networks.
#'
#' @usage
#' construct_ligand_target_matrix(lr_network, sig_network, gr_network)
#'
#' @param lr_network A list.
#' @param sig_network A list.
#' @param gr_network A list.
#'
#' @return A matrix containing the probability scores of ligand-target links.
#'
#' @importFrom igraph page_rank ajdacency_to_dataframe
#'
#' @example
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' construct_ligand_target_matrix(lr_network, sig_network, gr_network)
#'
#' @export
#'
construct_ligand_target_matrix <- function(lr_network, sig_network, gr_network) {
  requireNamespace("igraph")
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")

  ligand_to_target
}
