#' @title Construct ligand-target matrix
#'
#' @description \code{evaluate_ligand_prediction} construct a ligand-target matrix from input (weighted) networks.
#'
#' @usage
#' evaluate_ligand_prediction(lr_network, sig_network, gr_network)
#'
#' @param lr_network A list.
#' @param sig_network A list.
#' @param gr_network A list.
#'
#' @return A matrix containing the probability scores of ligand-target links.
#'
#' @importFrom igraph page_rank
#'
#' @examples
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' evaluate_ligand_prediction(lr_network, sig_network, gr_network)
#'
#' @export
#'
evaluate_ligand_prediction <- function(lr_network, sig_network, gr_network) {
  # requireNamespace("igraph")
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")

  ligand_to_target = lr_network
}
#' @title Construct ligand-target matrix
#'
#' @description \code{evaluate_target_prediction} construct a ligand-target matrix from input (weighted) networks.
#'
#' @usage
#' evaluate_target_prediction(lr_network, sig_network, gr_network)
#'
#' @param lr_network A list.
#' @param sig_network A list.
#' @param gr_network A list.
#'
#' @return A matrix containing the probability scores of ligand-target links.
#'
#' @importFrom igraph page_rank
#'
#' @examples
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' evaluate_target_prediction(lr_network, sig_network, gr_network)
#'
#' @export
#'
evaluate_target_prediction <- function(lr_network, sig_network, gr_network) {
  # requireNamespace("igraph")
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")

  ligand_to_target = lr_network
}
