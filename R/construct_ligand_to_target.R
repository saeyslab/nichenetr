#' @title Construct weighted layer-specific networks
#'
#' @description \code{construct_weighted_networks} construct layer-specific weighted integrated networks from input source networks via weighted aggregation.
#'
#' @usage
#' construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#'
#' @param lr_network A data frame / tibble containing ligand-receptor interactions (required columns: from, to, source)
#' @param sig_network A data frame / tibble containing signaling interactions (required columns: from, to, source)
#' @param gr_network A data frame / tibble containing gene regulatory interactions (required columns: from, to, source)
#' @param source_weights_df A data frame / tibble containing the weights associated to each individual data source. Sources with higher weights will contribute more to the final model performance (required columns: source, weight). Note that only interactions described by sources included here, will be retained during model construction.
#'
#' @return A list containing 3 elements: the integrated weighted ligand-receptor, signaling and gene regulatory networks in data frame / tibble format.
#'
#'
#' @examples
#' ## Generate the weighted networks from input source networks
#' library(tidyverse)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#'
#'
#' @export
#'
construct_weighted_networks = function(lr_network, sig_network, gr_network,source_weights_df) {
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.data.frame(source_weights_df))
    stop("source_weights_df must be a data frame or tibble object")
  # remove data sources with weight equals 0
  source_weights_df = source_weights_df %>% filter(weight > 0)
  # perform weighted network aggregation
  lr_network_w = lr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight))
  sig_network_w = sig_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight))
  gr_network_w = gr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight))

  weighted_networks = list(lr = lr_network_w, sig = sig_network_w, gr = gr_network_w)
}

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
#' @importFrom igraph page_rank
#'
#' @examples
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' construct_ligand_target_matrix(lr_network, sig_network, gr_network)
#'
#' @export
#'
construct_ligand_target_matrix = function(lr_network, sig_network, gr_network) {
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
