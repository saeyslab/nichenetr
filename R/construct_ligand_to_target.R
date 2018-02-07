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
#' @return A list containing 3 elements (lr, sig, gr): the integrated weighted ligand-receptor, signaling and gene regulatory networks in data frame / tibble format.
#'
#'
#' @examples
#' ## Generate the weighted networks from input source networks
#' #' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
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
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1 ")

  requireNamespace("tidyverse")
  # remove data sources for which weight equals 0
  source_weights_df = source_weights_df %>% filter(weight > 0)
  # perform weighted network aggregation
  lr_network_w = lr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
  sig_network_w = sig_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
  gr_network_w = gr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()

  weighted_networks = list(lr = lr_network_w, sig = sig_network_w, gr = gr_network_w)
}
#' @title Add a new data source to the model
#'
#' @description \code{add_new_datasource} adds a new data source to one of the ligand-receptor, signaling and gene regulatory data sources.
#'
#' @usage
#' add_new_datasource(new_source, network, new_weight,source_weights_df)
#'
#' @param new_source A data frame / tibble containing novel interactions of the ligand-receptor, signaling or gene regulatory layer (required columns: from, to, source)
#' @param network NULL or a data frame / tibble containing the base network to which you want to add the new data source (required columns: from, to, source)
#' @param new_weight a weight value between 0 and 1 to assign to your new data source
#' @param source_weights_df A data frame / tibble containing the weights associated to each already included individual data source (required columns: source, weight).
#'
#' @return A list containing 2 elements (network and source_weights_df): the updated network containing your data source and the updated source_weights_df containing the weight of the newly added data source.
#'
#'
#' @examples
#' ## Update the lr_network with a new ligand-receptor data source
#' library(tidyverse)
#' lr_toy = tibble(from = "A", to = "B", source = "toy")
#' new_lr_network = add_new_datasource(lr_toy, lr_network,1,source_weights_df)
#'
#' @export
#'
add_new_datasource = function(new_source, network, new_weight,source_weights_df) {
  # input check
  if (!is.data.frame(new_source))
    stop("new_source must be a data frame or tibble object")
  if (!is.data.frame(network) && !is.null(network))
    stop("network must be a data frame or tibble object")
  if (!is.numeric(new_weight) && length(new_weight) != 1 && (new_weight >= 0 | new_weight <= 1))
    stop("new_weight must be a one number between 0 and 1")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0 || new_weight > 1)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")


  requireNamespace("tidyverse")

  if (is.null(network)){
    return(list(network = new_source, source_weights_df = tibble(source = new_source$source %>% unique(),weight = new_weight)))
  }
  network = network %>% bind_rows(new_source)
  source_weights_df = source_weights_df %>% bind_rows(tibble(source = new_source$source %>% unique(),weight = new_weight))

  list(network = network, source_weights_df = source_weights_df)
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

  requireNamespace("tidyverse")

  ligand_to_target = lr_network
}
