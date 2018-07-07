#' @title Get the complete ligand-receptor, signaling and gene regulatory network.
#'
#' @description \code{get_all_networks} will return a list containing the complete ligand-receptor network, signaling network and gene regulatory network.
#'
#' @usage
#' get_all_networks()
#'
#' @return a list containing the following elements: $lr: The true complete ligand-receptor network, $sig: the true complete signaling network, $gr: the true complete gene regulatory network.
#'
#' @examples
#' \dontrun{
#' all_networks = get_all_networks()
#' }
#' @export
#'
get_all_networks = function(){
  load("R/sysdata.rda")
  return(list(lr = lr_network, sig = sig_network_real, gr = gr_network_real))
}
#' @title Get the complete set of ligand-treatment expression datasets.
#'
#' @description \code{get_expression_settings_validation} will return the complete set of ligand-treatment expression datasets.
#'
#' @usage
#' get_expression_settings_validation()
#'
#' @return The true complete set of ligand-treatment expression datasets.
#'
#' @examples
#' \dontrun{
#' expression_settings_validation = get_expression_settings_validation()
#' }
#' @export
#'
get_expression_settings_validation = function(){
  load("R/sysdata.rda")
  return(expression_settings_validation_real)
}
