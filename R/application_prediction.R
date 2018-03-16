#' @title Convert cluster assignment to settings format suitable for target gene prediction.
#'
#' @description \code{convert_cluster_to_settings} Convert cluster assignment to settings format suitable for target gene prediction.
#'
#' @usage
#' convert_cluster_to_settings = function(i, cluster_vector, setting_name, setting_from, background = NULL)
#'
#' @param i The cluster number of the cluster of interest to which genes should belong
#' @param cluster_vector Named vector containing the cluster number to which every gene belongs
#' @param setting_name Base name of the setting
#' @param setting_from Active ligands for the specific setting
#' @param background NULL or a character vector of genes belonging to the background. When NULL: the background will be formed by genes belonging to other clusters that the cluster of interest. Default NULL.
#'
#' @return A list with following elements: $name (indicating the cluster id), $from, $response. $response is a gene-named logical vector indicating whether the gene is part of the respective cluster.
#'
#' @examples
#'
#' genes_clusters = c("TGFB1" = 1,"TGFB2" = 1,"TGFB3" = 2)
#' cluster_settings = lapply(seq(length(unique(genes_clusters))), convert_cluster_to_settings, cluster_vector = genes_clusters, setting_name = "example", setting_from = "BMP2")
#'
#' @export
#'
convert_cluster_to_settings = function(i, cluster_vector, setting_name, setting_from, background = NULL){

  # input check
  if(!is.numeric(i) | length(i) != 1 | i <= 0)
    stop("i should be a number higher than 0")
  if(!is.numeric(cluster_vector) | is.null(names(cluster_vector)))
    stop("cluster_vector should be a named numeric vector")
  if(!is.character(setting_name))
    stop("setting_name should be a character vector")
  if(!is.character(setting_from))
    stop("setting_from should be a character vector")
  if(!is.character(background) & !is.null(background))
    stop("background should be a character vector or NULL")

  requireNamespace("dplyr")

  genes_cluster_oi = cluster_vector[cluster_vector == i] %>% names()
  if (is.null(background)){
    response = names(cluster_vector) %in% genes_cluster_oi
    names(response) = names(cluster_vector)
  } else {
    background_logical = rep(FALSE,times = length(background))
    names(background_logical) = background
    cluster_logical = rep(TRUE,times = length(genes_cluster_oi))
    names(cluster_logical) = genes_cluster_oi
    response = c(background_logical,cluster_logical)
  }
  return(list(name = paste0(setting_name,"_cluster_",i), from = setting_from, response = response))
}
