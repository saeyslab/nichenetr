#' @title Randomize a network
#'
#' @description \code{randomize_network} Randomizes a network of interest by edge swapping to preserve the degree of the nodes
#'
#' @usage
#' randomize_network(network,output_weighted = FALSE)
#'
#' @param network A data frame / tibble containing gene-gene interactions (required columns: $from, $to)
#' @param output_weighted Indicate whether the output network should be made weighted by assigning a weight of 1 to each interaction.
#'
#' @return A randomized network ($from, $to; and $weight = 1 if output_weighted == TRUE).
#'
#' @importFrom igraph graph_from_adjacency_matrix degree sample_degseq get.data.frame
#'
#' @examples
#' \dontrun{
#' random_lr = randomize_network(lr_network)
#' }
#'
#' @export
#'
randomize_network = function(network, output_weighted = FALSE){

  if (!is.data.frame(network))
    stop("network must be a data frame or tibble object")

  requireNamespace("dplyr")

  allgenes = unique(c(network$from, network$to)) %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_entrez_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  entrez2allgenes = mapper(allgenes_entrez_tbl,"allgenes_integer","allgenes")
  allgenes2entrez = mapper(allgenes_entrez_tbl,"allgenes","allgenes_integer")

  network = network %>% mutate(from_allgenes = entrez2allgenes[from], to_allgenes = entrez2allgenes[to])

  network_matrix = Matrix::sparseMatrix(network$from_allgenes %>% as.integer, network$to_allgenes %>% as.integer, x=1 %>% as.numeric, dims = c(length(allgenes), length(allgenes))) # cast to sparse matrix
  network_igraph = igraph::graph_from_adjacency_matrix(network_matrix, mode="directed")
  network_random_graph = igraph::sample_degseq(igraph::degree(network_igraph,mode = "out"), igraph::degree(network_igraph,mode = "in"))
  network_random_df = igraph::get.data.frame(network_random_graph) %>% as_tibble() %>% mutate(from = as.character(allgenes2entrez[from]), to = as.character(allgenes2entrez[to]))

  if (output_weighted == TRUE){
    network_random_df = network_random_df %>% mutate(weight = 1)
  }
  return(network_random_df)
}
#' @title Randomize a network of a particular data source.
#'
#' @description \code{randomize_datasource_network} Randomizes a network of a data source of interest by edge swapping to preserve the degree of the nodes.
#'
#' @usage
#' randomize_datasource_network(datasource, network)
#'
#' @param datasource The name of the data source for which the interactions need to be shuffled.
#' @param network A data frame / tibble containing gene-gene interactions (required columns: $from, $to, $source)
#'
#' @return A randomized network ($from, $to, $source)
#'
#' @importFrom igraph graph_from_adjacency_matrix degree sample_degseq get.data.frame
#'
#' @examples
#' \dontrun{
#' datasource_lr = lr_network$source[1]
#' lr_randomized_source = randomize_datasource_network(datasource_lr, lr_network)
#'}
#'
#' @export
#'
randomize_datasource_network = function(datasource,network){

  requireNamespace("dplyr")

  network = network %>% filter(source == datasource)

  allgenes = unique(c(network$from, network$to)) %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_entrez_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  entrez2allgenes = mapper(allgenes_entrez_tbl,"allgenes_integer","allgenes")
  allgenes2entrez = mapper(allgenes_entrez_tbl,"allgenes","allgenes_integer")

  network = network %>% mutate(from_allgenes = entrez2allgenes[from], to_allgenes = entrez2allgenes[to])

  network_matrix = Matrix::sparseMatrix(network$from_allgenes %>% as.integer, network$to_allgenes %>% as.integer, x=1 %>% as.numeric, dims = c(length(allgenes), length(allgenes))) # cast to sparse matrix
  network_igraph = igraph::graph_from_adjacency_matrix(network_matrix, mode="directed")
  network_random_graph = igraph::sample_degseq(igraph::degree(network_igraph,mode = "out"), igraph::degree(network_igraph,mode = "in"))
  network_random_df = igraph::get.data.frame(network_random_graph) %>% as_tibble() %>% mutate(from = as.character(allgenes2entrez[from]), to = as.character(allgenes2entrez[to])) %>% mutate(source = datasource)

}
#' @title Randomize an integrated network by shuffling its source networks
#'
#' @description \code{randomize_complete_network_source_specific} Randomizes an integrated network of interest by edge swapping in the source-specific networks to preserve the degree of the nodes
#'
#' @usage
#' randomize_complete_network_source_specific(network)
#'
#' @param network A data frame / tibble containing gene-gene interactions (required columns: $from, $to, $source)
#'
#' @return A randomized network ($from, $to, $source).
#'
#' @importFrom igraph graph_from_adjacency_matrix degree sample_degseq get.data.frame
#'
#' @examples
#'\dontrun{
#' random_lr = randomize_complete_network_source_specific(lr_network)
#'}
#'
#' @export
#'
randomize_complete_network_source_specific = function(network){
  network = lapply(unique(network$source),randomize_datasource_network,network) %>% bind_rows()
}
