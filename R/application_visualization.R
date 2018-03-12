#' @title Get ligand-target signaling paths between ligand(s) and target gene(s) of interest
#'
#' @description \code{get_ligand_signaling_path} Extract possible signaling paths between a ligand and target gene of interest. The most highly weighted path(s) will be extracted.
#'
#' @usage
#' get_ligand_signaling_path(ligand_tf_matrix, ligands_all, targets_all, top_n_regulators = 4, weighted_networks, ligands_position = "cols")
#'
#' @param ligand_tf_matrix A matrix of ligand-regulator probability scores
#' @param ligands_all A character vector of one or more ligands of interest
#' @param targets_all A character vector of one or more target genes of interest
#' @param top_n_regulators The number of top regulators that should be included in the ligand-target signaling network. Top regulators are regulators that score both high for being upstream of the target gene(s) and high for being downstream of the ligand. Default: 4.
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols".
#'
#' @return A list containing 2 elements (sig and gr): the integrated weighted ligand-signaling and gene regulatory networks data frame / tibble format with columns: from, to, weight
#'
#'
#' @examples
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' @export
#'
get_ligand_signaling_path = function(ligand_tf_matrix, ligands_all, targets_all, top_n_regulators = 4, weighted_networks, ligands_position = "cols"){

  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")
  if(!is.matrix(ligand_tf_matrix))
    stop("ligand_tf_matrix should be a matrix")
  if (ligands_position == "cols"){
    if(sum((ligands_all %in% colnames(ligand_tf_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_tf_matrix")
  } else if (ligands_position == "rows") {
    if(sum((ligands_all %in% rownames(ligand_tf_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_tf_matrix")
  }
  if(sum((targets_all %in% unique(c(weighted_networks$gr$to))) == FALSE) > 0)
    stop("target genes should be in gene regulatory network")
  if(!is.numeric(k) | length(k) != 1 | top_n_regulators <= 0)
    stop("k should be a number higher than 0")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  requireNamespace("dplyr")

  final_combined_df  = construct_ligand_signaling_df(ligands_all,targets_all,top_n_regulators,weighted_networks,ligand_tf_matrix)

  signaling_network_all = weighted_networks$lr_sig %>% mutate(weight = 1/weight) # inverse weight to prepare for SPL
  signaling_igraph = igraph::graph_from_data_frame(signaling_network_all, directed = TRUE)

  tf_nodes = lapply(ligands_all,get_shortest_path_signaling, final_combined_df,signaling_igraph) %>% unlist() %>% unique()

  tf_signaling = weighted_networks$lr_sig %>%
    filter(from %in% c(ligands_all,tf_nodes) & to %in% tf_nodes) %>% group_by(from,to) %>% mutate(weight = sum(weight)) %>% ungroup() %>% distinct()
  tf_regulatory = weighted_networks$gr %>%
    filter(from %in% final_combined_df$TF & to %in% targets_all)  %>% ungroup() %>% distinct()

  return(list(sig = tf_signaling, gr = tf_regulatory))
}
#' @title Prepare extracted ligand-target signaling network for visualization with DiagrammeR.
#'
#' @description \code{diagrammer_format_signaling_graph} Prepare extracted ligand-target signaling network for visualization with DiagrammeR.
#'
#' @usage
#' diagrammer_format_signaling_graph(signaling_graph_list, ligands_all,targets_all, sig_color = "steelblue", gr_color = "orange")
#'
#' @param signaling_graph_list A list of two elements: sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param sig_color The color for ligand-signaling edges and the ligand node. Default: steelblue.
#' @param gr_color The color for the gene regulatory edges and the target node. Default: orange.
#' @inheritParams get_ligand_signaling_path
#'
#' @return A DiagrammeR Graph object ready for visualization via DiagrammeR::render_graph.
#'
#' @importFrom DiagrammeR create_node_df create_edge_df create_graph
#'
#' @examples
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' graph = diagrammer_format_signaling_graph(ligand_target_signaling_list, all_ligands,all_targets)
#' # DiagrammeR::render_graph(graph, layout = "tree")
#' @export
#'
diagrammer_format_signaling_graph = function(signaling_graph_list, ligands_all,targets_all, sig_color = "steelblue", gr_color = "orange"){
  if (!is.list(signaling_graph_list))
    stop("signaling_graph_list must be a list object")
  if (!is.data.frame(signaling_graph_list$sig))
    stop("signaling_graph_list$sig must be a data frame or tibble object")
  if (!is.data.frame(signaling_graph_list$gr))
    stop("signaling_graph_list$gr must be a data frame or tibble object")

  if(sum((ligands_all %in% unique(c(signaling_graph_list$sig$from))) == FALSE) > 0)
    stop("ligands should be in signaling_graph_list")
  if(sum((targets_all %in% unique(c(signaling_graph_list$gr$to))) == FALSE) > 0)
    stop("target genes should be in signaling_graph_list")

  if(!is.character(sig_color) | length(sig_color) != 1)
    stop("sig_color should be a character vector of length 1, denoting a color of interest")
  if(!is.character(gr_color) | length(gr_color) != 1)
    stop("gr_color should be a character vector of length 1, denoting a color of interest")

  requireNamespace("dplyr")

  tf_signaling = signaling_graph_list$sig
  tf_regulatory = signaling_graph_list$gr

  edge_list = bind_rows(tf_signaling %>% mutate(edge_color = sig_color),tf_regulatory %>% mutate(edge_color = gr_color))
  nodes = unique(c(edge_list$from,edge_list$to))
  nodes_ids = 1:length(nodes)
  names(nodes_ids) = nodes

  true_colors = rep("grey50",times = length(nodes))

  true_colors[which(nodes %in% ligands_all)] = sig_color
  true_colors[which(nodes %in% targets_all)] = gr_color
  # Create a node data frame
  nodes_ = DiagrammeR::create_node_df(n = length(nodes_ids),
                                      nodes = nodes_ids,
                                      label = names(nodes_ids),
                                      style = "filled",
                                      width = 0.75,
                                      fillcolor = true_colors,
                                      fontcolor = "white")

  edge_list_ = edge_list %>% mutate(from = nodes_ids[from], to = nodes_ids[to])


  edges_ = DiagrammeR::create_edge_df(
    from = edge_list_$from,
    to = edge_list_$to,
    penwidth = edge_list$weight,
    color = edge_list$edge_color
  )

  graph = DiagrammeR::create_graph(
    nodes_df = nodes_,
    edges_df = edges_
  )
  return(graph)
}
#' @title Get the data sources that support the specific interactions in the extracted ligand-target signaling subnetwork
#'
#' @description \code{infer_supporting_datasources} Get the data sources that support the specific interactions in the extracted ligand-target signaling subnetwork
#'
#' @usage
#' infer_supporting_datasources(signaling_graph_list,lr_network, sig_network, gr_network)
#'
#' @inheritParams construct_weighted_networks
#' @inheritParams diagrammer_format_signaling_graph
#'
#' @return A tibble with columns from, to, source and layer
#'
#' @examples
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' data_source_info_network = infer_supporting_datasources(ligand_target_signaling_list, lr_network, sig_network, gr_network)
#'
#' @export
#'
infer_supporting_datasources = function(signaling_graph_list,lr_network, sig_network , gr_network){

  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(signaling_graph_list))
    stop("signaling_graph_list must be a list object")
  if (!is.data.frame(signaling_graph_list$sig))
    stop("signaling_graph_list$sig must be a data frame or tibble object")
  if (!is.data.frame(signaling_graph_list$gr))
    stop("signaling_graph_list$gr must be a data frame or tibble object")


  requireNamespace("dplyr")

  tf_signaling = signaling_graph_list$sig
  tf_regulatory = signaling_graph_list$gr

  signaling_filtered = tf_signaling %>% dplyr::select(from,to) %>% distinct()
  regulatory_filtered = tf_regulatory %>% dplyr::select(from,to) %>% distinct()

  bind_rows(inner_join(regulatory_filtered, gr_network, by = c("from","to")) %>% mutate(layer = "regulatory"), inner_join(signaling_filtered,bind_rows(lr_network, sig_network), by = c("from","to")) %>% mutate(layer = "ligand_signaling"))
}


