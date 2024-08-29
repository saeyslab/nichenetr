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
#' @param minmax_scaling Indicate whether the weights of both dataframes should be min-max scaled between 0.75 and 1. Default: FALSE.
#'
#' @return A list containing 2 elements (sig and gr): the integrated weighted ligand-signaling and gene regulatory networks data frame / tibble format with columns: from, to, weight
#'
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' }
#' @export
#'
get_ligand_signaling_path = function(ligand_tf_matrix, ligands_all, targets_all, top_n_regulators = 4, weighted_networks, ligands_position = "cols", minmax_scaling = FALSE){

  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

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
  if(!is.numeric(top_n_regulators) | length(top_n_regulators) != 1 | top_n_regulators <= 0)
    stop("top_n_regulators should be a number higher than 0")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  requireNamespace("dplyr")

  final_combined_df  = construct_ligand_signaling_df(ligands_all,targets_all, top_n_regulators, weighted_networks,ligand_tf_matrix)

  signaling_network_all = weighted_networks$lr_sig %>% mutate(weight = 1/weight) # inverse weight to prepare for SPL
  signaling_igraph = igraph::graph_from_data_frame(signaling_network_all, directed = TRUE)

  tf_nodes = lapply(ligands_all,get_shortest_path_signaling, final_combined_df,signaling_igraph) %>% unlist() %>% unique()

  tf_signaling = weighted_networks$lr_sig %>%
    filter(from %in% c(ligands_all,tf_nodes) & to %in% tf_nodes) %>% group_by(from,to) %>% mutate(weight = sum(weight)) %>% ungroup() %>% distinct()
  tf_regulatory = weighted_networks$gr %>%
    filter(from %in% final_combined_df$TF & to %in% targets_all)  %>% ungroup() %>% distinct()

  if (minmax_scaling){
    tf_signaling <- tf_signaling %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
    tf_regulatory <- tf_regulatory %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  }

  return(list(sig = tf_signaling, gr = tf_regulatory))
}
#' @title Get ligand-target signaling paths between ligand(s), receptors, and target gene(s) of interest
#'
#' @description \code{get_ligand_signaling_path_with_receptor} Extract possible signaling paths between a ligand(s), receptor(s) and target gene(s) of interest. The most highly weighted path(s) will be extracted.
#'
#' @usage
#' get_ligand_signaling_path_with_receptor(ligand_tf_matrix, ligands_all, receptors_all, targets_all, top_n_regulators = 3, weighted_networks, ligands_position = "cols")
#'
#' @inheritParams get_ligand_signaling_path
#' @param receptors_all A character vector of one or more receptors of interest
#'
#' @return A list containing 2 elements (sig and gr): the integrated weighted ligand-signaling and gene regulatory networks data frame / tibble format with columns: from, to, weight
#'
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_receptors = c("BMPR2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path_with_receptor(ligand_tf_matrix,all_ligands,all_receptors, all_targets,top_n_regulators,weighted_networks)
#' }
#'
#' @export
#'
get_ligand_signaling_path_with_receptor = function (ligand_tf_matrix, ligands_all, receptors_all, targets_all, top_n_regulators = 3, weighted_networks, ligands_position = "cols") {
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")
  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")
  if (!is.matrix(ligand_tf_matrix))
    stop("ligand_tf_matrix should be a matrix")
  if (ligands_position == "cols") {
    if (sum((ligands_all %in% colnames(ligand_tf_matrix)) ==
            FALSE) > 0)
      stop("ligands should be in ligand_tf_matrix")
  }
  else if (ligands_position == "rows") {
    if (sum((ligands_all %in% rownames(ligand_tf_matrix)) ==
            FALSE) > 0)
      stop("ligands should be in ligand_tf_matrix")
  }
  if (sum((targets_all %in% unique(c(weighted_networks$gr$to))) ==
          FALSE) > 0)
    stop("target genes should be in gene regulatory network")
  if (!is.numeric(top_n_regulators) | length(top_n_regulators) !=
      1 | top_n_regulators <= 0)
    stop("top_n_regulators should be a number higher than 0")
  if (ligands_position != "cols" & ligands_position !=
      "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  requireNamespace("dplyr")
  final_combined_df = construct_ligand_signaling_df(ligands_all,
                                                    targets_all, top_n_regulators, weighted_networks, ligand_tf_matrix)
  signaling_network_all = weighted_networks$lr_sig %>% mutate(weight = 1/weight)
  signaling_igraph = igraph::graph_from_data_frame(signaling_network_all,
                                                   directed = TRUE)
  tf_nodes = lapply(ligands_all, get_shortest_path_signaling,
                    final_combined_df, signaling_igraph) %>% unlist() %>%
    unique()
  tf_signaling = weighted_networks$lr_sig %>% filter(from %in% c(ligands_all, receptors_all, tf_nodes) & to %in% tf_nodes) %>% group_by(from, to) %>% mutate(weight = sum(weight)) %>% ungroup() %>%
    distinct() %>% bind_rows(weighted_networks$lr_sig %>% filter(from %in% ligands_all & to %in% receptors_all) %>% group_by(from, to) %>% mutate(weight = sum(weight)) %>% ungroup() %>%
                               distinct()) %>%
    distinct()
  tf_regulatory = weighted_networks$gr %>% filter(from %in%
                                                    final_combined_df$TF & to %in% targets_all) %>% ungroup() %>%
    distinct()
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
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' graph = diagrammer_format_signaling_graph(ligand_target_signaling_list, all_ligands,all_targets)
#' # DiagrammeR::render_graph(graph, layout = "tree")
#' }
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
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' all_ligands = c("BMP2")
#' all_targets = c("HEY1")
#' top_n_regulators = 2
#' ligand_target_signaling_list = get_ligand_signaling_path(ligand_tf_matrix,all_ligands,all_targets,top_n_regulators,weighted_networks)
#' data_source_info_network = infer_supporting_datasources(ligand_target_signaling_list, lr_network, sig_network, gr_network)
#' }
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
#' @title Make a ggplot heatmap object from an input matrix (2-color).
#'
#' @description \code{make_heatmap_ggplot} Make a ggplot heatmap object from an input matrix containing continuous values. Two-color scale from white to color of choice.
#'
#' @usage
#' make_heatmap_ggplot(matrix, y_name, x_name, y_axis = TRUE,x_axis = TRUE, x_axis_position = "top", legend_position = "top", color = "blue", legend_title = "score",...)
#'
#' @param matrix Matrix with continuous values to plot in heatmap
#' @param y_name Title of the y-axis
#' @param x_name Title of the x-axis
#' @param y_axis Should y-axis label names and titles be displayed? TRUE or FALSE. Default: TRUE.
#' @param x_axis Should x-axis label names and titles be displayed? TRUE or FALSE. Default: TRUE.
#' @param x_axis_position X-axis position: "top" or "bottomm"; only relevant if x_axis == TRUE. Default:"top".
#' @param legend_position Legend position: "top", "bottom", "left", "right" or "none". Default: "top".
#' @param color Color for highest continuous value in heatmap. Color gradient will go from "whitesmoke" to this color. Default: "blue".
#' @param legend_title Title of the legend.
#' @param ... Optional arguments passed to element_text(); used to set font type and size of axis labels and titles.
#'
#' @return A ggplot object displaying a heatmap
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' p = make_heatmap_ggplot(ligand_target_matrix[1:50,] %>% t(), y_name = "ligand", x_name = "target")
#' }
#' @export
#'
make_heatmap_ggplot = function(matrix, y_name, x_name, y_axis = TRUE,x_axis = TRUE, x_axis_position = "top", legend_position = "top", color = "blue",legend_title = "score", ...){

  # input checks
  if(!is.matrix(matrix))
    stop("matrix should be a matrix")
  if(!is.character(y_name) | length(y_name) != 1)
    stop("y_name should be a character vector of length 1")
  if(!is.character(x_name) | length(x_name) != 1)
    stop("x_name should be a character vector of length 1")
  if(!is.logical(y_axis) | length(y_axis) != 1)
    stop("y_axis should be a TRUE or FALSE")
  if(!is.logical(x_axis) | length(x_axis) != 1)
    stop("x_axis should be a TRUE or FALSE")
  if((x_axis_position %in% c("top","bottom")) == FALSE)
    stop("x_axis_position should be top or bottom")
  if((legend_position %in% c("top","bottom","left","right","none")) == FALSE)
    stop("legend_position should be top, bottom, left, right or none")
  if(!is.character(color) |  length(color) != 1)
    stop("color should be character vector of length 1")


  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  matrix_df_vis = matrix %>% data.frame(check.names = FALSE) %>% rownames_to_column("y") %>% as_tibble() %>% gather(x,"score", -y) %>% mutate(y = factor(y, levels = rownames(matrix), ordered = TRUE), x = factor(x, levels = colnames(matrix), ordered = TRUE))

  plot_object = matrix_df_vis %>% ggplot(aes(x,y,fill = score)) + geom_tile(color = "white", size = 0.5) + scale_fill_gradient(low = "whitesmoke", high = color) + theme_minimal()

  if (x_axis == FALSE){
    if(y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object  + ylab(paste0(y_name))
    } else if (y_axis == FALSE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object
    }

  } else if (x_axis == TRUE) {
    if (y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name)) + ylab(paste0(y_name))
    } else if (y_axis == FALSE) {

      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name))
    }
  }
  plot_object = plot_object + labs(fill = legend_title)
}
#' @title Make a ggplot heatmap object from an input matrix (3-color).
#'
#' @description \code{make_threecolor_heatmap_ggplot} Make a ggplot heatmap object from an input matrix containing continuous values. Three-color scale with colors of choice. Ideal for plotting log fold change expression.
#'
#' @usage
#' make_threecolor_heatmap_ggplot(matrix, y_name, x_name, y_axis = TRUE,x_axis = TRUE, x_axis_position = "top", legend_position = "top", low_color = "blue",mid_color = "whitesmoke", high_color = "red",mid = 0,legend_title = "score",...)
#'
#' @param low_color Color for lowest continuous value in heatmap. Color gradient will go from "whitesmoke" to this color. Default: "blue".
#' @param mid_color Color for the "mid" value as defined by that parameter. Default: "whitesmoke".
#' @param high_color Color for highest continuous value in heatmap. Color gradient will go from "whitesmoke" to this color. Default: "red".
#' @param mid Continuous value that will receive the "mid_color" color. Default: 0
#' @inheritParams make_heatmap_ggplot
#'
#' @return A ggplot object displaying a heatmap
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' p = make_threecolor_heatmap_ggplot(ligand_target_matrix[1:50,] %>% t(), y_name = "ligand", x_name = "target")
#' }
#' @export
#'
make_threecolor_heatmap_ggplot = function(matrix, y_name, x_name, y_axis = TRUE,x_axis = TRUE, x_axis_position = "top", legend_position = "top", low_color = "blue",mid_color = "whitesmoke", high_color = "red",mid = 0, legend_title = "score",...){

  # input checks
  if(!is.matrix(matrix))
    stop("matrix should be a matrix")
  if(!is.character(y_name) | length(y_name) != 1)
    stop("y_name should be a character vector of length 1")
  if(!is.character(x_name) | length(x_name) != 1)
    stop("x_name should be a character vector of length 1")
  if(!is.logical(y_axis) | length(y_axis) != 1)
    stop("y_axis should be a TRUE or FALSE")
  if(!is.logical(x_axis) | length(x_axis) != 1)
    stop("x_axis should be a TRUE or FALSE")
  if((x_axis_position %in% c("top","bottom")) == FALSE)
    stop("x_axis_position should be top or bottom")
  if((legend_position %in% c("top","bottom","left","right","none")) == FALSE)
    stop("legend_position should be top, bottom, left, right or none")
  if(!is.character(low_color) | !is.character(mid_color) | !is.character(high_color) | length(low_color) != 1 | length(mid_color) != 1 | length(high_color) != 1)
    stop("low_color, mid_color and high_color should be character vectors of length 1")
  if(!is.numeric(mid) | length(mid) != 1)
    stop("mid should be a numeric vector of length 1")




  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  matrix_df_vis = matrix %>% data.frame(check.names = FALSE) %>% rownames_to_column("y") %>% as_tibble() %>% gather(x,"score", -y) %>% mutate(y = factor(y, levels = rownames(matrix), ordered = TRUE), x = factor(x, levels = colnames(matrix), ordered = TRUE))

  plot_object = matrix_df_vis %>% ggplot(aes(x,y,fill = score)) + geom_tile(color = "white", size = 0.5) + scale_fill_gradient2(low = low_color, mid = mid_color,high = high_color, midpoint = mid) + theme_minimal()

  if (x_axis == FALSE){
    if(y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object  + ylab(paste0(y_name))
    } else if (y_axis == FALSE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object
    }

  } else if (x_axis == TRUE) {
    if (y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name)) + ylab(paste0(y_name))
    } else if (y_axis == FALSE) {

      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name))
    }
  }
  plot_object = plot_object + labs(fill = legend_title)
}
#' @title Make a ggplot heatmap object from an input ligand-target matrix.
#'
#' @description \code{make_heatmap_bidir_lt_ggplot} Make a ggplot heatmap object from an input ligand-target matrix in which it is indicated whether a gene is a top target of a ligand ("top-target"), the ligand is a top ligand of the gene ("top-ligand") or both ("top") or none ("none").
#'
#' @usage
#' make_heatmap_bidir_lt_ggplot(matrix, y_name, x_name, y_axis = TRUE, x_axis = TRUE, x_axis_position = "top", legend_position = "top", ...)
#' #'
#' @param matrix Matrix with continuous values to plot in heatmap
#' @param y_name Title of the y-axis
#' @param x_name Title of the x-axis
#' @param y_axis Should y-axis label names and titles be displayed? TRUE or FALSE. Default: TRUE.
#' @param x_axis Should x-axis label names and titles be displayed? TRUE or FALSE. Default: TRUE.
#' @param x_axis_position X-axis position: "top" or "bottomm"; only relevant if x_axis == TRUE. Default:"top".
#' @param legend_position Legend position: "top", "bottom", "left", "right" or "none". Default: "top"
#' @param ... Optional arguments passed to element_text(); used to set font type and size of axis labels and titles.
#'
#' @return A ggplot object displaying a heatmap
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#'
#' ligand_target_matrix_vis_genedirection = ligand_target_matrix %>% apply(1,scaling_modified_zscore) %>% .[,1:50]
#' ligand_target_matrix_vis_genedirection[ligand_target_matrix_vis_genedirection < 2] = 0
#' ligand_target_matrix_vis_genedirection[ligand_target_matrix_vis_genedirection != 0] = 1
#'
#' ligand_target_matrix_vis_liganddirection = ligand_target_matrix %>% apply(2,scaling_modified_zscore) %>% .[1:50, ] %>% t()
#' ligand_target_matrix_vis_liganddirection[ligand_target_matrix_vis_liganddirection < 2] = 0
#' ligand_target_matrix_vis_liganddirection[ligand_target_matrix_vis_liganddirection != 0] = 2
#'
#' bidirectional_ligand_target_matrix_vis = ligand_target_matrix_vis_genedirection + ligand_target_matrix_vis_liganddirection
#' bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 0] = "none"
#' bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 1] = "top-ligand"
#' bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 2] = "top-target"
#' bidirectional_ligand_target_matrix_vis[bidirectional_ligand_target_matrix_vis == 3] = "top"

#' p = make_heatmap_bidir_lt_ggplot(bidirectional_ligand_target_matrix_vis, y_name = "ligand", x_name = "target")
#' }
#' @export
#'
make_heatmap_bidir_lt_ggplot = function(matrix, y_name, x_name, y_axis = TRUE, x_axis = TRUE, x_axis_position = "top", legend_position = "top", ...){

  # input checks
  if(!is.matrix(matrix))
    stop("matrix should be a matrix")
  if(!is.character(y_name) | length(y_name) != 1)
    stop("y_name should be a character vector of length 1")
  if(!is.character(x_name) | length(x_name) != 1)
    stop("x_name should be a character vector of length 1")
  if(!is.logical(y_axis) | length(y_axis) != 1)
    stop("y_axis should be a TRUE or FALSE")
  if(!is.logical(x_axis) | length(x_axis) != 1)
    stop("x_axis should be a TRUE or FALSE")
  if((x_axis_position %in% c("top","bottom")) == FALSE)
    stop("x_axis_position should be top or bottom")
  if((legend_position %in% c("top","bottom","left","right","none")) == FALSE)
    stop("legend_position should be top, bottom, left, right or none")


  requireNamespace("dplyr")
  requireNamespace("ggplot2")

  matrix_df_vis = matrix %>% data.frame(stringsAsFactors = FALSE) %>% rownames_to_column("y") %>% as_tibble() %>% gather(x,"score", -y) %>% mutate(y = factor(y, levels = rownames(matrix) %>% make.names(), ordered = TRUE), x = factor(x, levels = colnames(matrix) %>% make.names(), ordered = TRUE))

  plot_object = matrix_df_vis %>% ggplot(aes(x,y,fill = score)) + geom_tile(color = "white", size = 0.5) +
    scale_fill_manual(values = c("top-ligand" = "indianred1", "top-target" = "lightskyblue1", "top" = "mediumpurple2", "none" = "whitesmoke")) + theme_minimal()

  if (x_axis == FALSE){
    if(y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object  + ylab(paste0(y_name))
    } else if (y_axis == FALSE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x =  element_blank(),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object
    }

  } else if (x_axis == TRUE) {
    if (y_axis == TRUE){
      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title = element_text(...), axis.text.y = element_text(...))
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name)) + ylab(paste0(y_name))
    } else if (y_axis == FALSE) {

      plot_object = plot_object + theme(panel.grid.minor = element_line(color = "transparent"), panel.grid.major = element_line(color = "transparent"), legend.position = legend_position, axis.ticks = element_line(size = 0), axis.text.x.top = element_text(angle = 90, hjust = 0,...), axis.text.x = element_text(angle = 90, hjust =1,...),  axis.title.y = element_blank(), axis.text.y = element_blank())
      plot_object = plot_object + scale_x_discrete(position = x_axis_position) + xlab(paste0(x_name))
    }
  }
}

#' @title Make a "mushroom plot" of ligand-receptor interactions
#'
#' @description \code{make_mushroom_plot} Make a plot in which each glyph consists of two semicircles corresponding to ligand- and receptor- information. The size of the semicircle is the percentage of cells that express the protein, while the saturation corresponds to the scaled average expression value.
#'
#' @param prioritization_table A prioritization table as generated by \code{\link{generate_prioritization_tables}}.
#' @param top_n An integer indicating how many ligand-receptor pairs to show
#' @param show_rankings A logical indicating whether to show the ranking of the ligand-receptor pairs (default: FALSE)
#' @param show_all_datapoints A logical indicating whether to show all ligand-receptor pairs (default: FALSE, if true they will be grayed out)
#' @param true_color_range A logical indicating whether to use the default color range as determined by ggplot (TRUE, default) or set the limits to a range of 0-1 (FALSE)
#' @param use_absolute_rank A logical indicating to whether use the absolute prioritization rank to filter the top_n ligand-receptor pairs (default: FALSE)
#' @param size A string indicating which column to use for the size of the semicircles (default: "scaled_avg_exprs"; use column name without "_ligand" or "_receptor" suffix)
#' @param color A string indicating which column to use for the color of the semicircles (default: "scaled_p_val_adapted"; use column name without "_ligand" or "_receptor" suffix)
#' @param ligand_fill_colors A vector of the low and high colors to use for the ligand semicircle fill gradient (default: c("#DEEBF7", "#08306B"))
#' @param receptor_fill_colors A vector of the low and high colors to use for the receptor semicircle fill gradient (default: c("#FEE0D2", "#A50F15"))
#' @param unranked_ligand_fill_colors A vector of the low and high colors to use for the unranked ligands when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))
#' @param unranked_receptor_fill_colors A vector of the low and high colors to use for the unranked receptors when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))
#' @param ... Additional arguments passed to \code{\link{ggplot2::theme}}. As there are often issues with the scales legend, it is recommended to change legend sizes and positions using this argument, i.e., \code{legend.key.height}, \code{legend.key.width}, \code{legend.title}, and \code{legend.text}.
#'
#' @details
#' If the values range of the column used as the "size" parameter is not between 0 and 1.001, an error will be thrown.
#'
#' The sender cell types can be ordered by encoding the "sender" column as a factor. If the "sender" column is not a factor, the sender cell types will be ordered alphabetically.
#'
#' By default, the top_n ligand-receptor pairs are shown despite their absolute ranking. So, if a receiver cell type has LR pairs that are only ranked from 31-40 and the top_n is set to 20, the LR pairs will be shown. If use_absolute_rank is set to TRUE, only LR pairs with absolute ranking from 1-20 will be shown.
#'
#'
#' @return A ggplot object
#'
#' @import ggplot2
#' @import ggforce
#' @import ggnewscale
#' @import shadowtext
#' @import cowplot
#'
#' @examples
#' \dontrun{
#' # Create a prioritization table
#' prior_table <- generate_prioritization_tables(processed_expr_table, processed_DE_table, ligand_activities, processed_condition_markers, prioritizing_weights)
#' make_mushroom_plot(prior_table)
#'
#'
#' # Show only top 20, and write rankings on the plot
#' make_mushroom_plot(prior_table, top_n = 20, show_rankings = TRUE)
#'
#' # Show all datapoints, and use true color range
#' make_mushroom_plot(prior_table, show_all_datapoints = TRUE, true_color_range = TRUE)
#'
#'
#' # Change the size and color columns
#' make_mushroom_plot(prior_table, size = "pct_expressed", color = "scaled_avg_exprs")
#'
#'
#' # For a prioritization table with multiple receiver cell types
#' make_mushroom_plot(prior_table_combined %>% filter(receiver == celltype_oi))
#'}
#'
#' @export
#'
make_mushroom_plot <- function(prioritization_table, top_n = 30, show_rankings = FALSE,
                              show_all_datapoints = FALSE, true_color_range = TRUE,
                              use_absolute_rank = FALSE,
                              size = "scaled_avg_exprs", color = "scaled_p_val_adapted",
                              ligand_fill_colors = c("#DEEBF7", "#08306B"),
                              receptor_fill_colors = c("#FEE0D2", "#A50F15"),
                              unranked_ligand_fill_colors = c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)),
                              unranked_receptor_fill_colors = c( alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)),
                              ...){
  size_ext <-  c("ligand", "receptor"); color_ext <- c("ligand", "receptor")
  if (size == "pct_expressed") size_ext <- c("sender", "receiver")
  if (color == "pct_expressed") color_ext <- c("sender", "receiver")

  cols_to_use <- c("sender", "ligand", "receptor", paste0(size, "_", size_ext), paste0(color, "_", color_ext))

  if (!all(cols_to_use %in% colnames(prioritization_table))){
    stop(paste(paste0("`", cols_to_use %>% .[!. %in% colnames(prioritization_table)], "`", collapse =", "), "column not in prioritization table"))
  }
  if(!is.logical(show_rankings) | length(show_rankings) != 1)
    stop("show_rankings should be a TRUE or FALSE")
  if(!is.logical(show_all_datapoints) | length(show_all_datapoints) != 1)
    stop("show_all_datapoints should be a TRUE or FALSE")
  if(!is.logical(true_color_range) | length(true_color_range) != 1)
    stop("true_color_range should be a TRUE or FALSE")
  if(!is.logical(use_absolute_rank) | length(use_absolute_rank) != 1)
    stop("use_absolute_rank should be a TRUE or FALSE")
  if(!is.numeric(top_n) | length(top_n) != 1)
    stop("top_n should be a numeric vector of length 1")
  if(length(ligand_fill_colors) != 2)
    stop("ligand_fill_colors should be a vector of length 2")
  if(length(receptor_fill_colors) != 2)
    stop("receptor_fill_colors should be a vector of length 2")
  if(length(unranked_ligand_fill_colors) != 2)
    stop("unranked_ligand_fill_colors should be a vector of length 2")
  if(length(unranked_receptor_fill_colors) != 2)
    stop("unranked_receptor_fill_colors should be a vector of length 2")

  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  requireNamespace("ggnewscale")
  requireNamespace("ggforce")
  requireNamespace("shadowtext")
  requireNamespace("cowplot")

  if (!"prioritization_rank" %in% colnames(prioritization_table)){
    prioritization_table <- prioritization_table %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score)))
  }
  # Add 'relative rank' column which is basically 1:n
  prioritization_table <- prioritization_table %>% dplyr::mutate(relative_rank = rank(desc(prioritization_score)))

  # If use_absolute_rank, use 'prioritization_rank' column to filter top_n
  rank_filter_col <- ifelse(use_absolute_rank, "prioritization_rank", "relative_rank")

  # Create a new column of ligand-receptor interactions, and filter table to
  # only include LR interactions that appear in the top_n
  filtered_table <- prioritization_table %>% dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = " - "))
  order_interactions <- unique(filtered_table %>% filter(.data[[rank_filter_col]] <= top_n) %>% pull(lr_interaction))
  filtered_table <- filtered_table %>% filter(lr_interaction %in% order_interactions) %>%
    mutate(lr_interaction = factor(lr_interaction, levels = rev(order_interactions)))

  # Check if filtered_table is empty
  if (nrow(filtered_table) == 0){
    stop("No ligand-receptor interactions found in the top_n. Please try use_absolute_rank = FALSE or increase top_n.")
  }

  # Keep order of senders, if present (if not, sort alphabetically)
  if (!is.factor(filtered_table$sender)){
    filtered_table$sender <- as.factor(filtered_table$sender)
  } else {
    # Drop levels that are not present in the filtered table
    filtered_table$sender <- droplevels(filtered_table$sender)
  }

  lr_interaction_vec <- 1:length(order_interactions) %>% setNames(order_interactions)

  # Make each ligand and receptor into separate rows (to draw 1 semicircle per row)
  filtered_table <- filtered_table %>% select(c("lr_interaction", all_of(cols_to_use), "prioritization_rank", "relative_rank")) %>%
    pivot_longer(c(ligand, receptor), names_to = "type", values_to = "protein") %>%
    mutate(size = ifelse(type == "ligand", get(paste0(size, "_", size_ext[1])), get(paste0(size, "_", size_ext[2]))),
           color = ifelse(type == "ligand", get(paste0(color, "_", color_ext[1])), get(paste0(color, "_",  color_ext[2])))) %>%
    select(-contains(c("_ligand", "_receptor", "_sender", "_receiver"))) %>%
    mutate(start = rep(c(-pi, 0), nrow(filtered_table))) %>%
    mutate(x = as.numeric(sender), y = lr_interaction_vec[lr_interaction])

  # Warning if size column is not scaled between 0 and 1.001
  if (any(filtered_table$size < 0) | any(filtered_table$size > 1.001)){
    stop("Size column is not scaled between 0 and 1. Please use this column as the color instead.")
  }

  # Rename size and color columns to be more human-readable
  keywords_adj <- c("LFC", "pval", "", "product", "mean", "adjusted", "expression") %>% setNames(c("lfc", "p", "val", "prod", "avg", "adj", "exprs"))
  size_title <- sapply(stringr::str_split(size, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>%  stringr::str_replace("^\\w{1}", toupper)
  color_title <- sapply(stringr::str_split(color, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>% stringr::str_replace("^\\w{1}", toupper)

  color_lims <- c(0,1)
  if (true_color_range) color_lims <- NULL

  scale <- 0.5

  ncelltypes <- length(unique(filtered_table$sender))
  n_interactions <- length(lr_interaction_vec)
  legend2_df <- data.frame(values = c(0.25, 0.5, 0.75, 1), x=(ncelltypes+2.5):(ncelltypes+5.5), y=rep(floor(n_interactions/3), 4), start=-pi)
  axis_rect <- data.frame(xmin=0, xmax=ncelltypes+1, ymin=0, ymax=n_interactions+1)
  panel_grid_y <- data.frame(x = rep(seq(from = 0.5, to = ncelltypes+0.5, by = 1), each=2),
                             y = c(n_interactions+1, 0), group = rep(1:(ncelltypes+1), each=2))
  panel_grid_x <- data.frame(y = rep(seq(from = 0.5, to = n_interactions+0.5, by = 1), each=2),
                             x = c(ncelltypes+1, 0), group = rep(1:(n_interactions+1), each=2))

  theme_args <- list(panel.grid.major = element_blank(),
                     legend.box = "horizontal",
                     panel.background = element_blank())

  theme_args[names(list(...))] <- list(...)

  # Check if legend.title is in the extra arguments
    # Multiply by ratio of 5/14: see https://stackoverflow.com/questions/25061822/ggplot-geom-text-font-size-control
  scale_legend_title_size <- ifelse("legend.title" %in% names(theme_args), theme_args$legend.title$size*(5/14), GeomLabel$default_aes$size)
  # Check if legend.text is in the extra arguments
  scale_legend_text_size <- ifelse("legend.text" %in% names(theme_args), theme_args$legend.text$size*(5/14), GeomLabel$default_aes$size)

  # Check if legend.justification is in the extra arguments
  if (!"legend.justification" %in% names(theme_args)) {
    theme_args$legend.justification <- c(1, 0.7)
  }

  # Check if legend.position is in the extra arguments
  if (!"legend.position" %in% names(theme_args)){
    theme_args$legend.position <- c(1, 0.7)
  }

  p1 <- ggplot() +
    # Draw ligand semicircle
    geom_arc_bar(data = filtered_table %>% filter(type=="ligand",  .data[[rank_filter_col]] <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    scale_fill_gradient(low = ligand_fill_colors[1] , high=ligand_fill_colors[2] ,
                        limits=color_lims, oob=scales::squish,
                        n.breaks = 3,
                        guide = guide_colorbar(order = 1),
                        name=paste0(color_title, " (", color_ext[1], ")") %>% stringr::str_wrap(width=15)) +
    # Create new fill scale for receptor semicircles
    new_scale_fill() +
    geom_arc_bar(data = filtered_table %>% filter(type=="receptor", .data[[rank_filter_col]] <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    # Size legend
    geom_arc_bar(data = legend2_df, aes(x0=x, y0=y, r0=0, r=sqrt(values)*scale, start=start, end=start+pi), fill="black") +
    geom_rect(data = legend2_df, aes(xmin=x-0.5, xmax=x+0.5, ymin=y-0.5, ymax=y+0.5), color="gray90", fill=NA) +
    geom_text(data = legend2_df, aes(label=values, x=x, y=y-0.6), vjust=1, size = scale_legend_text_size) +
    geom_text(data = data.frame(x = (ncelltypes+4), y = floor(n_interactions/3)+1,
                                label = size_title %>% stringr::str_wrap(width=15)),
              aes(x=x, y=y, label=label), size = scale_legend_title_size, vjust=0, lineheight = .75) +
    # Panel grid
    geom_line(data = panel_grid_y, aes(x=x, y=y, group=group), color = "gray90") +
    geom_line(data = panel_grid_x, aes(x=x, y=y, group=group), color = "gray90") +
    # Draw box to represent x and y "axis"
    geom_rect(data = axis_rect, aes(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), color = "black", fill = "transparent") +
    # Other plot information
    scale_fill_gradient(low = receptor_fill_colors[1], high=receptor_fill_colors[2] , limits=color_lims, oob=scales::squish, n.breaks = 3,
                        name=paste0(color_title,  " (", color_ext[2], ")") %>% stringr::str_wrap(width=15),
                        guide = guide_colorbar(order = 2)) +
    scale_y_continuous(breaks=n_interactions:1, labels=names(lr_interaction_vec), expand = expansion(add=c(0,0))) +
    scale_x_continuous(breaks=1:ncelltypes, labels=levels(filtered_table$sender), position="top", expand = expansion(add=c(0,0))) +
    xlab("Sender cell types") + ylab("Ligand-receptor interaction") +
    coord_fixed() +
    do.call(theme, theme_args)

  # Add unranked ligand and receptor semicircles if requested
  if (show_all_datapoints){

    # Limits will depend on true_color_range
    unranked_ligand_lims <- c(0,1); unranked_receptor_lims <- c(0,1)
    if (true_color_range){
      # Follow limits of the top_n lr pairs
      unranked_ligand_lims <- filtered_table %>% filter(type=="ligand",  .data[[rank_filter_col]] <= top_n) %>%
        select(color) %>% range
      unranked_receptor_lims <- filtered_table %>% filter(type=="receptor",  .data[[rank_filter_col]] <= top_n) %>%
        select(color) %>% range
    }

    p1 <- p1 + new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="ligand", .data[[rank_filter_col]] > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low = unranked_ligand_fill_colors[1], high=unranked_ligand_fill_colors[2],
                          limits=unranked_ligand_lims, oob = scales::oob_squish,
                          guide = "none") +
      new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="receptor", .data[[rank_filter_col]] > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low=unranked_receptor_fill_colors[1], high=unranked_receptor_fill_colors[2],
                          limits=unranked_receptor_lims, oob = scales::oob_squish,
                          guide = "none")
  }

  # Add ranking numbers if requested
  if (show_rankings){
    p1 <- p1 + geom_shadowtext(data = filtered_table %>% filter(.data[[rank_filter_col]] <= top_n),
                               aes(x=x, y=y, label=prioritization_rank))
  }

  p1
}


## Circos plot functions
#' @title Assign ligands to cell types
#' @description Assign ligands to a sender cell type, based on the strongest expressing cell type of that ligand. Ligands are only assigned to a cell type if that cell type is the only one to show an expression that is higher than the average + SD. Otherwise, it is assigned to "General".
#' @param seuratObj Seurat object
#' @param ligands Vector of ligands to assign to cell types
#' @param celltype_col Metadata column name in the Seurat object that contains the cell type information
#' @param func.agg Function to use to aggregate the expression of a ligand across all cells in a cell type (default = mean)
#' @param func.assign Function to use to assign a ligand to a cell type (default = mean + SD)
#' @param condition_oi Condition of interest to subset the Seurat object (default = NULL)
#' @param condition_col Metadata column name in the Seurat object that contains the condition of interest (default = NULL)
#' @param ... Arguments passed to Seurat::GetAssayData, e.g., for the slot/layer to use (default: data)
#' @return A data frame of two columns, the cell type the ligand has been assigned to (\code{ligand_type}) and the ligand name (\code{ligand})
#' @details If the provided slot/layer is "data",  the normalized counts are first exponentiated before aggregation is performed
#' @export
#' @examples \dontrun{
#' assign_ligands_to_celltype(seuratObj = seuratObj, ligands = best_upstream_ligands[1:20],
#'                           celltype_col = "celltype", func.agg = mean, func.assign = function(x) {mean(x)+sd(x)},
#'                           condition_oi = "LCMV", condition_col = "aggregate", slot = "data")
#' }
#'
assign_ligands_to_celltype <- function(seuratObj, ligands, celltype_col, func.agg = mean, func.assign = function(x) {mean(x)+sd(x)},
                                        condition_oi = NULL, condition_col = NULL, ...) {
  # Check that if condition_oi is given, then so is condition_oi, and vice versa
  if (any(!is.na(condition_col), !is.na(condition_oi)) & !all(!is.na(condition_col), !is.na(condition_oi))){
    stop("Please input both condition_colname and condition_oi")
  }

  # Check that all ligands are in the seurat object
  if (any(!ligands %in% rownames(seuratObj))){
    stop("Not all ligands are in the Seurat object")
  }

  slot <- "data"
  # Check if slot or layer is provided
  if (length(list(...)) > 0) {
    if (any(grepl("slot|layer", names(list(...))))){
      slot <- list(...)[[which(grepl("slot|layer", names(list(...))))]]
    } else {
      warning("No slot/layer provided even though extra argument was provided, using default slot = 'data'")
    }
  }

  seuratObj_subset <- subset(seuratObj, features = ligands)

  # Calculate average ligand expression in sender cells
  if (!is.null(condition_oi)){
    seuratObj_subset <- seuratObj_subset[, seuratObj_subset[[condition_col]] == condition_oi ]
  }
  avg_expression_ligands <- lapply(unique(seuratObj_subset[[celltype_col, drop=TRUE]]), function (celltype) {
    if (slot == "data"){
      # Exponentiate-1 and calculate in non-log space
      expm1(GetAssayData(seuratObj_subset[, seuratObj_subset[[celltype_col]] == celltype], ...)) %>%
        apply(1, func.agg)

    } else {
      apply(GetAssayData(seuratObj_subset[, seuratObj_subset[[celltype_col]] == celltype], ...), 1, func.agg)

    }
    }) %>% setNames(unique(seuratObj_subset[[celltype_col, drop=TRUE]])) %>%
    do.call(cbind, .)

  sender_ligand_assignment <- avg_expression_ligands %>% apply(1, function(ligand_expression){
    ligand_expression > func.assign(ligand_expression)
  }) %>% t()
  sender_ligand_assignment <- sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})

  all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
  unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
  general_ligands = ligands %>% setdiff(unique_ligands)

  ligand_type_indication_df <- lapply(names(sender_ligand_assignment), function(sender) {
    unique_ligands_sender <- names(sender_ligand_assignment[[sender]]) %>% setdiff(general_ligands)

    if (length(unique_ligands_sender) > 0) {
      return(data.frame(ligand_type = sender, ligand = unique_ligands_sender))
    }

  }) %>% bind_rows()

  ligand_type_indication_df <- bind_rows(ligand_type_indication_df,
                                         data.frame(ligand = general_ligands) %>% mutate(ligand_type = "General"))

  return(ligand_type_indication_df)
}

#' @title Get ligand-target links of interest
#' @usage get_ligand_target_links_oi(ligand_type_indication_df, active_ligand_target_links_df, cutoff = 0.40)
#' @description Filter ligand-target links based on a cutoff
#' @param ligand_type_indication_df Dataframe with column names \code{ligand_type} and \code{ligand}, from the function \code{\link{assign_ligands_to_celltype}}
#' @param active_ligand_target_links_df Dataframe with weighted ligand-target links from the function \code{\link{get_ligand_target_links}}, and an additional column \code{target_type} that indicates the grouping of target genes
#' @param cutoff Quantile to filter ligand-target links (default = 0.40, meaning 40\% of the lowest weighted ligand-target links are removed)
#' @return A dataframe with ligand-target links with weights above a certain cutoff. This dataframe also contains the attribute \code{cutoff_include_all_ligands}, which is the cutoff value of regulatory potential used at \code{cutoff} quantile.
#' @export
#' @examples \dontrun{
#' active_ligand_target_links_df <- lapply(best_upstream_ligands, get_weighted_ligand_target_links,
#'                                          geneset = geneset_oi,
#'                                          ligand_target_matrix = ligand_target_matrix,
#'                                          n = 200)
#' active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df))
#' ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj = seuratObj, ligands = best_upstream_ligands[1:20])
#' circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
#'                                            active_ligand_target_links_df %>% mutate(target_type = "LCMV-DE"),
#'                                            cutoff = 0.40)
#' attr(circos_links, "cutoff_include_all_ligands") # This is the cutoff value of regulatory potential used
#' }
#'
get_ligand_target_links_oi <- function(ligand_type_indication_df, active_ligand_target_links_df, cutoff = 0.40){
  # Check that ligand_type_indication_df has the correct colnames
  if (!all(c("ligand_type", "ligand") %in% colnames(ligand_type_indication_df))){
    stop("ligand_type_indication_df must have columns ligand_type and ligand")
  }

  # Check that active_ligand_target_links_df has the correct colnames
  if (!all(c("ligand", "target", "weight", "target_type") %in% colnames(active_ligand_target_links_df))){
    stop("active_ligand_target_links_df must have columns ligand, target, weight, and target_type")
  }

  # Check that cutoff is between 0 and 1
  if (cutoff < 0 | cutoff > 1){
    stop("cutoff must be between 0 and 1")
  }

  active_ligand_target_links_df <- active_ligand_target_links_df %>% inner_join(ligand_type_indication_df)
  cutoff_include_all_ligands <- active_ligand_target_links_df$weight %>% quantile(cutoff)
  active_ligand_target_links_df_circos <- active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
  ligands_to_remove <- setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
  targets_to_remove <- setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
  circos_links <- active_ligand_target_links_df %>% filter(!target %in% targets_to_remove & !ligand %in% ligands_to_remove)

  # Add this as an attribute
  attr(circos_links, "cutoff_include_all_ligands") <- cutoff_include_all_ligands
  return(circos_links)
}

#' @title Prepare circos visualization
#' @usage prepare_circos_visualization(circos_links, ligand_colors = NULL, target_colors = NULL, widths = NULL, celltype_order = NULL)
#' @description Prepare the data for the circos visualization by incorporating the colors and order of the links, as well as gaps between different cell types
#' @param circos_links Dataframe from the function \code{\link{get_ligand_target_links_oi}} containing weighted ligand-target links, cell type expressing the ligand, and target gene goruping
#' @param ligand_colors Named vector of colors for each cell type (default: NULL, where colors follow the ggplot default color scheme)
#' @param target_colors Named vector of colors for each target gene grouping (default: NULL, where colors follow the ggplot default color scheme)
#' @param widths Named list of widths for the different types groupings, including:
#' \itemize{
#' \item width_same_cell_same_ligand_type: Width of the links between ligands of the same cell type (default: 0.5)
#' \item width_different_cell: Width of the links between different cell types, or between different target gene groups (default: 6)
#' \item width_ligand_target: Width of the links between ligands and targets (default: 15)
#' \item width_same_cell_same_target_type: Width of the links between target genes of the same group (default: 0.5)
#' }
#' @param celltype_order Order of the cell types (default: NULL, where cell types are ordered alphabetically, followed by "General"). Cell types are then drawn counter-clockwise in the circos plot.
#' @return A list of four objects, including:
#' \itemize{
#' \item circos_links: Dataframe of weighted ligand-target links
#' \item ligand_colors: Named vector of ligands and their colors
#' \item order: Vector of order of the ligands and target genes
#' \item gaps: Vector of gaps between the different groupings
#' }
#' @examples \dontrun{
#' celltype_order <- c("General", "NK", "B", "DC", "Mono")
#' ligand_colors <- c("General" = "lawngreen", "NK" = "royalblue", "B" = "darkgreen", "Mono" = "violet", "DC" = "steelblue2")
#' target_colors <- c("LCMV-DE" = "tomato")
#' vis_circos_obj <- prepare_circos_visualization(circos_links, ligand_colors, target_colors, celltype_order = celltype_order)
#' }
#'
#' @export
prepare_circos_visualization <- function(circos_links, ligand_colors = NULL, target_colors = NULL, widths = NULL, celltype_order = NULL) {
  # Check that circos_links has the correct colnames
  if (!all(c("ligand", "target", "weight", "target_type", "ligand_type") %in% colnames(circos_links))){
    stop("circos_links must have columns ligand, target, weight, target_type, and ligand_type")
  }

  # If ligand_colors and/or target_colors is NULL, set to default ggplot colors (equally spaced colors around the color wheel)
  if (is.null(ligand_colors) | is.null(target_colors)){
    n_ligands <- is.null(ligand_colors)*length(unique(circos_links$ligand_type))
    n_targets <- is.null(target_colors)*length(unique(circos_links$target_type))
    n_total <- n_ligands + n_targets
    hues <- seq(15, 375, length = n_total + 1)

    if (is.null(target_colors)){
      target_colors <- hcl(h = hues, l = 65, c = 100)[(n_ligands+1):n_total]
      names(target_colors) <- unique(circos_links$target_type)
    }

    if (is.null(ligand_colors)){
      ligand_colors <- hcl(h = hues, l = 65, c = 100)[1:n_ligands]
      names(ligand_colors) <- unique(circos_links$ligand_type)
    }
  }

  # Check that ligand colors contains all ligand types
  if (!all(unique(circos_links$ligand_type) %in% names(ligand_colors))){
    stop("ligand_colors must contain all cell types in circos_links$ligand_type")
  }

    # If ligand colors contain additional cell types, intersect
  if (length(setdiff(names(ligand_colors), unique(circos_links$ligand_type))) > 0){
    warning("ligand_colors contains additional cell types not in circos_links$ligand_type, these will be removed")
    ligand_colors <- ligand_colors %>% .[names(.) %in% unique(circos_links$ligand_type)]
  }

  # Check that target colors contains all target types
  if (!all(unique(circos_links$target_type) %in% names(target_colors))){
    stop("target_colors must contain all target groupings in circos_links$target_type")
  }

  # If target colors contain additional target types, intersect
  if (length(setdiff(names(target_colors), unique(circos_links$target_type))) > 0){
    warning("target_colors contains additional target types not in circos_links$target_type, these will be removed")
    target_colors <- target_colors %>% .[names(.) %in% unique(circos_links$target_type)]
  }

  # Check that celltype_order contains all cell types
  if (!is.null(celltype_order) & !all(unique(circos_links$ligand_type) %in% celltype_order)){
    stop("celltype_order must contain all cell types in circos_links$ligand_type")
  }

  # If celltype_order contains additional cell types, intersect
  if (!is.null(celltype_order) & length(setdiff(celltype_order,  unique(circos_links$ligand_type))) > 0){
    warning("celltype_order contains additional cell types not in circos_links$ligand_type, these will be removed")
    celltype_order <- celltype_order %>% .[. %in% unique(circos_links$ligand_type)]
  }

  # If width is null, set default widths
  if (is.null(widths)){
    widths <- list(width_same_cell_same_ligand_type = 0.5,
                  width_different_cell = 6,
                  width_ligand_target = 15,
                  width_same_cell_same_target_type = 0.5)
  }

  # Check that widths contains all widths
  if (!all(c("width_same_cell_same_ligand_type", "width_different_cell", "width_ligand_target", "width_same_cell_same_target_type") %in% names(widths))){
    stop("widths must contain all four width names")
  }

  # Check that all widths are numeric
  if (!all(is.numeric(unlist(widths)))){
    stop("all widths must be numeric")
  }


  #  give each segment of ligands and targets a specific color and order
  grid_col_tbl_ligand <- tibble(ligand_type = ligand_colors %>% names(), color_ligand_type = ligand_colors)
  grid_col_tbl_target <- tibble(target_type = target_colors %>% names(), color_target_type = target_colors)

  circos_links <- circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
  circos_links <- circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
  links_circle <- circos_links %>% select(ligand, target, weight)

  ligand_color <- circos_links %>% distinct(ligand,color_ligand_type)
  grid_ligand_color <- ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
  target_color <- circos_links %>% distinct(target,color_target_type)
  grid_target_color <- target_color$color_target_type %>% set_names(target_color$target)

  grid_col <- c(grid_ligand_color, grid_target_color)

  # Prepare the circos visualization: order ligands and targets
  target_type_order <- circos_links %>% arrange(target_type) %>% pull(target_type) %>% unique()
  target_order <- circos_links %>% arrange(target_type, target) %>% pull(target) %>% unique()

  if (is.null(celltype_order)){
    circos_links_arranged <- circos_links %>% mutate(
      ligand_type_order = case_when(
        ligand_type == "General" ~ 1,
        TRUE ~ 2
      )) %>% arrange(ligand_type_order, desc(ligand_type)) %>%
      select(-ligand_type_order)
    ligand_type_order <- circos_links_arranged %>% pull(ligand_type) %>% unique
    ligand_order <- circos_links_arranged %>% pull(ligand) %>% unique
  } else {
    ligand_type_order <- celltype_order
    # Arrange circos_links according to celltype_order
    ligand_order <- circos_links %>% arrange(factor(ligand_type, levels = celltype_order), ligand) %>% pull(ligand) %>% unique

  }

  order <- c(ligand_order,target_order)

  # Prepare the circos visualization: define the gaps between the different segments
  gaps_sender_cell_types <- unlist(lapply(seq_along(ligand_type_order), function(i) {
    c(rep(widths$width_same_cell_same_ligand_type,
        times = (circos_links %>% filter(ligand_type == ligand_type_order[i]) %>% distinct(ligand) %>% nrow()-1)),
      if (i < length(ligand_type_order)) widths$width_different_cell)
    }))

  gaps_target_types <- unlist(lapply(seq_along(target_type_order), function(i) {
    c(rep(widths$width_same_cell_same_target_type,
          times = (circos_links %>% filter(target_type == target_type_order[i]) %>% distinct(target) %>% nrow()-1)),
      if (i < length(target_type_order)) widths$width_different_cell)
  }))

  gaps <- c(
    gaps_sender_cell_types,
    widths$width_ligand_target,
    gaps_target_types,
    widths$width_ligand_target
  )

  return(list(links_circle = links_circle, ligand_colors = grid_col, order=order, gaps = gaps))

}

#' @title Draw a circos plot
#' @usage make_circos_plot(vis_circos_obj, transparency = FALSE, args.circos.text = list(), ...)
#' @description Draw a circos plot
#' @param vis_circos_obj Object returned by \code{\link{prepare_circos_visualization}}
#' @param transparency Logical indicating whether the transparency of the links will correspond to the ligand-target potential score (default: FALSE)
#' @param args.circos.text List of arguments to pass to \code{\link{circos.text}} (by default, the text size is set to 1)
#' @param ... Additional arguments to pass to \code{\link{chordDiagram}}
#' @return A circos plot
#' @export
#' @examples
#' \dontrun{
#' # Default
#' make_circos_plot(vis_circos_obj, transparency = FALSE)
#'
#' # Transparency
#' make_circos_plot(vis_circos_obj, transparency = TRUE)
#'
#' # Make text smaller
#' make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5))
#'
#' # Don't sort links of each ligand based on widths (not recommended)
#' make_circos_plot(vis_circos_obj, transparency = TRUE, args.circos.text = list(cex = 0.5), link.sort = FALSE)
#' }
make_circos_plot <- function(vis_circos_obj, transparency = FALSE, args.circos.text = list(), ...){
  # Check that transparency is a logical
  if (!is.logical(transparency)) stop("transparency should be a logical")

  # Check that vis_circos_obj contains the required elements
  if (!all(c("links_circle", "ligand_colors", "order", "gaps") %in% names(vis_circos_obj))) stop("vis_circos_obj should contain the elements 'links_circle', 'ligand_colors', 'order' and 'gaps'")

  # Check that all elements of args.circos.text is part of the arguments of circos.text
  if (!all(names(args.circos.text) %in% names(formals(circos.text)))) {
    warning("args.circos.text contain element(s) that are not part of the arguments of circos.text")
  }

  # Check that all elements of ... is part of the arguments of chordDiagram
  if (!all(names(list(...)) %in% names(formals(chordDiagram)))) {
    warning("extra arguments contain element(s) that are not part of the arguments of chordDiagram")
  }

  # give the option that links in the circos plot will be transparant ~ ligand-target potential score
  if (!transparency){
    transparency_val <- 0
  } else if (transparency) {
    transparency_val <- vis_circos_obj$links_circle %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency
  }

  default_params <- list(x = vis_circos_obj$links_circle,
                         order=vis_circos_obj$order,
                         grid.col = vis_circos_obj$ligand_colors,
                         transparency = transparency_val,
                         directional = 1, link.sort = TRUE,
                         link.decreasing = FALSE,
                         diffHeight = 0.005,
                         direction.type = c("diffHeight", "arrows"),
                         link.arr.type = "big.arrow",
                         link.visible = vis_circos_obj$links_circle$weight >= attr(vis_circos_obj$links_circle, "cutoff_include_all_ligands"),
                         annotationTrack = "grid",
                         preAllocateTracks = list(track.height = 0.075)
                         )

  # Replace this with user arguments
  default_params[names(list(...))] =  list(...)

  circos_text_default_params <- list(
    facing = "clockwise",
    niceFacing = TRUE,
    adj = c(0, 0.55),
    cex = 1
  )
  circos_text_default_params[names(args.circos.text)] <-  args.circos.text

  # Only the widths of the blocks that indicate each target gene is proportional the ligand-target regulatory potential (~prior knowledge supporting the regulatory interaction).
  circos.par(gap.degree = vis_circos_obj$gaps)
  do.call(chordDiagram, default_params)

  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    do.call(circos.text, c(list(x=CELL_META$xcenter, y=CELL_META$ylim[1], label=CELL_META$sector.index), circos_text_default_params)
           )
  }, bg.border = NA)
  circos.clear()
}

#' @title Make a line plot
#' @usage make_line_plot(ligand_activities, potential_ligands, ranking_range = c(1, 20), agnostic_color = "tomato", focused_color = "black", tied_color = "gray75", inset_scale = 1)
#' @description Make a line plot comparing the ranking of ligands based on their activities in the sender-agnostic and sender-focused approaches
#' @param ligand_activities Dataframe containing the ligand activities from the sender-agnostic approach
#' @param potential_ligands Character vector containing the ligands that are expressed in the sender cell type (i.e. the ligands that are used in the sender-focused approach)
#' @param ranking_range Numeric vector of length 2 indicating the range of the rankings to be displayed (default: c(1, 20))
#' @param agnostic_color Color representing ligands only inthe sender-agnostic approach (default: "tomato")
#' @param focused_color Color representing expressed ligands from the sender-focused approach (default: "black")
#' @param tied_color Color to shade ligands that are tied in the same rank (default: "gray75")
#' @param inset_scale Numeric value indicating the size of the points and text in the inset (default: 1)
#' @return A ggplot object showing the distribution of sender-focused ligands, as well as a line plot inset comparing the rankings between the two approaches
#' @examples \dontrun{
#' # Default
#' make_line_plot(ligand_activities, potential_ligands)
#' }
#' @export

make_line_plot <- function(ligand_activities, potential_ligands, ranking_range = c(1, 20),
                           agnostic_color = "tomato", focused_color = "black", tied_color = "gray75",
                           inset_scale = 1) {

  inset_text_size <- ggplot2::GeomLabel$default_aes$size*inset_scale
  axis_text_size <- ggplot2::GeomLabel$default_aes$size*0.75*inset_scale
  axis_title_size <- ggplot2::GeomLabel$default_aes$size*inset_scale
  point_size <- ggplot2::GeomPoint$default_aes$size*inset_scale
  segment_linewidth <- ggplot2::GeomSegment$default_aes$linewidth*inset_scale
  nudge_x <- 0.05/inset_scale

  # Check if all potential ligands are in ligand_activities
  if (!all(potential_ligands %in% ligand_activities$test_ligand)) stop("Not all potential ligands are in ligand_activities")

  # Check if ranking_range is small -> large
  if (ranking_range[1] >= ranking_range[2]) stop("Starting range should be smaller than ending range")

  # Add rank to ligand_activities if it doesn't exist
  if (!"rank" %in% names(ligand_activities)) ligand_activities <- ligand_activities %>% mutate(rank = rank(desc(aupr_corrected)))

  # x position of "sender-agnostic" ligands
  agnostic_x <- 3.25
  focused_x <- agnostic_x+(1.5*inset_scale)

  # Create dataframe of the two approaches
  rankings_df <- bind_rows(ligand_activities %>% select(test_ligand, rank) %>% mutate(type = "agnostic"),
                          ligand_activities %>% filter(test_ligand %in% potential_ligands) %>%
                            select(test_ligand, rank) %>% mutate(type = "focused")) %>%
    group_by(type) %>% mutate(new_rank = 1:n(),
                              x = case_when(type == "agnostic" ~ agnostic_x,
                                            type == "focused" ~ focused_x)) %>%
    # Ligands that are expressed
    group_by(test_ligand) %>% mutate(expressed = (n() > 1)) %>% ungroup()

  # Define some variables
  start_n <- ranking_range[1]
  end_n <- ranking_range[2]
  n_ligands <- max(rankings_df$new_rank)
  margin <- 1/10                                                    # Leave 1/10 of plot empty at top and bottom
  by_n <- ((n_ligands*margin*9)-(n_ligands*margin))/(end_n-start_n) # Space between each rank
  #cutoff <- by_n+((end_n-start_n+0.5)*by_n) # Doesn't work for all cases

  # Set index to 1 for the start_n ligand
  rankings_df <- rankings_df %>%
    group_by(type) %>%
          mutate(index = (-start_n+2):(n()-start_n+1),
          y = (n_ligands*margin)+(by_n*(index-1)))

  cutoff <- (rankings_df %>% filter(new_rank == end_n, type == "agnostic") %>% pull(y)) + (by_n*0.25)

  # Line segments that go beyond the inset
  line_df <- rankings_df %>% filter(expressed) %>% select(-c(rank, new_rank, index, expressed)) %>%
    pivot_wider(names_from = type, values_from = c(x, y)) %>%
    rename(x1 = x_agnostic, y1 = y_agnostic, x2 = x_focused, y2 = y_focused) %>%
    # Use equation of a line to find the x value at the cutoff
    # Different lines for the top and bottom cutoff
    mutate(m = (y2-y1)/(x2-x1), x0 = case_when(y2 > by_n ~ ((cutoff-y1)/m)+x1,
                                              y2 <= by_n ~ (((n_ligands*margin)-y1)/m)+x1)) %>%
    filter(m != 0)

  # Highlight ties
  ties_df <- rankings_df %>% group_by(type, rank) %>%
    mutate(ties = n() > 1) %>% filter(ties == TRUE) %>%
    ungroup() %>% group_split(type) %>% lapply(., function(group) {
      group %>% split(f = .$rank) %>%
        sapply(., function (k) data.frame(range(k$y))) %>% bind_cols %>% t() %>% data.frame() %>%
        `colnames<-`(c("ystart", "yend")) %>% `rownames<-`(NULL) %>%
        mutate(xstart = case_when(unique(group$type) == "agnostic" ~ agnostic_x,
                                  unique(group$type) == "focused" ~ focused_x),
              xend = xstart)
    }) %>% bind_rows()

  if (nrow(ties_df) > 0){
    # Clip the lines to the cutoff
    ties_df <- ties_df %>% filter(ystart < cutoff, yend > by_n) %>%
      mutate(yend = case_when(yend > cutoff ~ cutoff,
                              TRUE ~ yend),
             ystart = case_when(ystart <= (n_ligands*margin) ~ (n_ligands*margin)-(by_n*0.5),
                                TRUE ~ ystart))
  } else {
    ties_df <- data.frame(ystart = numeric(), yend = numeric(), xstart = numeric(), xend = numeric())
  }

  # Subset the dataframe to the range of interest
  rankings_df_subset <- rankings_df %>% filter(new_rank <= end_n, new_rank >= start_n)

  ggplot() +
    # BAR PLOT
    # Base + expressed ligands drawn as line segments
    geom_rect(aes(ymin=1, ymax=n_ligands, xmin=0.5, xmax=1.5), fill = agnostic_color) +
    geom_segment(data = rankings_df %>% filter(type == "agnostic", expressed),
                aes(y=new_rank, yend=new_rank, x=0.5, xend=1.5), color = focused_color, linewidth = segment_linewidth) +
    # y-axis + ticks
    geom_segment(aes(y=0, yend=max(labeling::extended(0, n_ligands, 5)), x=0.3, xend=0.3)) +
    geom_segment(data = (axis_df <- data.frame(x=0.3, xend=0.25, y=labeling::extended(0, n_ligands, 5),
                                  yend=labeling::extended(0, n_ligands, 5))),
                aes(y=y, yend=yend, x=x, xend=xend)) +
    # y-axis ticklabels + title
    geom_text(data=axis_df, aes(x=xend-0.05, y=y, label=y, hjust=1), size=axis_text_size) +
    geom_text(aes(x=0, y=n_ligands/2, label="Ligand rankings", vjust=-1/inset_scale), angle=90, size=axis_title_size) +
    # Title
    #geom_text(aes(x=1, y=0, label = "Distribution of expressed ligands\nacross all sender-agnostic ligands"), nudge_y=by_n) +
    ggtitle("Distribution of expressed ligands\nacross all sender-agnostic ligands") +
    # ELBOW CONNECTORS
    # Top, vertical line, bottom, horizontal line connecting to inset
    geom_segment(data = data.frame(x = c(1.6, 1.65, 1.65, 1.65),
                                  xend = c(1.65, 1.65, 1.60, agnostic_x-1.25),
                                  y=c(min(rankings_df_subset$new_rank), min(rankings_df_subset$new_rank), max(rankings_df_subset$new_rank), ((end_n-start_n)/2)+start_n-1),
                                  yend=c(min(rankings_df_subset$new_rank), max(rankings_df_subset$new_rank), max(rankings_df_subset$new_rank), ((end_n-start_n)/2)+start_n-1)),
                aes(x=x, xend=xend, y=y, yend=yend)) +
    # LINE PLOT
    # Ties
    geom_segment(data = ties_df, aes(x=xstart, y=ystart, xend=xend, yend=yend),
                color = tied_color, linewidth=3, lineend="round") +
    # Points
    geom_point(data = rankings_df_subset, aes(x=x, y=y, color = expressed), size = point_size) +
    # Line segment from focused -> agnostic
    geom_segment(data = line_df %>% filter(x0 > agnostic_x, y2 < cutoff, y2 >= (n_ligands*margin)), aes(x=x2, y=y2, xend=x0, yend=cutoff), linewidth = segment_linewidth) +
    # Line segment from agnostic -> focused
    geom_segment(data = line_df %>% filter(x0 > agnostic_x, y1 < cutoff, y1 >= (n_ligands*margin)), aes(x=x1, y=y1, xend=x0, yend=(n_ligands*margin)), linewidth = segment_linewidth) +
    # Line segment within range
    geom_line(data = rankings_df_subset, aes(x=x, y=y, group = test_ligand), linewidth = segment_linewidth) +
    # Ligand names
    geom_text(data = rankings_df_subset %>% filter(type == "agnostic"), aes(x=x, y=y, label = test_ligand, hjust = "right", color = expressed), nudge_x = -nudge_x, size=inset_text_size) +
    geom_text(data = rankings_df_subset %>% filter(type == "focused"), aes(x=x, y=y, label = test_ligand, hjust = "left",  color = expressed), nudge_x = nudge_x, size=inset_text_size) +
    # Ranking labels
    geom_text(data = rankings_df_subset, aes(x=agnostic_x-0.75, y=y, label=new_rank), size=inset_text_size) +
    # Outer rectangle
    geom_rect(aes(ymin=0, ymax=n_ligands*(margin*19/2), xmin=agnostic_x-1.25, xmax = agnostic_x+2.25), fill=NA, color="black") +
     # Heading
    geom_text(data = data.frame(x = c(agnostic_x-0.75, agnostic_x, focused_x),
                                y = n_ligands*(margin/2), label = c("Rank", "Agnostic", "Focused")),
                                aes(x=x, y=y, label=label), size=inset_text_size) +
    # PLOT SETTINGS
    scale_color_manual(values = c("TRUE" = focused_color, "FALSE" = agnostic_color), breaks=c(TRUE, FALSE), labels = c("Only agnostic", "Tied")) +
    scale_y_reverse() +
    xlim(0, agnostic_x+2.5) +
    labs(y = "Ligand rankings") +
    guides(color = guide_legend(override.aes = list(shape = c(19, 15), size = c(2, 4), color = c(agnostic_color, tied_color)))) +
    theme_classic() +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.line = element_blank(), axis.ticks = element_blank(),
          legend.title = element_blank(), legend.direction = "horizontal",
          legend.position = c(0.5, 0.05),
          legend.background = element_blank(),
          legend.text = element_text(size = 12),
          plot.title = element_text(hjust=0.11, margin=margin(0, 0, -20, 0)))

}

