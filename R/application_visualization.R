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
get_ligand_signaling_path = function(ligand_tf_matrix, ligands_all, targets_all, top_n_regulators = 4, weighted_networks, ligands_position = "cols"){

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

  matrix_df_vis = matrix %>% data.frame() %>% rownames_to_column("y") %>% as_tibble() %>% gather(x,"score", -y) %>% mutate(y = factor(y, levels = rownames(matrix), ordered = TRUE), x = factor(x, levels = colnames(matrix), ordered = TRUE))

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

  matrix_df_vis = matrix %>% data.frame() %>% rownames_to_column("y") %>% as_tibble() %>% gather(x,"score", -y) %>% mutate(y = factor(y, levels = rownames(matrix), ordered = TRUE), x = factor(x, levels = colnames(matrix), ordered = TRUE))

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
#' @usage
#' make_mushroom_plot(prioritization_table, top_n = 30, show_ranking = FALSE, show_all_datapoints = FALSE, true_color_range = FALSE, size = "scaled_avg_exprs", color = "scaled_lfc",
#'  ligand_fill_colors = c("#DEEBF7", "#08306B"), receptor_fill_colors = c("#FEE0D2", "#A50F15"),
#'  unranked_ligand_fill_colors = c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)), unranked_receptor_fill_colors = c( alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))
#'
#' @param prioritization_table A prioritization table as generated by \code{\link{generate_prioritization_tables}}
#' @param top_n An integer indicating how many ligand-receptor pairs to show
#' @param show_ranking A logical indicating whether to show the ranking of the ligand-receptor pairs (default: FALSE)
#' @param show_all_datapoints A logical indicating whether to show all ligand-receptor pairs (default: FALSE, if true they will be grayed out)
#' @param true_color_range A logical indicating whether to use the true color range for the ligand-receptor pairs (default: FALSE; range 0-1 is used)
#' @param size A string indicating which column to use for the size of the semicircles (default: "scaled_avg_exprs"; use column name without "_ligand" or "_receptor" suffix)
#' @param color A string indicating which column to use for the color of the semicircles (default: "scaled_lfc"; use column name without "_ligand" or "_receptor" suffix)
#' @param ligand_fill_colors A vector of the low and high colors to use for the ligand semicircle fill gradient (default: c("#DEEBF7", "#08306B"))
#' @param receptor_fill_colors A vector of the low and high colors to use for the receptor semicircle fill gradient (default: c("#FEE0D2", "#A50F15"))
#' @param unranked_ligand_fill_colors A vector of the low and high colors to use for the unranked ligands when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))
#' @param unranked_receptor_fill_colors A vector of the low and high colors to use for the unkraed receptors when show_all_datapoints is TRUE (default: c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)))
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
#' # Show only top 20, and write rankings on the plot
#' make_mushroom_plot(prior_table, top_n = 20, show_ranking = TRUE)
#' 
#' # Show all datapoints, and use true color range
#' make_mushroom_plot(prior_table, show_all_datapoints = TRUE, true_color_range = TRUE)
#' 
#' # Change the size and color columns
#' make_mushroom_plot(prior_table, size = "pct_expressed", color = "scaled_avg_exprs")
#' }
#' @export
#'
make_mushroom_plot = function(prioritization_table, top_n = 30, show_rankings = FALSE,
                              show_all_datapoints = FALSE, true_color_range = FALSE,
                              size = "scaled_avg_exprs", color = "scaled_lfc",
                              ligand_fill_colors = c("#DEEBF7", "#08306B"),
                              receptor_fill_colors = c("#FEE0D2", "#A50F15"),
                              unranked_ligand_fill_colors = c(alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2)),
                              unranked_receptor_fill_colors = c( alpha("#FFFFFF", alpha=0.2), alpha("#252525", alpha=0.2))){
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

  # Filter to top_n, create a new column of ligand-receptor interactions
  filtered_table <- prioritization_table %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>%
    dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = " - "))
  order_interactions <- unique(filtered_table %>% filter(prioritization_rank <= top_n) %>% pull(lr_interaction))
  filtered_table <- filtered_table %>% filter(lr_interaction %in% order_interactions) %>%
    mutate(lr_interaction = factor(lr_interaction, levels = rev(order_interactions)))

  celltypes_vec <- 1:length(unique(filtered_table$sender)) %>% setNames(sort(unique(filtered_table$sender)))
  lr_interaction_vec <- 1:length(order_interactions) %>% setNames(order_interactions)

  # Make each ligand and receptor into separate rows (to draw 1 semicircle per row)
  filtered_table <- filtered_table %>% select(c("lr_interaction", all_of(cols_to_use), "prioritization_rank")) %>%
    pivot_longer(c(ligand, receptor), names_to = "type", values_to = "protein") %>%
    mutate(size = ifelse(type == "ligand", get(paste0(size, "_", size_ext[1])), get(paste0(size, "_", size_ext[2]))),
           color = ifelse(type == "ligand", get(paste0(color, "_", color_ext[1])), get(paste0(color, "_",  color_ext[2])))) %>%
    select(-contains(c("_ligand", "_receptor", "_sender", "_receiver"))) %>%
    mutate(start = rep(c(-pi, 0), nrow(filtered_table))) %>%
    mutate(x = celltypes_vec[sender], y = lr_interaction_vec[lr_interaction])

  # Rename size and color columns to be more human-readable
  keywords_adj <- c("LFC", "p-val", "product", "mean", "adjusted", "expression") %>% setNames(c("lfc", "pval", "prod", "avg", "adj", "exprs"))
  size_title <- sapply(stringr::str_split(size, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>%  R.utils::capitalize()
  color_title <- sapply(stringr::str_split(color, "_")[[1]], function(k) ifelse(is.na(keywords_adj[k]), k, keywords_adj[k])) %>%
    paste0(., collapse = " ") %>% R.utils::capitalize()

  color_lims <- c(0,1)
  if (true_color_range) color_lims <- NULL

  scale <- 0.5
  p1 <- ggplot() +
    # Draw ligand semicircle
    geom_arc_bar(data = filtered_table %>% filter(type=="ligand",  prioritization_rank <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    scale_fill_gradient(low = ligand_fill_colors[1] , high=ligand_fill_colors[2] ,
                        limits=color_lims, oob=scales::squish,
                        name=paste0(color_title, " (", color_ext[1], ")") %>% str_wrap(width=15)) +
    # Create new fill scale for receptor semicircles
    new_scale_fill() +
    geom_arc_bar(data = filtered_table %>% filter(type=="receptor",  prioritization_rank <= top_n),
                 aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                     start = start, end = start + pi, fill=color),
                 color = "white") +
    scale_fill_gradient(low = receptor_fill_colors[1], high=receptor_fill_colors[2] , limits=color_lims, oob=scales::squish,
                        name=paste0(color_title,  " (", color_ext[2], ")") %>% str_wrap(width=15)) +
    # Other plot information
    scale_y_continuous(breaks=length(lr_interaction_vec):1, labels=names(lr_interaction_vec)) +
    scale_x_continuous(breaks=1:length(celltypes_vec), labels=names(celltypes_vec), position="top") +
    xlab("Sender cell types") + ylab("Ligand-receptor interaction") +
    coord_fixed() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          legend.box = "horizontal")


  # Add unranked ligand and receptor semicircles if requested
  if (show_all_datapoints){
    p1 <- p1 + new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="ligand", prioritization_rank > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low = unranked_ligand_fill_colors[1], high=unranked_ligand_fill_colors[2],
                          limits=color_lims, oob=scales::squish, guide = "none") +
      new_scale_fill() +
      geom_arc_bar(data = filtered_table %>% filter(type=="receptor", prioritization_rank > top_n),
                   aes(x0 = x, y0 = y, r0 = 0, r = sqrt(size)*scale,
                       start = start, end = start + pi, fill=color),
                   color = "white") +
      scale_fill_gradient(low=unranked_receptor_fill_colors[1], high=unranked_receptor_fill_colors[2],
                          limits=color_lims, oob=scales::squish, guide = "none")
  }

  # Add ranking numbers if requested
  if (show_rankings){
    p1 <- p1 + geom_shadowtext(data = filtered_table %>% filter(prioritization_rank <= top_n),
                               aes(x=x, y=y, label=prioritization_rank))
  }

  legend1 <- ggpubr::as_ggplot(ggpubr::get_legend(p1))

  # For the size legend, create a new plot
  legend2 <- ggplot(data.frame(values = c(0.25, 0.5, 0.75, 1), x=1:4, y=1, start=-pi)) +
    geom_rect(aes(xmin=x-0.5, xmax=x+0.5, ymin=y-0.5, ymax=y+0.5), color="gray80", fill=NA) +
    geom_arc_bar(aes(x0=x, y0=y, r0=0, r=sqrt(values)*scale, start=start, end=start+pi), fill="black") +
    geom_text(aes(label=values, x=x, y=y-0.6), vjust=1) +
    labs(tag = size_title) +
    scale_x_continuous(breaks = 1:4, labels=c(0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(limits=c(-0.5, 1.5)) +
    labs(x="Percent expressed") +
    coord_fixed() +  theme_classic() +
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          plot.margin = margin(0, 0, 10, 0),
          plot.tag.position = "top",
          plot.tag = element_text(margin=margin(0, 0, 5,0), size=10),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank())

  # Combine the two legends
  legends <- cowplot::plot_grid(NULL, legend1, legend2, NULL, nrow=4, scale=c(1,1,0.5,1),
                                rel_heights = c(2, 1, 2, 2), align = "v", axis="tb")
  cowplot::plot_grid(p1 + theme(legend.position="none"), legends)
}



