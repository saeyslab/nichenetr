#' @title Convert cluster assignment to settings format suitable for target gene prediction.
#'
#' @description \code{convert_cluster_to_settings} Convert cluster assignment to settings format suitable for target gene prediction.
#'
#' @usage
#' convert_cluster_to_settings(i, cluster_vector, setting_name, setting_from, background = NULL)
#'
#' @param i The cluster number of the cluster of interest to which genes should belong
#' @param cluster_vector Named vector containing the cluster number to which every gene belongs
#' @param setting_name Base name of the setting
#' @param setting_from Active ligands for the specific setting
#' @param background NULL or a character vector of genes belonging to the background. When NULL: the background will be formed by genes belonging to other clusters that the cluster of interest. Default NULL. If not NULL and genes present in the cluster of interest are in this vector of background gene names, these genes will be removed from the background.
#'
#' @return A list with following elements: $name (indicating the cluster id), $from, $response. $response is a gene-named logical vector indicating whether the gene is part of the respective cluster.
#'
#' @examples
#' \dontrun{
#' genes_clusters = c("TGFB1" = 1,"TGFB2" = 1,"TGFB3" = 2)
#' cluster_settings = lapply(seq(length(unique(genes_clusters))), convert_cluster_to_settings, cluster_vector = genes_clusters, setting_name = "example", setting_from = "BMP2")
#' }
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
    background = background[(background %in% genes_cluster_oi) == FALSE]
    background_logical = rep(FALSE,times = length(background))
    names(background_logical) = background
    cluster_logical = rep(TRUE,times = length(genes_cluster_oi))
    names(cluster_logical) = genes_cluster_oi
    response = c(background_logical,cluster_logical)
  }
  return(list(name = paste0(setting_name,"_cluster_",i), from = setting_from, response = response))
}
#' @title Predict activities of ligands in regulating expression of a gene set of interest
#'
#' @description \code{predict_ligand_activities} Predict activities of ligands in regulating expression of a gene set of interest. Ligand activities are defined as how well they predict the observed transcriptional response (i.e. gene set) according to the NicheNet model.
#'
#' @usage
#' predict_ligand_activities(geneset, background_expressed_genes,ligand_target_matrix, potential_ligands, single = TRUE,...)
#'
#' @param geneset Character vector of the gene symbols of genes of which the expression is potentially affected by ligands from the interacting cell.
#' @param background_expressed_genes Character vector of gene symbols of the background, non-affected, genes (can contain the symbols of the affected genes as well).
#' @param ligand_target_matrix The NicheNet ligand-target matrix denoting regulatory potential scores between ligands and targets (ligands in columns).
#' @param potential_ligands Character vector giving the gene symbols of the potentially active ligands you want to define ligand activities for.
#' @param single TRUE if you want to calculate ligand activity scores by considering every ligand individually (recommended). FALSE if you want to calculate ligand activity scores as variable importances of a multi-ligand classification model.
#' @param ... Additional parameters for get_multi_ligand_importances if single = FALSE.
#'
#' @return A tibble giving several ligand activity scores. Following columns in the tibble: $test_ligand, $auroc, $aupr and $pearson.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' ligand_activities = predict_ligand_activities(geneset = geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#' }
#'
#' @export
#'
predict_ligand_activities = function(geneset,background_expressed_genes,ligand_target_matrix, potential_ligands, single = TRUE,...){
  setting = list(geneset) %>%
    lapply(convert_gene_list_settings_evaluation, name = "gene set", ligands_oi = potential_ligands, background = background_expressed_genes)
  if (single == TRUE){
    settings_ligand_prediction = setting %>%
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE) %>% bind_rows()

  } else {
    settings_ligand_prediction = setting %>%
      convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = FALSE)
    ligand_importances = settings_ligand_prediction %>% lapply(get_multi_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE, ...) %>% bind_rows()

  }
  return(ligand_importances %>% select(test_ligand,auroc,aupr,aupr_corrected, pearson))
}
#' @title Infer weighted active ligand-target links between a possible ligand and target genes of interest
#'
#' @description \code{get_weighted_ligand_target_links} Infer active ligand target links between possible ligands and genes belonging to a gene set of interest: consider the intersect between the top n targets of a ligand and the gene set.
#'
#' @usage
#' get_weighted_ligand_target_links(ligand, geneset,ligand_target_matrix,n = 250)
#'
#' @param geneset Character vector of the gene symbols of genes of which the expression is potentially affected by ligands from the interacting cell.
#' @param ligand Character vector giving the gene symbols of the potentially active ligand for which you want to find target genes.
#' @param n The top n of targets per ligand that will be considered. Default: 250.
#' @inheritParams predict_ligand_activities
#'
#' @return A tibble with columns ligand, target and weight (i.e. regulatory potential score).
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligand = "TNF"
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' active_ligand_target_links_df = get_weighted_ligand_target_links(ligand = potential_ligand, geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 250)
#' }
#'
#' @export
#'
get_weighted_ligand_target_links = function(ligand, geneset,ligand_target_matrix,n = 250){
  top_n_score = ligand_target_matrix[,ligand] %>% sort(decreasing = T) %>% head(n) %>% min()
  targets = intersect(ligand_target_matrix[,ligand] %>% .[. >= top_n_score ] %>% names(),geneset)
  if (length(targets) == 0){
    ligand_target_weighted_df = tibble(ligand = ligand, target = NA, weight = NA)
  } else {
    ligand_target_weighted_df = tibble(ligand = ligand, target = targets, weight = ligand_target_matrix[targets,ligand])
  }
  return(ligand_target_weighted_df)
}
#' @title Prepare heatmap visualization of the ligand-target links starting from a ligand-target tibble.
#'
#' @description \code{prepare_ligand_target_visualization} Prepare heatmap visualization of the ligand-target links starting from a ligand-target tibble. Get regulatory potential scores between all pairs of ligands and targets documented in this tibble. For better visualization, we propose to define a quantile cutoff on the ligand-target scores.
#'
#' @usage
#' prepare_ligand_target_visualization(ligand_target_df, ligand_target_matrix, cutoff = 0.25)
#'
#' @param cutoff Quantile cutoff on the ligand-target scores of the input weighted ligand-target network. Scores under this cutoff will be set to 0.
#' @param ligand_target_df Tibble with columns 'ligand', 'target' and 'weight' to indicate ligand-target regulatory potential scores of interest.
#' @inheritParams predict_ligand_activities
#'
#' @return A matrix giving the ligand-target regulatory potential scores between ligands of interest and their targets genes part of the gene set of interest.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' active_ligand_target_links_df = potential_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
#' active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
#' }
#'
#' @export
#'
prepare_ligand_target_visualization = function(ligand_target_df, ligand_target_matrix, cutoff = 0.25){

  # define a cutoff on the ligand-target links
  cutoff_include_all_ligands = ligand_target_df$weight %>% quantile(cutoff)

  # give a score of 0 to ligand-target links not higher than the defined cutoff
  ligand_target_matrix_oi = ligand_target_matrix
  ligand_target_matrix_oi[ligand_target_matrix_oi < cutoff_include_all_ligands] = 0

  # consider only targets belonging to the top250 targets of individual ligands and with at least one ligand-link with score higher than the defined cutoff
  ligand_target_vis = ligand_target_matrix_oi[ligand_target_df$target %>% unique(),ligand_target_df$ligand %>% unique()]
  dim(ligand_target_vis) = c(length(ligand_target_df$target %>% unique()), length(ligand_target_df$ligand %>% unique()))
  all_targets = ligand_target_df$target %>% unique()
  all_ligands = ligand_target_df$ligand %>% unique()
  rownames(ligand_target_vis) = all_targets
  colnames(ligand_target_vis) = all_ligands

  keep_targets = all_targets[ligand_target_vis %>% apply(1,sum) > 0]
  keep_ligands = all_ligands[ligand_target_vis %>% apply(2,sum) > 0]


  ligand_target_vis_filtered = ligand_target_vis[keep_targets,keep_ligands]


  if(is.matrix(ligand_target_vis_filtered)){
    rownames(ligand_target_vis_filtered) = keep_targets
    colnames(ligand_target_vis_filtered) = keep_ligands

  } else {
    dim(ligand_target_vis_filtered) = c(length(keep_targets), length(keep_ligands))
    rownames(ligand_target_vis_filtered) = keep_targets
    colnames(ligand_target_vis_filtered) = keep_ligands
  }

  if(nrow(ligand_target_vis_filtered) > 1 & ncol(ligand_target_vis_filtered) > 1){
    distoi = dist(1-cor(t(ligand_target_vis_filtered)))
    hclust_obj = hclust(distoi, method = "ward.D2")
    order_targets = hclust_obj$labels[hclust_obj$order]

    distoi_targets = dist(1-cor(ligand_target_vis_filtered))
    hclust_obj = hclust(distoi_targets, method = "ward.D2")
    order_ligands = hclust_obj$labels[hclust_obj$order]

  } else {
    order_targets = rownames(ligand_target_vis_filtered)
    order_ligands = colnames(ligand_target_vis_filtered)
  }

  vis_ligand_target_network = ligand_target_vis_filtered[order_targets,order_ligands]
  dim(vis_ligand_target_network) = c(length(order_targets), length(order_ligands))
  rownames(vis_ligand_target_network) = order_targets
  colnames(vis_ligand_target_network) = order_ligands
  return(vis_ligand_target_network)

}

#' @title Get the weighted ligand-receptor links between a possible ligand and its receptors
#' @description \code{get_weighted_ligand_receptor_links} Get the weighted ligand-receptor links between a possible ligand and its receptors
#'
#' @param best_upstream_ligands Character vector containing ligands of interest.
#' @param expressed_receptors Character vector of receptors expressed in the cell type of interest.
#' @param lr_network A data frame with two columns, \code{from} and \code{to}, containing the ligand-receptor interactions.
#' @param weighted_networks_lr_sig A data frame with three columns, \code{from}, \code{to} and \code{weight}, containing the ligand-receptor interactions and their weights.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   ligand_receptor_links_df <- get_weighted_ligand_receptor_links(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks$lr_sig)
#' }
#'
get_weighted_ligand_receptor_links = function(best_upstream_ligands, expressed_receptors, lr_network, weighted_networks_lr_sig) {

  lr_network <- lr_network %>% distinct(from, to)
  weighted_networks_lr <- inner_join(weighted_networks_lr_sig, lr_network, by = c("from","to"))

  lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_long <- weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  return(lr_network_top_df_long)

}

#' @title Prepare ligand-receptor visualization
#' @description \code{prepare_ligand_receptor_visualization} Prepare a matrix of ligand-receptor interactions for visualization.
#'
#' @param lr_network_top_df_long A data frame with three columns, \code{from}, \code{to} and \code{weight}, containing the ligand-receptor interactions and their weights.
#' @param best_upstream_ligands Character vector of ligands of interest. This will only be used if \code{order_hclust = "receptors"} or \code{order_hclust = "none"}.
#' @param order_hclust "both", "ligands", "receptors", or "none". If "both", the ligands and receptors are ordered by hierarchical clustering. If "ligands" or "receptors" only the ligands or receptors are ordered hierarchically. If "none", no hierarchical clustering is performed, and the ligands are ordered based on \code{best_upstream_ligands}, and the receptors are ordered alphabetically.
#'
#' @return A matrix of ligand-receptor interactions for visualization.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ligand_receptor_network = prepare_ligand_receptor_visualization(best_upstream_ligands = best_upstream_ligands, expressed_receptors = expressed_receptors, lr_network = lr_network, weighted_networks_lr_sig = weighted_networks_lr_sig, order_hclust = TRUE)
#' }
#'
prepare_ligand_receptor_visualization = function(lr_network_top_df_long, best_upstream_ligands, order_hclust = "both") {

  lr_network_top_df <- lr_network_top_df_long %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  # Check if order_hclust is valid
  if (!(order_hclust %in% c("both", "ligands", "receptors", "none"))) {
    stop("order_hclust must be one of 'both', 'ligands', 'receptors', or 'none'")
  }

  if (order_hclust == "both" | order_hclust == "receptors") {
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
  }

  if (order_hclust == "both" | order_hclust == "ligands") {
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  }

  if (order_hclust == "none" | order_hclust == "receptors") {
    order_ligands_receptor = rev(best_upstream_ligands)
  }

  if (order_hclust == "none" | order_hclust == "ligands") {
    order_receptors = rownames(lr_network_top_matrix)
  }

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) <- order_receptors
  colnames(vis_ligand_receptor_network) <- order_ligands_receptor

  return(vis_ligand_receptor_network)

}

#' @title Assess probability that a target gene belongs to the geneset based on a multi-ligand random forest model
#'
#' @description \code{assess_rf_class_probabilities} Assess probability that a target gene belongs to the geneset based on a multi-ligand random forest model (with cross-validation). Target genes and background genes will be split in different groups in a stratified way.
#'
#' @usage
#' assess_rf_class_probabilities(round,folds,geneset,background_expressed_genes,ligands_oi,ligand_target_matrix)
#'
#' @param ligands_oi Character vector giving the gene symbols of the ligands you want to build the multi-ligand with.
#' @param round Integer describing which fold of the cross-validation scheme it is.
#' @param folds Integer describing how many folds should be used.
#' @inheritParams predict_ligand_activities
#'
#' @return A tibble with columns: $gene, $response, $prediction. Response indicates whether the gene belongs to the geneset of interest, prediction gives the probability this gene belongs to the geneset according to the random forest model.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' fold1_rf_prob = assess_rf_class_probabilities(round = 1,folds = 2,geneset = geneset,background_expressed_genes = background_expressed_genes ,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)
#' }
#'
#' @export
#'
assess_rf_class_probabilities = function(round,folds,geneset,background_expressed_genes,ligands_oi, ligand_target_matrix){
  set.seed(round)
  geneset_shuffled = sample(geneset, size = length(geneset))
  geneset_grouped = split(geneset_shuffled,1:folds)

  strict_background_expressed_genes = background_expressed_genes[!background_expressed_genes %in% geneset]
  set.seed(round)
  strict_background_expressed_genes_shuffled = sample(strict_background_expressed_genes, size = length(strict_background_expressed_genes))
  strict_background_expressed_genes_grouped = split(strict_background_expressed_genes_shuffled,1:folds)

  geneset_predictions_all = seq(length(geneset_grouped)) %>% lapply(rf_target_prediction,geneset_grouped,strict_background_expressed_genes_grouped,ligands_oi,ligand_target_matrix) %>% bind_rows()
  geneset_predictions_all = geneset_predictions_all %>% mutate(response = gsub("\\.","",response) %>% as.logical())
}
#' @title Assess how well classification predictions accord to the expected response
#'
#' @description \code{classification_evaluation_continuous_pred_wrapper} Assess how well classification predictions accord to the expected response.
#'
#' @usage
#' classification_evaluation_continuous_pred_wrapper(response_prediction_tibble)
#'
#' @param response_prediction_tibble Tibble with columns "response" and "prediction" (e.g. output of function `assess_rf_class_probabilities`)
#'
#' @return A tibble showing several classification evaluation metrics.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' fold1_rf_prob = assess_rf_class_probabilities(round = 1,folds = 2,geneset = geneset,background_expressed_genes = background_expressed_genes ,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)
#  classification_evaluation_continuous_pred_wrapper(fold1_rf_prob)
#' }
#'
#' @export
#'
classification_evaluation_continuous_pred_wrapper = function(response_prediction_tibble) {
  prediction_performances = classification_evaluation_continuous_pred(response_prediction_tibble$prediction, response_prediction_tibble$response, iregulon = FALSE)
  return(prediction_performances)
}
#' @title Find which genes were among the top-predicted targets genes in a specific cross-validation round and see whether these genes belong to the gene set of interest as well.
#'
#' @description \code{get_top_predicted_genes} Find which genes were among the top-predicted targets genes in a specific cross-validation round and see whether these genes belong to the gene set of interest as well.
#'
#' @usage
#' get_top_predicted_genes(round,gene_prediction_list, quantile_cutoff = 0.95)
#'
#' @param gene_prediction_list List with per round of cross-validation: a tibble with columns "gene", "prediction" and "response" (e.g. output of function `assess_rf_class_probabilities`)
#' @param round Integer describing which fold of the cross-validation scheme it is.
#' @param quantile_cutoff Quantile of which genes should be considered as top-predicted targets. Default: 0.95, thus considering the top 5 percent predicted genes as predicted targets.
#'
#' @return A tibble indicating for every gene whether it belongs to the geneset and whether it belongs to the top-predicted genes in a specific cross-validation round.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' gene_predictions_list = seq(2) %>% lapply(assess_rf_class_probabilities,2, geneset = geneset,background_expressed_genes = background_expressed_genes,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)
#' seq(length(gene_predictions_list))  %>% lapply(get_top_predicted_genes,gene_predictions_list)
#' }
#'
#' @export
#'
get_top_predicted_genes = function(round,gene_prediction_list, quantile_cutoff = 0.95){
  affected_gene_predictions = gene_prediction_list[[round]]
  predicted_positive = affected_gene_predictions %>%
    arrange(-prediction) %>%
    mutate(predicted_top_target = prediction >= quantile(prediction,quantile_cutoff)) %>%
    filter(predicted_top_target) %>% rename(true_target = response) %>%
    select(gene,true_target,predicted_top_target)
  colnames(predicted_positive) = c("gene","true_target",paste0("predicted_top_target_round",round))
  return(predicted_positive)
}
#' @title Determine the fraction of genes belonging to the geneset or background and to the top-predicted genes.
#'
#' @description \code{calculate_fraction_top_predicted} Defines the fraction of genes belonging to the geneset or background and to the top-predicted genes.
#'
#' @usage
#' calculate_fraction_top_predicted(affected_gene_predictions, quantile_cutoff = 0.95)
#'
#' @param affected_gene_predictions Tibble with columns "gene", "prediction" and "response" (e.g. output of function `assess_rf_class_probabilities`)
#' @param quantile_cutoff Quantile of which genes should be considered as top-predicted targets. Default: 0.95, thus considering the top 5 percent predicted genes as predicted targets.
#'
#' @return A tibble indicating the number of genes belonging to the gene set of interest or background (true_target column), the number and fraction of genes of these gruops that were part of the top predicted targets in a specific cross-validation round.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' gene_predictions_list = seq(2) %>% lapply(assess_rf_class_probabilities,2, geneset = geneset,background_expressed_genes = background_expressed_genes,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)
#' target_prediction_performances_discrete_cv = gene_predictions_list %>% lapply(calculate_fraction_top_predicted) %>% bind_rows() %>% ungroup() %>% mutate(round=rep(1:length(gene_predictions_list), each = 2))

#' }
#'
#' @export
#'
calculate_fraction_top_predicted = function(affected_gene_predictions, quantile_cutoff = 0.95){
  predicted_positive = affected_gene_predictions %>% arrange(-prediction) %>% filter(prediction >= quantile(prediction,quantile_cutoff)) %>% group_by(response) %>% count() %>% rename(positive_prediction = n) %>% rename(true_target = response)
  all = affected_gene_predictions %>% arrange(-prediction) %>% rename(true_target = response) %>% group_by(true_target) %>% count()
  inner_join(all,predicted_positive, by = "true_target") %>% mutate(fraction_positive_predicted = positive_prediction/n)
}
#' @title Perform a Fisher's exact test to determine whether genes belonging to the gene set of interest are more likely to be part of the top-predicted targets.
#'
#' @description \code{calculate_fraction_top_predicted_fisher} Performs a Fisher's exact test to determine whether genes belonging to the gene set of interest are more likely to be part of the top-predicted targets.
#'
#' @usage
#' calculate_fraction_top_predicted_fisher(affected_gene_predictions, quantile_cutoff = 0.95, p_value_output = TRUE)
#'
#' @param p_value_output Should total summary or p-value be returned as output? Default: TRUE.
#' @inheritParams calculate_fraction_top_predicted
#'
#' @return Summary of the Fisher's exact test or just the p-value
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' geneset = c("SOCS2","SOCS3", "IRF1")
#' background_expressed_genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' gene_predictions_list = seq(2) %>% lapply(assess_rf_class_probabilities,2, geneset = geneset,background_expressed_genes = background_expressed_genes,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)
#' target_prediction_performances_fisher_pval = gene_predictions_list %>% lapply(calculate_fraction_top_predicted_fisher) %>% unlist() %>% mean()
#' }
#'
#' @export
#'
calculate_fraction_top_predicted_fisher = function(affected_gene_predictions, quantile_cutoff = 0.95, p_value_output = TRUE){
  predicted_positive = affected_gene_predictions %>% arrange(-prediction) %>% filter(prediction >= quantile(prediction,quantile_cutoff)) %>% group_by(response) %>% count() %>% rename(positive_prediction = n)
  all = affected_gene_predictions %>% arrange(-prediction)  %>% group_by(response) %>% count()
  results_df = left_join(all, predicted_positive, by="response") %>% mutate(positive_prediction = replace_na(positive_prediction, 0))

  tp = results_df %>% filter(response == TRUE) %>% .$positive_prediction
  fp = results_df %>% filter(response == FALSE) %>% .$positive_prediction
  fn = (results_df %>% filter(response == TRUE) %>% .$n) - (results_df %>% filter(response == TRUE) %>% .$positive_prediction)
  tn = (results_df %>% filter(response == FALSE) %>% .$n) - (results_df %>% filter(response == FALSE) %>% .$positive_prediction)

  contingency_table = matrix(c(tp,fp,fn,tn), nrow = 2,dimnames = list(c("geneset", "background"), c("top-predicted", "no-top-predicted")))
  summary = fisher.test(contingency_table, alternative = "greater")
  if(p_value_output == TRUE){
    return(summary$p.value)
  } else {
    return(summary)
  }
}

#' @title Run ligand activity analysis with bootstrap
#' @description \code{bootstrap_ligand_activity_analysis} Randomly sample a gene set from all expressed genes in the receiver cell type, then perform ligand activity analysis on this gene set. This 'null gene set' has equal length to the gene set of interest.
#' @usage bootstrap_ligand_activity_analysis(expressed_genes_receiver, geneset_oi, background_expressed_genes, ligand_target_matrix, potential_ligands, n_iter = 10, n_cores = 1, parallel_func = "mclapply")
#' @param expressed_genes_receiver Genes expressed in the receiver cell type
#' @inheritParams predict_ligand_activities
#' @param n_iter Number of iterations to perform (Default: 10)
#' @param geneset_oi Character vector of the gene symbols of genes of which the expression is potentially affected by ligands from the interacting cell.
#' @param n_cores Number of cores to use for parallelization (Default: 1)
#' @param parallel_func Parallelization function to use from "mclapply", "pbmclapply", or "parLapply" (Default: "mclapply")
#' @return List of n_iter elements, each element containing the output of predict_ligand_activities for a random gene set
#' @examples
#' \dontrun{
#' permutations <- bootstrap_ligand_activity_analysis(expressed_genes_receiver, geneset_oi, background_expressed_genes,
#'                            ligand_target_matrix, potential_ligands, n_iter = 10, n_cores = 1, parallel_func = "mclapply")
#' }
#' @export
bootstrap_ligand_activity_analysis <- function(expressed_genes_receiver, geneset_oi, background_expressed_genes,
                                               ligand_target_matrix, potential_ligands, n_iter = 10,
                                               n_cores=1, parallel_func = "mclapply"){

  # Check if parallel function is valid
  if (!(parallel_func %in% c("mclapply", "pbmclapply", "parLapply"))) {
    stop("parallel_func must be one of 'mclapply', 'pbmclapply', or 'parLapply'")
  }

  # Check if parallel function is mclapply/pbmclapply that they are not on Windows
  if (Sys.info()[['sysname']] == "Windows" && parallel_func %in% c("mclapply", "pbmclapply")){
    if (n_cores > 1){
      warning("Setting 'n_cores' to 1 as Windows OS cannot make use of mclapply and pbmclapply parallelization. Consider using 'parLapply' instead.")
      n_cores = 1
    }
  }

  if (parallel_func %in% c("mclapply", "pbmclapply")){
    parFunc <- ifelse(parallel_func == "mclapply", parallel::mclapply, pbmcapply::pbmclapply)

    perms <- parFunc(1:n_iter, function (i){
      random_geneset <- sample(expressed_genes_receiver, size = length(geneset_oi))
      predict_ligand_activities(geneset = random_geneset, background_expressed_genes = background_expressed_genes,
                                ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    }, mc.cores = n_cores)

  } else if (parallel_func == "parLapply"){
    if(Sys.info()[['sysname']] == "Windows"){
      # STILL HAS TO BE TESTED ON WINDOWS
      clust <- parallel::makeCluster(n_cores)
      parallel::clusterExport(clust, c("expressed_genes_receiver", "geneset_oi","background_expressed_genes","ligand_target_matrix","potential_ligands"), envir = environment())
      parallel::clusterEvalQ(clust, library(nichenetr))
      parallel::clusterEvalQ(clust, library(tidyverse))

    } else {
      clust <- parallel::makeCluster(n_cores, type="FORK")
    }

    on.exit(parallel::stopCluster(clust))
    perms <- parallel::parLapply(clust, 1:n_iter, function (i){
      random_geneset <- sample(expressed_genes_receiver, size = length(geneset_oi))
      predict_ligand_activities(geneset = random_geneset, background_expressed_genes = background_expressed_genes,
                                ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    })

  }

  return(perms)

}


#' @title Calculate ligand p-values from the bootstrapped ligand activity analysis
#' @description \code{calculate_p_value_bootstrap} Calculate the p-value for each ligand from the bootstrapped ligand activity analysis
#' @usage calculate_p_value_bootstrap(bootstrap_results, ligand_activities, metric = "aupr_corrected")
#' @param bootstrap_results Output of \code{\link{bootstrap_ligand_activity_analysis}}
#' @param ligand_activities Output of \code{\link{predict_ligand_activities}}
#' @param metric Metric to use (Default: "aupr_corrected")
#' @return Data frame with the ligand name, the number of times the bootstrapped value had a higher score than the observed (\code{total}), and the p-value for each ligand, calculated as \code{total/length(bootstrap_results)}
#' @examples
#' \dontrun{
#' permutations <- bootstrap_ligand_activity_analysis(expressed_genes_receiver, geneset_oi, background_expressed_genes,
#'                            ligand_target_matrix, potential_ligands, n_iter = 10)
#' p_values <- calculate_p_value_bootstrap(permutations, ligand_activities, metric = "aupr_corrected")
#' }
#' @export
calculate_p_value_bootstrap <- function(bootstrap_results, ligand_activities, metric = "aupr_corrected"){
  n_iter <- length(bootstrap_results)
  bootstrap_results %>% bind_rows(.id = "round") %>% group_by(test_ligand) %>%
    select(round, test_ligand, all_of(metric)) %>% arrange(test_ligand) %>% rename(null_metric = metric) %>%
    # Merge with observed values
    inner_join(ligand_activities %>% select(test_ligand, all_of(metric)) %>% rename(observed_metric = metric), by = "test_ligand") %>%
    # Calculate fraction of permutations with higher AUPR
    mutate(null_is_larger = null_metric > observed_metric) %>%
    summarise(total = sum(null_is_larger)) %>%
    mutate(pval = total/n_iter)

}


#' @title Cut off outer quantiles and rescale to a [0, 1] range
#'
#' @description \code{scale_quantile} Cut off outer quantiles and rescale to a [0, 1] range
#'
#' @usage
#' scale_quantile(x, outlier_cutoff = .05)
#'
#' @param x A numeric vector, matrix or data frame.
#' @param outlier_cutoff The quantile cutoff for outliers (default 0.05).
#'
#' @return The centered, scaled matrix or vector. The numeric centering and scalings used are returned as attributes.
#'
#' @examples
#' \dontrun{
#' ## Generate a matrix from a normal distribution
#' ## with a large standard deviation, centered at c(5, 5)
#' x <- matrix(rnorm(200*2, sd = 10, mean = 5), ncol = 2)
#'
#' ## Scale the dataset between [0,1]
#' x_scaled <- scale_quantile(x)
#'
#' ## Show ranges of each column
#' apply(x_scaled, 2, range)
#' }
#' @export
scale_quantile <- function(x, outlier_cutoff = .05) {
  # same function as scale_quantile from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the amazing (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)!
  if (is.null(dim(x))) {
    sc <- scale_quantile(matrix(x, ncol = 1), outlier_cutoff = outlier_cutoff)
    out <- sc[,1]
    names(out) <- names(x)
    attr(out, "addend") <- attr(sc, "addend")
    attr(out, "multiplier") <- attr(sc, "multiplier")
    out
  } else {
    quants <- apply(x, 2, stats::quantile, c(outlier_cutoff, 1 - outlier_cutoff), na.rm = TRUE)

    addend <- -quants[1,]
    divisor <- apply(quants, 2, diff)
    divisor[divisor == 0] <- 1

    apply_quantile_scale(x, addend, 1 / divisor)
  }
}
#' @title Prepare single-cell expression data to perform ligand activity analysis
#'
#' @description \code{convert_single_cell_expression_to_settings} Prepare single-cell expression data to perform ligand activity analysis
#'
#' @usage
#' convert_single_cell_expression_to_settings(cell_id, expression_matrix, setting_name, setting_from, regression = FALSE)
#'
#' @param cell_id Identity of the cell of interest
#' @param setting_name Name of the dataset
#' @param expression_matrix Gene expression matrix of single-cells
#' @param setting_from Character vector giving the gene symbols of the potentially active ligands you want to define ligand activities for.
#' @param regression Perform regression-based ligand activity analysis (TRUE) or classification-based ligand activity analysis (FALSE) by considering the genes expressed higher than the 0.975 quantiles as genes of interest. Default: FALSE.
#'
#' @return A list with slots $name, $from and $response respectively containing the setting name, potentially active ligands and the response to predict (whether genes belong to gene set of interest; i.e. most strongly expressed genes in a cell)
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' cell_ids = c("cell1","cell2")
#' expression_scaled = matrix(rnorm(length(genes)*2, sd = 0.5, mean = 0.5), nrow = 2)
#' rownames(expression_scaled) = cell_ids
#' colnames(expression_scaled) = genes
#' settings = convert_single_cell_expression_to_settings(cell_id = cell_ids[1], expression_matrix = expression_scaled, setting_name = "test", setting_from = potential_ligands)
#' }
#'
#' @export
#'
convert_single_cell_expression_to_settings = function(cell_id, expression_matrix, setting_name, setting_from, regression = FALSE){
  # input check
  requireNamespace("dplyr")

  if (regression == TRUE){
    response = expression_matrix[cell_id,]
  } else {
    response_continuous = expression_matrix[cell_id,]
    response = response_continuous >= quantile(response_continuous,0.975)
  }
  return(list(name = paste0(setting_name,"_",cell_id), from = setting_from, response = response))
}
#' @title Single-cell ligand activity prediction
#'
#' @description \code{predict_single_cell_ligand_activities} For every individual cell of interest, predict activities of ligands in regulating expression of genes that are stronger expressed in that cell compared to other cells (0.975 quantile). Ligand activities are defined as how well they predict the observed transcriptional response (i.e. gene set) according to the NicheNet model.
#'
#' @usage
#' predict_single_cell_ligand_activities(cell_ids, expression_scaled,ligand_target_matrix, potential_ligands, single = TRUE,...)
#'
#' @param cell_ids Identities of cells for which the ligand activities should be calculated.
#' @param expression_scaled Scaled expression matrix of single-cells (scaled such that high values indicate that a gene is stronger expressed in that cell compared to others)
#' @param ligand_target_matrix The NicheNet ligand-target matrix denoting regulatory potential scores between ligands and targets (ligands in columns).
#' @param potential_ligands Character vector giving the gene symbols of the potentially active ligands you want to define ligand activities for.
#' @param single TRUE if you want to calculate ligand activity scores by considering every ligand individually (recommended). FALSE if you want to calculate ligand activity scores as variable importances of a multi-ligand classification model.
#' @param ... Additional parameters for get_multi_ligand_importances if single = FALSE.
#'
#' @return A tibble giving several ligand activity scores for single cells. Following columns in the tibble: $setting, $test_ligand, $auroc, $aupr and $pearson.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' cell_ids = c("cell1","cell2")
#' expression_scaled = matrix(rnorm(length(genes)*2, sd = 0.5, mean = 0.5), nrow = 2)
#' rownames(expression_scaled) = cell_ids
#' colnames(expression_scaled) = genes
#' ligand_activities = predict_single_cell_ligand_activities(cell_ids = cell_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#' }
#'
#' @export
#'
predict_single_cell_ligand_activities = function(cell_ids, expression_scaled,ligand_target_matrix, potential_ligands, single = TRUE,...){
  settings_single_cell_ligand_pred = cell_ids %>% lapply(convert_single_cell_expression_to_settings, expression_scaled, "", potential_ligands)
  if (single == TRUE){
    settings_ligand_prediction = settings_single_cell_ligand_pred %>% convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = TRUE)

    ligand_importances = settings_ligand_prediction %>% lapply(get_single_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE) %>% bind_rows() %>% mutate(setting = gsub("^_","",setting))

  } else {
    settings_ligand_prediction = settings_single_cell_ligand_pred %>% convert_settings_ligand_prediction(all_ligands = potential_ligands, validation = FALSE, single = FALSE)

    ligand_importances = settings_ligand_prediction %>% lapply(get_multi_ligand_importances,ligand_target_matrix = ligand_target_matrix, known = FALSE, ...) %>% bind_rows() %>% mutate(setting = gsub("^_","",setting))

  }
  return(ligand_importances %>% select(setting,test_ligand,auroc,aupr,pearson))
}
#' @title Normalize single-cell ligand activities
#'
#' @description \code{normalize_single_cell_ligand_activities} Normalize single-cell ligand activities to make ligand activities over different cells comparable.
#' @usage
#' normalize_single_cell_ligand_activities(ligand_activities)
#'
#' @param ligand_activities Output from the function `predict_single_cell_ligand_activities`.
#'
#' @return A tibble giving the normalized ligand activity scores for single cells. Following columns in the tibble: $cell, $ligand, $pearson, which is the normalized ligand activity value.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' cell_ids = c("cell1","cell2")
#' expression_scaled = matrix(rnorm(length(genes)*2, sd = 0.5, mean = 0.5), nrow = 2)
#' rownames(expression_scaled) = cell_ids
#' colnames(expression_scaled) = genes
#' ligand_activities = predict_single_cell_ligand_activities(cell_ids = cell_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#' normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
#' }
#'
#' @export
#'
normalize_single_cell_ligand_activities = function(ligand_activities){
  single_ligand_activities_aupr_norm = ligand_activities %>%
    group_by(setting) %>%
    mutate(aupr = nichenetr::scaling_modified_zscore(aupr)) %>%
    ungroup() %>%
    rename(cell = setting, ligand = test_ligand) %>%
    distinct(cell,ligand,aupr)

  single_ligand_activities_aupr_norm_df = single_ligand_activities_aupr_norm %>%
    spread(cell, aupr,fill = min(.$aupr))

  single_ligand_activities_aupr_norm_matrix = single_ligand_activities_aupr_norm_df  %>%
    select(-ligand) %>%
    t() %>%
    magrittr::set_colnames(single_ligand_activities_aupr_norm_df$ligand)

  single_ligand_activities_aupr_norm_df = single_ligand_activities_aupr_norm_matrix %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    as_tibble()
}
#' @title Perform a correlation and regression analysis between cells' ligand activities and property scores of interest
#'
#' @description \code{single_ligand_activity_score_regression} Performs a correlation and regression analysis between cells' ligand activities and property scores of interest.
#' @usage
#' single_ligand_activity_score_regression(ligand_activities, scores_tbl)
#'
#' @param ligand_activities Output from the function `normalize_single_cell_ligand_activities`.
#' @param scores_tbl a tibble containing scores for every cell (columns: $cell and $score). The score should correspond to the property of interest
#'
#' @return A tibble giving for every ligand, the correlation/regression coefficients giving information about the relation between its activity and the property of interest.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' cell_ids = c("cell1","cell2")
#' expression_scaled = matrix(rnorm(length(genes)*2, sd = 0.5, mean = 0.5), nrow = 2)
#' rownames(expression_scaled) = cell_ids
#' colnames(expression_scaled) = genes
#' ligand_activities = predict_single_cell_ligand_activities(cell_ids = cell_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#' normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
#' cell_scores_tbl = tibble(cell = cell_ids, score = c(1,4))
#' regression_analysis_output = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)
#' }
#'
#' @export
#'
single_ligand_activity_score_regression = function(ligand_activities, scores_tbl){
  combined = inner_join(scores_tbl,ligand_activities)
  output = lapply(combined %>% select(-cell, -score), function(activity_prediction, combined){
    geneset_score = combined$score
    metrics = regression_evaluation(activity_prediction,geneset_score)
  }, combined)
  ligands = names(output)
  output_df = output %>% bind_rows() %>% mutate(ligand = ligands)
  return(output_df)
}
#' @title Assess how well cells' ligand activities predict a binary property of interest of cells.
#'
#' @description \code{single_ligand_activity_score_classification} Evaluates classification performances: it assesses how well cells' ligand activities can predict a binary property of interest.
#' @usage
#' single_ligand_activity_score_classification(ligand_activities, scores_tbl)
#'
#' @param ligand_activities Output from the function `normalize_single_cell_ligand_activities`.
#' @param scores_tbl a tibble indicating for every cell whether the property of interests holds TRUE or FALSE (columns: $cell: character vector with cell ids and $score: logical vector according to property of interest).
#'
#' @return A tibble giving for every ligand, the classification performance metrics giving information about the relation between its activity and the property of interest.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2","IL4")
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' potential_ligands = c("TNF","BMP2","IL4")
#' genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
#' cell_ids = c("cell1","cell2")
#' expression_scaled = matrix(rnorm(length(genes)*2, sd = 0.5, mean = 0.5), nrow = 2)
#' rownames(expression_scaled) = cell_ids
#' colnames(expression_scaled) = genes
#' ligand_activities = predict_single_cell_ligand_activities(cell_ids = cell_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#' normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
#' cell_scores_tbl = tibble(cell = cell_ids, score = c(TRUE,FALSE))
#' classification_analysis_output = single_ligand_activity_score_classification(normalized_ligand_activities,cell_scores_tbl)
#' }
#'
#' @export
#'
single_ligand_activity_score_classification = function(ligand_activities, scores_tbl){
  combined = inner_join(scores_tbl, ligand_activities)
  output = lapply(combined %>% select(-cell, -score), function(activity_prediction,
                                                               combined) {
    geneset_score = combined$score

    metrics = classification_evaluation_continuous_pred(activity_prediction,
                                                        geneset_score, iregulon = F)
  }, combined)


  ligands = names(output)
  output_df = output %>% bind_rows() %>% mutate(ligand = ligands)
  return(output_df)
}
single_ligand_activity_score_regression = function(ligand_activities, scores_tbl){
  combined = inner_join(scores_tbl,ligand_activities)
  output = lapply(combined %>% select(-cell, -score), function(activity_prediction, combined){
    geneset_score = combined$score
    metrics = regression_evaluation(activity_prediction,geneset_score)
  }, combined)
  ligands = names(output)
  output_df = output %>% bind_rows() %>% mutate(ligand = ligands)
  return(output_df)
}
#' @title Perform NicheNet analysis on Seurat object: explain DE between conditions
#'
#' @description \code{nichenet_seuratobj_aggregate} Perform NicheNet analysis on Seurat object: explain differential expression (DE) in a receiver celltype between two different conditions by ligands expressed by sender cells
#'
#' @param receiver Name of cluster identity/identities of cells that are presumably affected by intercellular communication with other cells
#' @param seurat_obj Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.
#' @param condition_colname Name of the column in the meta data dataframe that indicates which condition/sample cells were coming from.
#' @param condition_oi Condition of interest in which receiver cells were presumably affected by other cells. Should be a name present in the `condition_colname` column of the metadata.
#' @param condition_reference The second condition (e.g. reference or steady-state condition). Should be a name present in the `condition_colname` column of the metadata.
#' @param sender Determine the potential sender cells. Name of cluster identity/identities of cells that presumably affect expression in the receiver cell type. In case you want to look at all possible sender cell types in the data, you can  give this argument the value "all". "all" indicates thus that all cell types in the dataset will be considered as possible sender cells. As final option, you could give this argument the value "undefined"."undefined" won't look at ligands expressed by sender cells, but at all ligands for which a corresponding receptor is expressed. This could be useful if the presumably active sender cell is not profiled. Default: "all".
#' @param expression_pct To determine ligands and receptors expressed by sender and receiver cells, we consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Default: 0.10
#' @param lfc_cutoff Cutoff on log fold change in the wilcoxon differential expression test. Default: 0.25.
#' @param geneset Indicate whether to consider all DE genes between condition 1 and 2 ("DE"), or only genes upregulated in condition 1 ("up"), or only genes downregulad in condition 1 ("down").
#' @param filter_top_ligands Indicate whether output tables for ligand-target and ligand-receptor networks should be done for a filtered set of top ligands (TRUE) or for all ligands (FALSE). Default: TRUE.
#' @param top_n_ligands Indicate how many ligands should be extracted as top-ligands after ligand activity analysis. Only for these ligands, target genes and receptors will be returned. Default: 30.
#' @param top_n_targets To predict active, affected targets of the prioritized ligands, consider only DE genes if they also belong to the a priori top n ("top_n_targets") targets of a ligand. Default = 200.
#' @param cutoff_visualization Because almost no ligand-target scores have a regulatory potential score of 0, we clarify the heatmap visualization by giving the links with the lowest scores a score of 0. The cutoff_visualization paramter indicates this fraction of links that are given a score of zero. Default = 0.33.
#' @param ligand_target_matrix The NicheNet ligand-target matrix of the organism of interest denoting regulatory potential scores between ligands and targets (ligands in columns).
#' @param lr_network The ligand-receptor network (columns that should be present: $from, $to) of the organism of interest.
#' @param weighted_networks The NicheNet weighted networks of the organism of interest denoting interactions and their weights/confidences in the ligand-signaling and gene regulatory network.
#' @param verbose Print out the current analysis stage. Default: TRUE.
#' @param assay_oi The assay to be used for calculating expressed genes and the DE analysis. If NULL, the default assay of the Seurat object will be used.
#'
#' @return A list with the following elements:
#' $ligand_activities: data frame with output ligand activity analysis;
#' $top_ligands: top_n ligands based on ligand activity;
#' $top_targets: active, affected target genes of these ligands;
#' $top_receptors: receptors of these ligands;
#' $ligand_target_matrix: matrix indicating regulatory potential scores between active ligands and their predicted targets;
#' $ligand_target_heatmap: heatmap of ligand-target regulatory potential;
#' $ligand_target_df: data frame showing regulatory potential scores of predicted active ligand-target network;
#' $ligand_activity_target_heatmap: heatmap showing both ligand activity scores and target genes of these top ligands;
#' $ligand_expression_dotplot: expression dotplot of the top ligands;
#' $ligand_differential_expression_heatmap = differential expression heatmap of the top ligands;
#' $ligand_receptor_matrix: matrix of ligand-receptor interactions;
#' $ligand_receptor_heatmap: heatmap showing ligand-receptor interactions;
#' $ligand_receptor_df: data frame of ligand-receptor interactions;
#' $geneset_oi: a vector containing the set of genes used as input for the ligand activity analysis;
#' $background_expressed_genes: the background of genes to which the geneset will be compared in the ligand activity analysis.
#'
#' @import Seurat
#' @import dplyr
#' @importFrom magrittr set_rownames set_colnames
#'
#' @examples
#' \dontrun{
#' seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
#' lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#' weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#' nichenet_seuratobj_aggregate(receiver = "CD8 T", seurat_obj = seuratObj, condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "Mono", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
#' }
#'
#' @export
#'
nichenet_seuratobj_aggregate = function(receiver, seurat_obj, condition_colname, condition_oi, condition_reference, sender = "all",ligand_target_matrix,lr_network,weighted_networks,
                                        assay_oi = NULL, expression_pct = 0.10, lfc_cutoff = 0.25, geneset = "DE", filter_top_ligands = TRUE ,top_n_ligands = 30,
                                        top_n_targets = 200, cutoff_visualization = 0.33,
                                        verbose = TRUE)
{
  requireNamespace("Seurat")
  requireNamespace("dplyr")

  if (is.null(assay_oi)){
    assay_oi <- DefaultAssay(seurat_obj)
  } else {
    DefaultAssay(seurat_obj) <- assay_oi
  }

  if (verbose){
    print(paste0("The ", assay_oi, " assay will be used for the analysis."))
  }

  obj_version <- as.numeric(substr(seurat_obj@version, 1, 1))

  if (assay_oi == "Spatial") {
    warning("You are going to apply NicheNet on a spatial seurat object. Be sure it's ok to use NicheNet the way you are planning to do it. So this means: you should have changes in gene expression in receiver cells caused by cell-cell interactions. Note that in the case of spatial transcriptomics, you are not dealing with single cells but with 'spots' containing multiple cells of the same of different cell types.")
  }

  if(assay_oi == "integrated"){
    warning("The used assay is a result of the Seurat integration workflow. Make sure that the way of defining expressed and differentially expressed genes in this wrapper is appropriate for your integrated data.")
  }

  # Input check
  # Version 5
  if (obj_version >= 5){
    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, layer = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, layer = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }


  } else if (obj_version < 5) {

    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, slot = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, slot = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }
  }


  if(!condition_colname %in% colnames(seurat_obj@meta.data))
    stop("Your column indicating the conditions/samples of interest should be in the metadata dataframe")
  if(sum(condition_oi %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_oi))
    stop("condition_oi should be in the condition-indicating column")
  if(sum(condition_reference %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_reference))
    stop("condition_reference should be in the condition-indicating column")
  if(sum(receiver %in% unique(Idents(seurat_obj))) != length(receiver))
    stop("The defined receiver cell type should be an identity class of your seurat object")
  if(length(sender) == 1){
    if(sender != "all" & sender != "undefined"){
      if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
        stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
      }
    }
  } else {
    if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
      stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
    }
  }
  if(geneset != "DE" & geneset != "up" & geneset != "down")
    stop("geneset should be 'DE', 'up' or 'down'")

  # Read in and process NicheNet networks, define ligands and receptors
  if (verbose == TRUE){print("Read in and process NicheNet's networks")}
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()
  if (verbose == TRUE){print("Define expressed ligands and receptors in receiver and sender cells")}

  # step1 nichenet analysis: get expressed genes in sender and receiver cells

  ## receiver
  list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver) = receiver %>% unique()
  expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

  ## sender
  if (length(sender) == 1){
    if (sender == "all"){
      sender_celltypes = Idents(seurat_obj) %>% levels()
      list_expressed_genes_sender = sender_celltypes %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

    } else if (sender == "undefined") {
      expressed_genes_sender = union(Seurat::GetAssayData(seurat_obj, assay=assay_oi) %>% rownames(),rownames(ligand_target_matrix)) %>% union(colnames(ligand_target_matrix))
    } else if (sender != "all" & sender != "undefined") {
      sender_celltypes = sender
      list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes %>% unique()
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
    }
  } else {
    sender_celltypes = sender
    list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
    names(list_expressed_genes_sender) = sender_celltypes %>% unique()
    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }

  # step2 nichenet analysis: define background and gene list of interest: here differential expression between two conditions of cell type of interest
  if (verbose == TRUE){print("Perform DE analysis in receiver cell")}

  seurat_obj_receiver= subset(seurat_obj, idents = receiver)
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[[condition_colname, drop=TRUE]])
  DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = expression_pct, assay = assay_oi) %>% rownames_to_column("gene")

  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_receiver)

  if(SeuratV4 == TRUE){
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC <= lfc_cutoff) %>% pull(gene)
    }
  } else {
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC <= lfc_cutoff) %>% pull(gene)
    }
  }


  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0){
    stop("No genes were differentially expressed")
  }
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  # step3 nichenet analysis: define potential ligands
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  if (length(expressed_ligands) == 0){
    stop("No ligands expressed in sender cell")
  }
  if (length(expressed_receptors) == 0){
    stop("No receptors expressed in receiver cell")
  }
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  if (length(potential_ligands) == 0){
    stop("No potentially active ligands")
  }


  if (verbose == TRUE){print("Perform NicheNet ligand activity analysis")}

  # step4 perform NicheNet's ligand activity analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))

  if(filter_top_ligands == TRUE){
    best_upstream_ligands = ligand_activities %>% top_n(top_n_ligands, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  } else {
    best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  }

  if (verbose == TRUE){print("Infer active target genes of the prioritized ligands")}

  # step5 infer target genes of the top-ranked ligands
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets) %>% bind_rows() %>% drop_na()
  if(nrow(active_ligand_target_links_df) > 0){
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff_visualization)
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

    order_targets = order_targets %>% intersect(rownames(active_ligand_target_links))
    order_ligands = order_ligands %>% intersect(colnames(active_ligand_target_links))

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands,drop=FALSE] %>% t()
    p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) #+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
  } else {
    vis_ligand_target = NULL
    p_ligand_target_network = NULL
    print("no highly likely active targets found for top ligands")
  }
  # combined heatmap: overlay ligand activities
  ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
  p_ligand_aupr

  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_aupr)+10, ncol(vis_ligand_target)))
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h")

  combined_plot = cowplot::plot_grid(figures_without_legend,
                                     legends,
                                     rel_heights = c(10,2), nrow = 2, align = "hv")

  # ligand-receptor plot
  # get the ligand-receptor network of the top-ranked ligands
  if (verbose == TRUE){print("Infer receptors of the prioritized ligands")}

  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  if (nrow(lr_network_top_matrix) > 1){
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
  } else {
    order_receptors = rownames(lr_network_top_matrix)
  }
  if (ncol(lr_network_top_matrix) > 1) {
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  } else {
    order_ligands_receptor = colnames(lr_network_top_matrix)
  }

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  dim(vis_ligand_receptor_network) = c(length(order_receptors), length(order_ligands_receptor))
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  # DE analysis for each sender cell type -- of course only possible when having sender cell types
  if (length(sender) > 1){
    are_there_senders = TRUE
  }
  if(length(sender) == 1){
    if(sender != "undefined"){
      are_there_senders = TRUE
    } else {
      are_there_senders = FALSE
    }
  }

  if (are_there_senders == TRUE){
    if (verbose == TRUE){print("Perform DE analysis in sender cells")}
    seurat_obj = subset(seurat_obj, features = potential_ligands)

    DE_table_all = Idents(seurat_obj) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = seurat_obj, condition_colname = condition_colname, condition_oi = condition_oi, condition_reference = condition_reference, expression_pct = expression_pct, celltype_col = NULL) %>% reduce(full_join, by = "gene") # use this if cell type labels are the identities of your Seurat object -- if not: indicate the celltype_col properly
    DE_table_all[is.na(DE_table_all)] = 0

    # Combine ligand activities with DE information
    ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene), by = "ligand")
    ligand_activities_de[is.na(ligand_activities_de)] = 0

    # make LFC heatmap
    lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
    rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

    order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
    vis_ligand_lfc = lfc_matrix[order_ligands,]
    vis_ligand_lfc = vis_ligand_lfc %>% as.matrix(ncol = length(Idents(seurat_obj) %>% levels() %>% intersect(sender_celltypes)))
    colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

    p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))

    # ligand expression Seurat dotplot
    real_makenames_conversion = lr_network$from %>% unique() %>% magrittr::set_names(lr_network$from %>% unique() %>% make.names())
    order_ligands_adapted = real_makenames_conversion[order_ligands]
    names(order_ligands_adapted) = NULL

    seurat_obj_subset = seurat_obj %>% subset(idents = sender_celltypes)
    seurat_obj_subset = SetIdent(seurat_obj_subset, value = seurat_obj_subset[[condition_colname, drop=TRUE]]) %>% subset(idents = condition_oi) ## only shows cells of the condition of interest
    rotated_dotplot = DotPlot(seurat_obj %>% subset(cells = Cells(seurat_obj_subset)), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
    rm(seurat_obj_subset)

    # combined plot
    figures_without_legend = cowplot::plot_grid(
      p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
      rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
      p_ligand_lfc + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()) + ylab(""),
      p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
      align = "hv",
      nrow = 1,
      rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc) + 7, ncol(vis_ligand_lfc) + 8, ncol(vis_ligand_target)))

    legends = cowplot::plot_grid(
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
      ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_lfc)),
      ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
      nrow = 1,
      align = "h", rel_widths = c(1.5, 1, 1, 1))

    combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
    combined_plot

  } else {
    rotated_dotplot = NULL
    p_ligand_lfc = NULL
  }

  return(list(
    ligand_activities = ligand_activities,
    top_ligands = best_upstream_ligands,
    top_targets = active_ligand_target_links_df$target %>% unique(),
    top_receptors = lr_network_top_df_large$to %>% unique(),
    ligand_target_matrix = vis_ligand_target,
    ligand_target_heatmap = p_ligand_target_network,
    ligand_target_df = active_ligand_target_links_df,
    ligand_expression_dotplot = rotated_dotplot,
    ligand_differential_expression_heatmap = p_ligand_lfc,
    ligand_activity_target_heatmap = combined_plot,
    ligand_receptor_matrix = vis_ligand_receptor_network,
    ligand_receptor_heatmap = p_ligand_receptor_network,
    ligand_receptor_df = lr_network_top_df_large %>% rename(ligand = from, receptor = to),
    geneset_oi = geneset_oi,
    background_expressed_genes = background_expressed_genes
  ))
}
#' @title Determine expressed genes of a cell type from an input object
#' @description Return the genes that are expressed in given cell cluster(s) based on the fraction of cells in the cluster(s) that should express the cell.
#' @param celltype_oi Character vector of cell type(s) to be considered. If input is a Seurat object, these must correspond to the cell identities from \code{Idents}.
#' @param object Input matrix with rows as genes and columns as cells
#' @param pct We consider genes expressed if they are expressed in at least a specific fraction of cells of the given cluster(s). This number indicates this fraction. Default: 0.10. Choice of this parameter is important and depends largely on the used sequencing platform. We recommend to require a lower fraction (like the default 0.10) for 10X data than for e.g. Smart-seq2 data.
#' @return A vector containing gene symbols of the expressed genes
#'
#' @rdname get_expressed_genes
#' @export
get_expressed_genes <- function(celltype_oi, object, ...) {
  UseMethod("get_expressed_genes", object)
}

#' @param celltype_annot Vector of cell type annotations
#' @import dplyr
#' @rdname get_expressed_genes
#' @examples
#' \dontrun{
#' # For sparse matrix
#' get_expressed_genes("CD8 T", GetAssayData(seuratObj), seuratObj$celltype, pct = 0.10)
#' }
#'
#' @export
get_expressed_genes.default <- function(celltype_oi, object, celltype_annot, pct = 0.1) {
  requireNamespace("dplyr")

  # Check that length of metadata is equal to number of cells
  if (length(celltype_annot) != ncol(object)) {
    stop("Length of metadata is not equal to number of cells")
  }

  # Check that celltype_oi is in metadata
  if (!celltype_oi %in% celltype_annot) {
    stop("Cell type of interest is not in celltype_annot")
  }

  exprs_mat <- object %>% .[, celltype_annot == celltype_oi]

  # Get cells of interest in matrix
  n_cells_oi_in_matrix <- ncol(exprs_mat)

  if (n_cells_oi_in_matrix < 5000) {
    genes <- exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi_in_matrix
    }) %>% .[. >= pct] %>% names()
  }
  else {
    splits <- split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes <- splits %>% lapply(function(genes_indices) {
      begin_i <- genes_indices[1]
      end_i <- genes_indices[length(genes_indices)]
      exprs <- exprs_mat[begin_i:end_i, , drop = FALSE]
      genes <- exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi_in_matrix
      }) %>% .[. >= pct] %>% names()
    }) %>% unlist() %>%
      unname()
  }
  return(genes)

}

#' @param seurat_obj Single-cell expression or spatial dataset as Seurat object
#' @param assay_oi If wanted: specify yourself which assay to look for. If not NULL, the \code{DefaultAssay} of the Seurat object is used.
#' @param ... additional parameters passed to \code{GetAssayData} (in case the slot/layer needs to be specified)
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # For Seurat object
#' get_expressed_genes(celltype_oi = "CD8 T", seurat_obj = seuratObj, pct = 0.10)
#' }
#'
#' @rdname get_expressed_genes
#' @export
#'
get_expressed_genes.Seurat = function(celltype_oi, seurat_obj, pct = 0.1, assay_oi = NULL, ...){
  requireNamespace("Seurat")
  requireNamespace("dplyr")

  # If assay_oi is not specified, use DefaultAssay
  if (is.null(assay_oi)){
    assay_oi <- DefaultAssay(seurat_obj)
  } else {
    # Check that assay exists in Seurat object
    if (!assay_oi %in% Seurat::Assays(seurat_obj)){
      stop("assay_oi not found in Seurat object")
    }
    # Set assay_oi to DefaultAssay
    DefaultAssay(seurat_obj) <- assay_oi

  }

  # If assay_oi is Spatial, give warning
  if (assay_oi == "Spatial"){
    warning("The Spatial assay will be used to define expressed gene. If the spatial data is spot-based (mixture of cells) and not single-cell resolution, we recommend against directly using nichenetr on because you want to look at cell-cell interactions, and not at spot-spot interactions! ;-) ")
  }

  # Check that celltype_oi is in Seurat object
  if (!all(celltype_oi %in% unique(Idents(seurat_obj)))) {
    stop("One or more provided cell clusters is not part of the 'Idents' of your Seurat object")
  }

  # Get cell identities of cluster of interest
  cells_oi <- Idents(seurat_obj) %>% .[. %in% celltype_oi] %>% names()

  # Get expression matrix from assay_oi
  #cells_oi_in_matrix <- intersect(colnames(seurat_obj[[assay_oi]]@data), cells_oi)
  #exprs_mat = seurat_obj[[assay_oi]]@data %>% .[, cells_oi_in_matrix]

  exprs_mat <- subset(seurat_obj, idents = celltype_oi) %>%
    GetAssayData(assay = assay_oi, ...)

  if (length(cells_oi) != ncol(exprs_mat)){
    warning("Not all cells of interest are in your expression matrix.")
  }

  # Use defined cells and exprs matrix to get expressed genes
  n_cells_oi_in_matrix <-  ncol(exprs_mat)
  if (n_cells_oi_in_matrix < 5000) {
    genes <- exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi_in_matrix
    }) %>% .[. >= pct] %>% names()
  }
  else {
    splits <- split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes <- splits %>% lapply(function(genes_indices) {
      begin_i <- genes_indices[1]
      end_i <- genes_indices[length(genes_indices)]
      exprs <- exprs_mat[begin_i:end_i, , drop = FALSE]
      genes <- exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi_in_matrix
      }) %>% .[. >= pct] %>% names()
    }) %>% unlist() %>%
      unname()
  }
  return(genes)
}
#' @title Perform NicheNet analysis on Seurat object: explain DE between two cell clusters
#'
#' @description \code{nichenet_seuratobj_cluster_de} Perform NicheNet analysis on Seurat object: explain differential expression (DE) between two 'receiver' cell clusters by ligands expressed by neighboring cells.
#' @usage
#' nichenet_seuratobj_cluster_de(seurat_obj, receiver_affected, receiver_reference, sender = "all",ligand_target_matrix,lr_network,weighted_networks,expression_pct = 0.10, lfc_cutoff = 0.25, geneset = "DE", filter_top_ligands = TRUE, top_n_ligands = 30,top_n_targets = 200, cutoff_visualization = 0.33,verbose = TRUE, assay_oi = NULL)
#'
#' @param seurat_obj Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.
#' @param receiver_reference Name of cluster identity/identities of "steady-state" cells, before they are affected by intercellular communication with other cells
#' @param receiver_affected Name of cluster identity/identities of "affected" cells that were presumably affected by intercellular communication with other cells
#' @param sender Determine the potential sender cells. Name of cluster identity/identities of cells that presumably affect expression in the receiver cell type. In case you want to look at all possible sender cell types in the data, you can  give this argument the value "all". "all" indicates thus that all cell types in the dataset will be considered as possible sender cells. As final option, you could give this argument the value "undefined"."undefined" won't look at ligands expressed by sender cells, but at all ligands for which a corresponding receptor is expressed. This could be useful if the presumably active sender cell is not profiled. Default: "all".
#' @param expression_pct To determine ligands and receptors expressed by sender and receiver cells, we consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Default: 0.10
#' @param lfc_cutoff Cutoff on log fold change in the wilcoxon differential expression test. Default: 0.25.
#' @param geneset Indicate whether to consider all DE genes between condition 1 and 2 ("DE"), or only genes upregulated in condition 1 ("up"), or only genes downregulad in condition 1 ("down").
#' @param filter_top_ligands Indicate whether output tables for ligand-target and ligand-receptor networks should be done for a filtered set of top ligands (TRUE) or for all ligands (FALSE). Default: TRUE.
#' @param top_n_ligands Indicate how many ligands should be extracted as top-ligands after ligand activity analysis. Only for these ligands, target genes and receptors will be returned. Default: 30.
#' @param top_n_targets To predict active, affected targets of the prioritized ligands, consider only DE genes if they also belong to the a priori top n ("top_n_targets") targets of a ligand. Default = 200.
#' @param cutoff_visualization Because almost no ligand-target scores have a regulatory potential score of 0, we clarify the heatmap visualization by giving the links with the lowest scores a score of 0. The cutoff_visualization paramter indicates this fraction of links that are given a score of zero. Default = 0.33.
#' @param ligand_target_matrix The NicheNet ligand-target matrix denoting regulatory potential scores between ligands and targets (ligands in columns).
#' @param lr_network The ligand-receptor network (columns that should be present: $from, $to).
#' @param weighted_networks The NicheNet weighted networks denoting interactions and their weights/confidences in the ligand-signaling and gene regulatory network.
#' @param verbose Print out the current analysis stage. Default: TRUE.
#' @param assay_oi The assay to be used for calculating expressed genes and the DE analysis. If NULL, the default assay of the Seurat object will be used.
#'
#' @return A list with the following elements:
#' $ligand_activities: data frame with output ligand activity analysis;
#' $top_ligands: top_n ligands based on ligand activity;
#' $top_targets: active, affected target genes of these ligands;
#' $top_receptors: receptors of these ligands;
#' $ligand_target_matrix: matrix indicating regulatory potential scores between active ligands and their predicted targets;
#' $ligand_target_heatmap: heatmap of ligand-target regulatory potential;
#' $ligand_target_df: data frame showing regulatory potential scores of predicted active ligand-target network;
#' $ligand_activity_target_heatmap: heatmap showing both ligand activity scores and target genes of these top ligands;
#' $ligand_expression_dotplot: expression dotplot of the top ligands;
#' $ligand_receptor_matrix: matrix of ligand-receptor interactions;
#' $ligand_receptor_heatmap: heatmap showing ligand-receptor interactions;
#' $ligand_receptor_df: data frame of ligand-receptor interactions;
#' $geneset_oi: a vector containing the set of genes used as input for the ligand activity analysis;
#' $background_expressed_genes: the background of genes to which the geneset will be compared in the ligand activity analysis.
#'
#' @import Seurat
#' @import dplyr
#' @importFrom magrittr set_rownames set_colnames
#'
#' @examples
#' \dontrun{
#' seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
#' lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#' weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#' # works, but does not make sense
#' nichenet_seuratobj_cluster_de(seurat_obj = seuratObj, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = "Mono", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
#' # type of analysis for which this would make sense
#' nichenet_seuratobj_cluster_de(seurat_obj = seuratObj, receiver_affected = "p-EMT-pos-cancer", receiver_reference = "p-EMT-neg-cancer", sender = "Fibroblast", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
#' }
#'
#' @export
#'
nichenet_seuratobj_cluster_de = function(seurat_obj, receiver_affected, receiver_reference, sender = "all",ligand_target_matrix,lr_network,weighted_networks,
                                         assay_oi = NULL, expression_pct = 0.10, lfc_cutoff = 0.25, geneset = "DE", filter_top_ligands = TRUE, top_n_ligands = 30,
                                        top_n_targets = 200, cutoff_visualization = 0.33,
                                        verbose = TRUE)
{
  requireNamespace("Seurat")
  requireNamespace("dplyr")

  if (is.null(assay_oi)){
    assay_oi <- DefaultAssay(seurat_obj)
  } else {
    DefaultAssay(seurat_obj) <- assay_oi
  }

  if (verbose){
    print(paste0("The ", assay_oi, " assay will be used for the analysis."))
  }

  obj_version <- as.numeric(substr(seurat_obj@version, 1, 1))

  if (assay_oi == "Spatial") {
    warning("You are going to apply NicheNet on a spatial seurat object. Be sure it's ok to use NicheNet the way you are planning to do it. So this means: you should have changes in gene expression in receiver cells caused by cell-cell interactions. Note that in the case of spatial transcriptomics, you are not dealing with single cells but with 'spots' containing multiple cells of the same of different cell types.")
  }

  if(assay_oi == "integrated"){
    warning("The used assay is a result of the Seurat integration workflow. Make sure that the way of defining expressed and differentially expressed genes in this wrapper is appropriate for your integrated data.")
  }

  # Input check
  # Version 5
  if (obj_version >= 5){
    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, layer = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, layer = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }


  } else if (obj_version < 5) {

    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, slot = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, slot = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }
  }


  if(sum(receiver_affected %in% unique(Idents(seurat_obj))) != length(receiver_affected))
    stop("The defined receiver_affected cell type should be an identity class of your seurat object")
  if(sum(receiver_reference %in% unique(Idents(seurat_obj))) != length(receiver_reference))
    stop("The defined receiver_reference cell type should be an identity class of your seurat object")
  if(length(sender) == 1){
    if(sender != "all" & sender != "undefined"){
      if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
        stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
      }
    }
  } else {
    if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
      stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
    }
  }
  if(geneset != "DE" & geneset != "up" & geneset != "down")
    stop("geneset should be 'DE', 'up' or 'down'")

  # Read in and process NicheNet networks, define ligands and receptors
  if (verbose == TRUE){print("Read in and process NicheNet's networks")}
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  if (verbose == TRUE){print("Define expressed ligands and receptors in receiver and sender cells")}

  # step1 nichenet analysis: get expressed genes in sender and receiver cells

  ## receiver
  # expressed genes: only in steady state population (for determining receptors)
  list_expressed_genes_receiver_ss = c(receiver_reference) %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver_ss) = c(receiver_reference) %>% unique()
  expressed_genes_receiver_ss = list_expressed_genes_receiver_ss %>% unlist() %>% unique()

  # expressed genes: both in steady state and affected population (for determining background of expressed genes)
  list_expressed_genes_receiver = c(receiver_reference,receiver_affected) %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver) = c(receiver_reference,receiver_affected) %>% unique()
  expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

  ## sender
  if (length(sender) == 1){
    if (sender == "all"){
      sender_celltypes = Idents(seurat_obj) %>% levels()
      list_expressed_genes_sender = sender_celltypes %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

    } else if (sender == "undefined") {
      expressed_genes_sender = union(Seurat::GetAssayData(seurat_obj, assay=assay_oi) %>% rownames(),rownames(ligand_target_matrix)) %>% union(colnames(ligand_target_matrix))
    } else if (sender != "all" & sender != "undefined") {
      sender_celltypes = sender
      list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes %>% unique()
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
    }
  } else {
    sender_celltypes = sender
    list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
    names(list_expressed_genes_sender) = sender_celltypes %>% unique()
    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }

  # step2 nichenet analysis: define background and gene list of interest: here differential expression between two conditions of cell type of interest
  if (verbose == TRUE){print("Perform DE analysis between two receiver cell clusters")}

  DE_table_receiver = FindMarkers(object = seurat_obj, ident.1 = receiver_affected, ident.2 = receiver_reference, min.pct = expression_pct, assay = assay_oi) %>% rownames_to_column("gene")

  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_receiver)

  if(SeuratV4 == TRUE){
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC <= lfc_cutoff) %>% pull(gene)
    }
  } else {
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC <= lfc_cutoff) %>% pull(gene)
    }
  }



  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0){
    stop("No genes were differentially expressed")
  }
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  # step3 nichenet analysis: define potential ligands
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  if (length(expressed_ligands) == 0){
    stop("No ligands expressed in sender cell")
  }
  if (length(expressed_receptors) == 0){
    stop("No receptors expressed in receiver cell")
  }
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  if (length(potential_ligands) == 0){
    stop("No potentially active ligands")
  }

  if (verbose == TRUE){print("Perform NicheNet ligand activity analysis")}

  # step4 perform NicheNet's ligand activity analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))

  if(filter_top_ligands == TRUE){
    best_upstream_ligands = ligand_activities %>% top_n(top_n_ligands, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  } else {
    best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  }
  if (verbose == TRUE){print("Infer active target genes of the prioritized ligands")}

  # step5 infer target genes of the top-ranked ligands
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets) %>% bind_rows() %>% drop_na()

  if(nrow(active_ligand_target_links_df) > 0){
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff_visualization)
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

    order_targets = order_targets %>% intersect(rownames(active_ligand_target_links))
    order_ligands = order_ligands %>% intersect(colnames(active_ligand_target_links))

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
    p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) #+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
  } else {
    vis_ligand_target = NULL
    p_ligand_target_network = NULL
    print("no highly likely active targets found for top ligands")
  }

  # combined heatmap: overlay ligand activities
  ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
  p_ligand_aupr

  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_aupr)+10, ncol(vis_ligand_target)))
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h")

  combined_plot = cowplot::plot_grid(figures_without_legend,
                                     legends,
                                     rel_heights = c(10,2), nrow = 2, align = "hv")

  # ligand-receptor plot
  # get the ligand-receptor network of the top-ranked ligands
  if (verbose == TRUE){print("Infer receptors of the prioritized ligands")}

  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  if (nrow(lr_network_top_matrix) > 1){
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
  } else {
    order_receptors = rownames(lr_network_top_matrix)
  }
  if (ncol(lr_network_top_matrix) > 1) {
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  } else {
    order_ligands_receptor = colnames(lr_network_top_matrix)
  }

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  dim(vis_ligand_receptor_network) = c(length(order_receptors), length(order_ligands_receptor))

  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  # ligand expression Seurat dotplot
  if (length(sender) > 1){
    are_there_senders = TRUE
  }
  if(length(sender) == 1){
    if(sender != "undefined"){
      are_there_senders = TRUE
    } else {
      are_there_senders = FALSE
    }
  }

  if (are_there_senders == TRUE){
    real_makenames_conversion = lr_network$from %>% unique() %>% magrittr::set_names(lr_network$from %>% unique() %>% make.names())
    order_ligands_adapted = real_makenames_conversion[order_ligands]
    names(order_ligands_adapted) = NULL
    rotated_dotplot = DotPlot(seurat_obj %>% subset(idents = sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

  } else {
    rotated_dotplot = NULL
  }


  return(list(
    ligand_activities = ligand_activities,
    top_ligands = best_upstream_ligands,
    top_targets = active_ligand_target_links_df$target %>% unique(),
    top_receptors = lr_network_top_df_large$to %>% unique(),
    ligand_target_matrix = vis_ligand_target,
    ligand_target_heatmap = p_ligand_target_network,
    ligand_target_df = active_ligand_target_links_df,
    ligand_expression_dotplot = rotated_dotplot,
    ligand_activity_target_heatmap = combined_plot,
    ligand_receptor_matrix = vis_ligand_receptor_network,
    ligand_receptor_heatmap = p_ligand_receptor_network,
    ligand_receptor_df = lr_network_top_df_large %>% rename(ligand = from, receptor = to),
    geneset_oi = geneset_oi,
    background_expressed_genes = background_expressed_genes

  ))
}
#' @title Perform NicheNet analysis on Seurat object: explain DE between two cell clusters from separate conditions
#'
#' @description \code{nichenet_seuratobj_aggregate_cluster_de} Perform NicheNet analysis on Seurat object: explain differential expression (DE) between two 'receiver' cell clusters coming from different conditions, by ligands expressed by neighboring cells.
#' @usage
#' nichenet_seuratobj_aggregate_cluster_de(seurat_obj, receiver_affected, receiver_reference, condition_colname, condition_oi, condition_reference, sender = "all",ligand_target_matrix,lr_network,weighted_networks,expression_pct = 0.10, lfc_cutoff = 0.25, geneset = "DE", filter_top_ligands = TRUE, top_n_ligands = 30,top_n_targets = 200, cutoff_visualization = 0.33,verbose = TRUE, assay_oi = NULL)
#'
#' @param seurat_obj Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.
#' @param receiver_reference Name of cluster identity/identities of "steady-state" cells, before they are affected by intercellular communication with other cells
#' @param receiver_affected Name of cluster identity/identities of "affected" cells that were presumably affected by intercellular communication with other cells
#' @param condition_colname Name of the column in the meta data dataframe that indicates which condition/sample cells were coming from.
#' @param condition_oi Condition of interest in which receiver cells were presumably affected by other cells. Should be a name present in the `condition_colname` column of the metadata.
#' @param condition_reference The second condition (e.g. reference or steady-state condition). Should be a name present in the `condition_colname` column of the metadata.
#' @param sender Determine the potential sender cells. Name of cluster identity/identities of cells that presumably affect expression in the receiver cell type. In case you want to look at all possible sender cell types in the data, you can  give this argument the value "all". "all" indicates thus that all cell types in the dataset will be considered as possible sender cells. As final option, you could give this argument the value "undefined"."undefined" won't look at ligands expressed by sender cells, but at all ligands for which a corresponding receptor is expressed. This could be useful if the presumably active sender cell is not profiled. Default: "all".
#' @param expression_pct To determine ligands and receptors expressed by sender and receiver cells, we consider genes expressed if they are expressed in at least a specific fraction of cells of a cluster. This number indicates this fraction. Default: 0.10
#' @param lfc_cutoff Cutoff on log fold change in the wilcoxon differential expression test. Default: 0.25.
#' @param geneset Indicate whether to consider all DE genes between condition 1 and 2 ("DE"), or only genes upregulated in condition 1 ("up"), or only genes downregulad in condition 1 ("down").
#' @param filter_top_ligands Indicate whether output tables for ligand-target and ligand-receptor networks should be done for a filtered set of top ligands (TRUE) or for all ligands (FALSE). Default: TRUE.
#' @param top_n_ligands Indicate how many ligands should be extracted as top-ligands after ligand activity analysis. Only for these ligands, target genes and receptors will be returned. Default: 30.
#' @param top_n_targets To predict active, affected targets of the prioritized ligands, consider only DE genes if they also belong to the a priori top n ("top_n_targets") targets of a ligand. Default = 200.
#' @param cutoff_visualization Because almost no ligand-target scores have a regulatory potential score of 0, we clarify the heatmap visualization by giving the links with the lowest scores a score of 0. The cutoff_visualization paramter indicates this fraction of links that are given a score of zero. Default = 0.33.
#' @param ligand_target_matrix The NicheNet ligand-target matrix of the organism of interest denoting regulatory potential scores between ligands and targets (ligands in columns).
#' @param lr_network The ligand-receptor network (columns that should be present: $from, $to) of the organism of interest.
#' @param weighted_networks The NicheNet weighted networks of the organism of interest denoting interactions and their weights/confidences in the ligand-signaling and gene regulatory network.
#' @param verbose Print out the current analysis stage. Default: TRUE.
#' @param assay_oi The assay to be used for calculating expressed genes and the DE analysis. If NULL, the default assay of the Seurat object will be used.
#'
#' @return A list with the following elements:
#' $ligand_activities: data frame with output ligand activity analysis;
#' $top_ligands: top_n ligands based on ligand activity;
#' $top_targets: active, affected target genes of these ligands;
#' $top_receptors: receptors of these ligands;
#' $ligand_target_matrix: matrix indicating regulatory potential scores between active ligands and their predicted targets;
#' $ligand_target_heatmap: heatmap of ligand-target regulatory potential;
#' $ligand_target_df: data frame showing regulatory potential scores of predicted active ligand-target network;
#' $ligand_activity_target_heatmap: heatmap showing both ligand activity scores and target genes of these top ligands;
#' $ligand_expression_dotplot: expression dotplot of the top ligands;
#' $ligand_receptor_matrix: matrix of ligand-receptor interactions;
#' $ligand_receptor_heatmap: heatmap showing ligand-receptor interactions;
#' $ligand_receptor_df: data frame of ligand-receptor interactions;
#' $geneset_oi: a vector containing the set of genes used as input for the ligand activity analysis;
#' $background_expressed_genes: the background of genes to which the geneset will be compared in the ligand activity analysis.
#'
#' @import Seurat
#' @import dplyr
#' @importFrom magrittr set_rownames set_colnames
#'
#' @examples
#' \dontrun{
#' seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
#' lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
#' weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
#' nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seuratObj, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", sender = "Mono", ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks)
#' }
#'
#' @export
#'
nichenet_seuratobj_aggregate_cluster_de = function(seurat_obj, receiver_affected, receiver_reference,
                                         condition_colname, condition_oi, condition_reference, sender = "all",
                                         ligand_target_matrix,lr_network,weighted_networks, assay_oi = NULL,
                                         expression_pct = 0.10, lfc_cutoff = 0.25, geneset = "DE", filter_top_ligands = TRUE, top_n_ligands = 30,
                                         top_n_targets = 200, cutoff_visualization = 0.33,
                                         verbose = TRUE)
{

  requireNamespace("Seurat")
  requireNamespace("dplyr")

  if (is.null(assay_oi)){
    assay_oi <- DefaultAssay(seurat_obj)
  } else {
    DefaultAssay(seurat_obj) <- assay_oi
  }

  if (verbose){
    print(paste0("The ", assay_oi, " assay will be used for the analysis."))
  }

  obj_version <- as.numeric(substr(seurat_obj@version, 1, 1))

  if (assay_oi == "Spatial") {
    warning("You are going to apply NicheNet on a spatial seurat object. Be sure it's ok to use NicheNet the way you are planning to do it. So this means: you should have changes in gene expression in receiver cells caused by cell-cell interactions. Note that in the case of spatial transcriptomics, you are not dealing with single cells but with 'spots' containing multiple cells of the same of different cell types.")
  }

  if(assay_oi == "integrated"){
    warning("The used assay is a result of the Seurat integration workflow. Make sure that the way of defining expressed and differentially expressed genes in this wrapper is appropriate for your integrated data.")
  }

  # Input check
  # Version 5
  if (obj_version >= 5){
    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, layer = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, layer = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }


  } else if (obj_version < 5) {

    if (sum(dim(GetAssayData(seurat_obj, assay = assay_oi, slot = "data"))) == 0){
      stop("Seurat object should contain normalized expression data (numeric matrix). Check 'GetAssayData(seurat_obj, assay = assay_oi, layer = 'data')'")
    }

    if (!class(GetAssayData(seurat_obj, assay = assay_oi, slot = "data")) %in% c("matrix", "dgCMatrix")) {
      warning("The normalized expression data should be a matrix object.")
    }
  }


  if(sum(receiver_affected %in% unique(Idents(seurat_obj))) != length(receiver_affected))
    stop("The defined receiver_affected cell type should be an identity class of your seurat object")
  if(sum(receiver_reference %in% unique(Idents(seurat_obj))) != length(receiver_reference))
    stop("The defined receiver_reference cell type should be an identity class of your seurat object")
  if(!condition_colname %in% colnames(seurat_obj@meta.data))
    stop("Your column indicating the conditions/samples of interest should be in the metadata dataframe")
  if(sum(condition_oi %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_oi))
    stop("condition_oi should be in the condition-indicating column")
  if(sum(condition_reference %in% c(seurat_obj[[condition_colname]] %>% unlist() %>% as.character() %>% unique())) != length(condition_reference))
    stop("condition_reference should be in the condition-indicating column")
  if(length(sender) == 1){
    if(sender != "all" & sender != "undefined"){
      if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
        stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
      }
    }
  } else {
      if(sum(sender %in% unique(Idents(seurat_obj))) != length(sender)){
        stop("The sender argument should be 'all' or 'undefined' or an identity class of your seurat object")
      }
  }
  if(geneset != "DE" & geneset != "up" & geneset != "down")
    stop("geneset should be 'DE', 'up' or 'down'")

  # Read in and process NicheNet networks, define ligands and receptors
  if (verbose == TRUE){print("Read in and process NicheNet's networks")}
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

  ligands = lr_network %>% pull(from) %>% unique()
  receptors = lr_network %>% pull(to) %>% unique()

  if (verbose == TRUE){print("Define expressed ligands and receptors in receiver and sender cells")}

  # step1 nichenet analysis: get expressed genes in sender and receiver cells

  ## receiver
  # expressed genes: only in steady state population (for determining receptors)
  list_expressed_genes_receiver_ss = c(receiver_reference) %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver_ss) = c(receiver_reference) %>% unique()
  expressed_genes_receiver_ss = list_expressed_genes_receiver_ss %>% unlist() %>% unique()

  # expressed genes: both in steady state and affected population (for determining background of expressed genes)
  list_expressed_genes_receiver = c(receiver_reference,receiver_affected) %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
  names(list_expressed_genes_receiver) = c(receiver_reference,receiver_affected) %>% unique()
  expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

  ## sender
  if (length(sender) == 1){
    if (sender == "all"){
      sender_celltypes = Idents(seurat_obj) %>% levels()
      list_expressed_genes_sender = sender_celltypes %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

    } else if (sender == "undefined") {
      expressed_genes_sender = union(Seurat::GetAssayData(seurat_obj, assay=assay_oi) %>% rownames(),rownames(ligand_target_matrix)) %>% union(colnames(ligand_target_matrix))
    } else if (sender != "all" & sender != "undefined") {
      sender_celltypes = sender
      list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
      names(list_expressed_genes_sender) = sender_celltypes %>% unique()
      expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
    }
  } else {
    sender_celltypes = sender
    list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, expression_pct, assay_oi)
    names(list_expressed_genes_sender) = sender_celltypes %>% unique()
    expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }

  # step2 nichenet analysis: define background and gene list of interest: here differential expression between two conditions of cell type of interest
  if (verbose == TRUE){print("Perform DE analysis between two receiver cell clusters")}

  seurat_obj_receiver_affected= subset(seurat_obj, idents = receiver_affected)
  seurat_obj_receiver_affected = SetIdent(seurat_obj_receiver_affected, value = seurat_obj_receiver_affected[[condition_colname, drop=TRUE]])
  seurat_obj_receiver_affected= subset(seurat_obj_receiver_affected, idents = condition_oi)

  seurat_obj_receiver_reference= subset(seurat_obj, idents = receiver_reference)
  seurat_obj_receiver_reference = SetIdent(seurat_obj_receiver_reference, value = seurat_obj_receiver_reference[[condition_colname, drop=TRUE]])
  seurat_obj_receiver_reference= subset(seurat_obj_receiver_reference, idents = condition_reference)

  seurat_obj_receiver = merge(seurat_obj_receiver_affected, seurat_obj_receiver_reference)

  if (obj_version >= 5 & inherits(seurat_obj[[assay_oi]], "Assay5")){
    # If an object was updated instead of being created in v5, it will not inherit Assay5
    seurat_obj_receiver <- SeuratObject::JoinLayers(seurat_obj_receiver, assay=assay_oi)
  }

  DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = expression_pct, assay = assay_oi) %>% rownames_to_column("gene")


  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_receiver)

  if(SeuratV4 == TRUE){
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_log2FC <= lfc_cutoff) %>% pull(gene)
    }
  } else {
    if (geneset == "DE"){
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "up") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC >= lfc_cutoff) %>% pull(gene)
    } else if (geneset == "down") {
      geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & avg_logFC <= lfc_cutoff) %>% pull(gene)
    }
  }

  geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0){
    stop("No genes were differentially expressed")
  }
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

  # step3 nichenet analysis: define potential ligands
  expressed_ligands = intersect(ligands,expressed_genes_sender)
  expressed_receptors = intersect(receptors,expressed_genes_receiver)
  if (length(expressed_ligands) == 0){
    stop("No ligands expressed in sender cell")
  }
  if (length(expressed_receptors) == 0){
    stop("No receptors expressed in receiver cell")
  }
  potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
  if (length(potential_ligands) == 0){
    stop("No potentially active ligands")
  }

  if (verbose == TRUE){print("Perform NicheNet ligand activity analysis")}

  # step4 perform NicheNet's ligand activity analysis
  ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  ligand_activities = ligand_activities %>%
    arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected)))

  if(filter_top_ligands == TRUE){
    best_upstream_ligands = ligand_activities %>% top_n(top_n_ligands, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  } else {
    best_upstream_ligands = ligand_activities %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()
  }
  if (verbose == TRUE){print("Infer active target genes of the prioritized ligands")}

  # step5 infer target genes of the top-ranked ligands
  active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = top_n_targets) %>% bind_rows() %>% drop_na()

  if(nrow(active_ligand_target_links_df) > 0){
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff_visualization)
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

    order_targets = order_targets %>% intersect(rownames(active_ligand_target_links))
    order_ligands = order_ligands %>% intersect(colnames(active_ligand_target_links))

    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
    p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) #+ scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.006,0.012))
  } else {
    vis_ligand_target = NULL
    p_ligand_target_network = NULL
    print("no highly likely active targets found for top ligands")
  }
  # combined heatmap: overlay ligand activities
  ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

  rownames(ligand_aupr_matrix) = rownames(ligand_aupr_matrix) %>% make.names()
  colnames(ligand_aupr_matrix) = colnames(ligand_aupr_matrix) %>% make.names()

  vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
  p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)") + theme(legend.text = element_text(size = 9))
  p_ligand_aupr

  figures_without_legend = cowplot::plot_grid(
    p_ligand_aupr + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv",
    nrow = 1,
    rel_widths = c(ncol(vis_ligand_aupr)+10, ncol(vis_ligand_target)))
  legends = cowplot::plot_grid(
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
    ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
    nrow = 1,
    align = "h")

  combined_plot = cowplot::plot_grid(figures_without_legend,
                                     legends,
                                     rel_heights = c(10,2), nrow = 2, align = "hv")

  # ligand-receptor plot
  # get the ligand-receptor network of the top-ranked ligands
  if (verbose == TRUE){print("Infer receptors of the prioritized ligands")}

  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

  lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  if (nrow(lr_network_top_matrix) > 1){
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
  } else {
    order_receptors = rownames(lr_network_top_matrix)
  }
  if (ncol(lr_network_top_matrix) > 1) {
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  } else {
    order_ligands_receptor = colnames(lr_network_top_matrix)
  }

  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  dim(vis_ligand_receptor_network) = c(length(order_receptors), length(order_ligands_receptor))

  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

  # ligand expression Seurat dotplot
  if (length(sender) > 1){
    are_there_senders = TRUE
  }
  if(length(sender) == 1){
    if(sender != "undefined"){
      are_there_senders = TRUE
    } else {
      are_there_senders = FALSE
    }
  }

  if (are_there_senders == TRUE){
    real_makenames_conversion = lr_network$from %>% unique() %>% magrittr::set_names(lr_network$from %>% unique() %>% make.names())
    order_ligands_adapted = real_makenames_conversion[order_ligands]
    names(order_ligands_adapted) = NULL
    rotated_dotplot = DotPlot(seurat_obj %>% subset(idents = sender_celltypes), features = order_ligands_adapted, cols = "RdYlBu") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots

  } else {
    rotated_dotplot = NULL
  }
  return(list(
    ligand_activities = ligand_activities,
    top_ligands = best_upstream_ligands,
    top_targets = active_ligand_target_links_df$target %>% unique(),
    top_receptors = lr_network_top_df_large$to %>% unique(),
    ligand_target_matrix = vis_ligand_target,
    ligand_target_heatmap = p_ligand_target_network,
    ligand_target_df = active_ligand_target_links_df,
    ligand_expression_dotplot = rotated_dotplot,
    ligand_activity_target_heatmap = combined_plot,
    ligand_receptor_matrix = vis_ligand_receptor_network,
    ligand_receptor_heatmap = p_ligand_receptor_network,
    ligand_receptor_df = lr_network_top_df_large %>% rename(ligand = from, receptor = to),
    geneset_oi = geneset_oi,
    background_expressed_genes = background_expressed_genes
  ))
}
#' @title Get log fold change values of genes in cell type of interest
#'
#' @description \code{get_lfc_celltype} Get log fold change of genes between two conditions in cell type of interest when using a Seurat single-cell object.
#'
#' @usage
#' get_lfc_celltype(celltype_oi, seurat_obj, condition_colname, condition_oi, condition_reference, celltype_col = "celltype")
#' #'
#' @param seurat_obj Single-cell expression dataset as Seurat object https://satijalab.org/seurat/.
#' @param celltype_oi Name of celltype of interest. Should be present in the celltype metadata dataframe.
#' @param condition_colname Name of the column in the meta data dataframe that indicates which condition/sample cells were coming from.
#' @param condition_oi Condition of interest. Should be a name present in the "condition_colname" column of the metadata.
#' @param condition_reference The second condition (e.g. reference or steady-state condition). Should be a name present in the "condition_colname" column of the metadata.
#' @param celltype_col Metadata colum name where the cell type identifier is stored. Default: "celltype". If this is NULL, the Idents() of the seurat object will be considered as your cell type identifier.
#' @param ... Additional arguments passed to \code{\link{FindMarkers}}.
#'
#' @return A tbl with the log fold change values of genes. Positive lfc values: higher in condition_oi compared to condition_reference.
#'
#' @import Seurat
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' requireNamespace("dplyr")
#' seuratObj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
#' get_lfc_celltype(seurat_obj = seuratObj, celltype_oi = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", celltype_col = "celltype", expression_pct = 0.10)
#' }
#' @export
#'
get_lfc_celltype = function(celltype_oi, seurat_obj, condition_colname, condition_oi, condition_reference, celltype_col = "celltype", ...){
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  if(!is.null(celltype_col)){
    seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[[celltype_col, drop=TRUE]])
    seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)

  } else {
    seuratObj_sender = subset(seurat_obj, idents = celltype_oi)

  }

  seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname, drop=TRUE]])
  DE_table_sender = FindMarkers(object = seuratObj_sender, ident.1 = condition_oi, ident.2 = condition_reference, ...) %>% rownames_to_column("gene")

  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_sender)

  if(SeuratV4 == TRUE){
    DE_table_sender = DE_table_sender %>% as_tibble() %>% select(-p_val) %>% select(gene, avg_log2FC)
  } else {
    DE_table_sender = DE_table_sender %>% as_tibble() %>% select(-p_val) %>% select(gene, avg_logFC)
  }

  colnames(DE_table_sender) = c("gene",celltype_oi)
  return(DE_table_sender)
}
