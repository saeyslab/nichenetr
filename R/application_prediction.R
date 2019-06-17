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
  return(ligand_importances %>% select(test_ligand,auroc,aupr,pearson))
}
#' @title Infer weighted active ligand-target links between a possible ligand and target genes of interest
#'
#' @description \code{get_weighted_ligand_target_links} Infer active ligand target links between possible lignands and genes belonging to a gene set of interest: consider the intersect between the top n targets of a ligand and the gene set.
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
  } else if (length(targets) == 1) {
    ligand_target_weighted_df = tibble(ligand = ligand, target = targets, weight = ligand_target_matrix[targets,ligand])
  } else {
    ligand_target_weighted_df = tibble(ligand = ligand, target = names(ligand_target_matrix[targets,ligand])) %>% inner_join(tibble(target = names(ligand_target_matrix[targets,ligand]), weight = ligand_target_matrix[targets,ligand]), by = "target")
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
#'  = c("TNF","BMP2","IL4")
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
  ligand_target_vis_filtered = ligand_target_vis[ligand_target_vis %>% apply(1,sum) > 0,ligand_target_vis %>% apply(2,sum) > 0]

  distoi = dist(1-cor(t(ligand_target_vis_filtered)))
  hclust_obj = hclust(distoi, method = "ward.D2")
  order_targets = hclust_obj$labels[hclust_obj$order]

  distoi_targets = dist(1-cor(ligand_target_vis_filtered))
  hclust_obj = hclust(distoi_targets, method = "ward.D2")
  order_ligands = hclust_obj$labels[hclust_obj$order]

  vis_ligand_target_network = ligand_target_vis_filtered[order_targets,order_ligands]
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
  predicted_positive = affected_gene_predictions %>% arrange(-prediction) %>% filter(prediction >= quantile(prediction,0.95)) %>% group_by(response) %>% count() %>% rename(positive_prediction = n) %>% rename(true_target = response)
  all = affected_gene_predictions %>% arrange(-prediction) %>% rename(true_target = response) %>% group_by(true_target) %>% count()
  inner_join(all,predicted_positive, by = "true_target") %>% mutate(fraction_positive_predicted = positive_prediction/n)
}
#' @title Perform a Fisher's exact test to determine whether genes belonging to the gene set of interest are more likely to be part of the top-predicted targets.
#'
#' @description \code{calculate_fraction_top_predicted_fisher} Performs a Fisher's exact test to determine whether genes belonging to the gene set of interest are more likely to be part of the top-predicted targets.
#'
#' @usage
#' calculate_fraction_top_predicted_fisher(affected_gene_predictions, quantile_cutoff = 0.95,output = "p-value")
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
  results_df = inner_join(all,predicted_positive, by = "response") %>% mutate(fraction_positive_predicted = positive_prediction/n)
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
  single_ligand_activities_pearson_norm = ligand_activities %>%
    group_by(setting) %>%
    mutate(pearson = nichenetr::scaling_modified_zscore(pearson)) %>%
    ungroup() %>%
    rename(cell = setting, ligand = test_ligand) %>%
    distinct(cell,ligand,pearson)

  single_ligand_activities_pearson_norm_df = single_ligand_activities_pearson_norm %>%
    spread(cell, pearson,fill = min(.$pearson))

  single_ligand_activities_pearson_norm_matrix = single_ligand_activities_pearson_norm_df  %>%
    select(-ligand) %>%
    t() %>%
    magrittr::set_colnames(single_ligand_activities_pearson_norm_df$ligand)

  single_cell_ligand_activities_pearson_norm_df = single_ligand_activities_pearson_norm_matrix %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    tbl_df()
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
