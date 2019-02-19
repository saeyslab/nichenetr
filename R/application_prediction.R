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
  targets = intersect(ligand_target_matrix[,ligand] %>% sort(decreasing = T) %>% head(n) %>% names(),geneset)
  print(ligand)
  print(targets)
  if(length(targets) == 0){
    stop("none of the specified possible targets belongs to the top n")
  }
  ligand_target_weighted_df = tibble(ligand = ligand, target = names(ligand_target_matrix[targets,ligand])) %>% inner_join(tibble(target = names(ligand_target_matrix[targets,ligand]), weight = ligand_target_matrix[targets,ligand]), by = "target")
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
