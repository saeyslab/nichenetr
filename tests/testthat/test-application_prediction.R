context("NicheNet analysis on Seurat objects")
test_that("Seurat wrapper works", {

  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
  seurat_object_lite = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "up")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "down")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "all", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse", filter_top_ligands = FALSE)
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "up")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "down")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "all", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse", filter_top_ligands = FALSE)
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "up")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse",geneset = "down")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = "all", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse", filter_top_ligands = FALSE)
  expect_type(nichenet_output,"list")

  seurat_object_lite@meta.data$aggregate = seurat_object_lite@meta.data$aggregate %>% as.factor()
  seurat_object_lite@meta.data$celltype = seurat_object_lite@meta.data$celltype %>% as.factor()

  nichenet_output = nichenet_seuratobj_aggregate(seurat_obj = seurat_object_lite, receiver = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = c("Mono"), ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_aggregate_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "CD8 T", condition_oi = "LCMV", condition_reference = "SS", condition_colname = "aggregate", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = seurat_object_lite, receiver_affected = "CD8 T", receiver_reference = "Mono", sender = "undefined", ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network, organism = "mouse")
  expect_type(nichenet_output,"list")

  lfc_output = get_lfc_celltype(seurat_obj = seurat_object_lite, celltype_oi = "CD8 T", condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", expression_pct = 0.10)
  expect_type(lfc_output,"list")


})

context("Target gene prediction functions for application")
test_that("Gene-cluster assignments can be converted to settings", {

  genes_clusters = c("TGFB1" = 1,"TGFB2" = 1,"TGFB3" = 2)
  cluster_settings = lapply(seq(length(unique(genes_clusters))), convert_cluster_to_settings, cluster_vector = genes_clusters, setting_name = "example", setting_from = "BMP2")

  expect_type(cluster_settings,"list")
  expect_type(cluster_settings[[1]]$response,"logical")
  expect_type(cluster_settings[[1]]$from,"character")
  expect_type(cluster_settings[[1]]$name,"character")
  expect_type(cluster_settings[[2]]$response,"logical")

})

context("Ligand activity prediction for application")
test_that("Ligand activities can be predicted well", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","IL4")
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  potential_ligands = c("TNF","IL4")
  geneset = c("VCAM1","SOCS3", "IRF1","LTBR","STAT4")
  background_expressed_genes = c("VCAM1","SOCS3","IRF1","ICAM1","ID1","ID2","ID3","LTBR",'STAT4')
  ligand_activities = predict_ligand_activities(geneset = geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

  expect_type(ligand_activities,"list")
  expect_type(ligand_activities$test_ligand,"character")
  expect_type(ligand_activities$pearson,"double")

  active_ligand_target_links_df = potential_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset, ligand_target_matrix = ligand_target_matrix, n = 500) %>% bind_rows()
  expect_type(active_ligand_target_links_df,"list")

  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.001)
  expect_type(active_ligand_target_links,"double")

})

context("Assessment target gene prediction by multi-ligand models")
test_that("Target gene prediction can be predicted by multi-ligand models", {

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","IL4")
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  potential_ligands = c("TNF","IL4")
  geneset = c("VCAM1","SOCS3", "IRF1","LTBR","STAT4")
  background_expressed_genes = c("VCAM1","SOCS3","IRF1","ICAM1","ID1","ID2","ID3","LTBR",'STAT4')
  gene_predictions_list = seq(10) %>% lapply(assess_rf_class_probabilities,2, geneset = geneset,background_expressed_genes = background_expressed_genes,ligands_oi = potential_ligands,ligand_target_matrix = ligand_target_matrix)

  expect_type(gene_predictions_list,"list")
  expect_type(gene_predictions_list[[1]]$gene,"character")
  expect_type(gene_predictions_list[[1]]$prediction,"double")

  target_prediction_performances_cv = gene_predictions_list %>% lapply(classification_evaluation_continuous_pred_wrapper) %>% bind_rows() %>% mutate(round=seq(1:nrow(.)))
  expect_type(target_prediction_performances_cv,"list")
  expect_type(target_prediction_performances_cv$auroc,"double")

  top_genes_df = seq(length(gene_predictions_list)) %>% lapply(get_top_predicted_genes,gene_predictions_list) %>% reduce(full_join, by = c("gene","true_target"))
  expect_type(top_genes_df,"list")

  target_prediction_performances_discrete_cv = gene_predictions_list %>% lapply(calculate_fraction_top_predicted,quantile_cutoff = 0.66) %>% bind_rows()
  expect_type(target_prediction_performances_discrete_cv,"list")

  target_prediction_performances_fisher_pval = gene_predictions_list %>% lapply(calculate_fraction_top_predicted_fisher,quantile_cutoff = 0.66) %>% unlist() %>% mean()
  expect_type(target_prediction_performances_fisher_pval,"double")

})

context("Single-cell ligand activity prediction functions")
test_that("Single-cell ligand activity prediction functions work a bit OK", {
  x = matrix(rnorm(200*2, sd = 10, mean = 5), ncol = 2)
  x_scaled = scale_quantile(x)
  expect_type(x_scaled,"double")
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","BMP2","IL4")
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
  potential_ligands = c("TNF","BMP2","IL4")
  genes = c("SOCS2","SOCS3","IRF1","ICAM1","ID1","ID2","ID3")
  cell_ids = c("cell1","cell2","cell3","cell4")
  set.seed(1)
  expression_scaled = matrix(rnorm(length(genes)*length(cell_ids), sd = 0.5, mean = 0.5), nrow = length(cell_ids))
  rownames(expression_scaled) = cell_ids
  colnames(expression_scaled) = genes

  settings = convert_single_cell_expression_to_settings(cell_id = cell_ids[1], expression_matrix = expression_scaled, setting_name = "test", setting_from = potential_ligands)
  expect_type(settings,"list")

  ligand_activities = predict_single_cell_ligand_activities(cell_ids = cell_ids, expression_scaled = expression_scaled, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
  expect_type(ligand_activities,"list")
  expect_type(ligand_activities$test_ligand,"character")
  expect_type(ligand_activities$setting,"character")
  expect_type(ligand_activities$pearson,"double")

  normalized_ligand_activities = normalize_single_cell_ligand_activities(ligand_activities)
  expect_type(normalized_ligand_activities,"list")

  cell_scores_tbl = tibble(cell = cell_ids, score = c(1,4,2,3))
  regression_analysis_output = single_ligand_activity_score_regression(normalized_ligand_activities,cell_scores_tbl)
  expect_type(regression_analysis_output,"list")

  cell_scores_tbl = tibble(cell = cell_ids, score = c(TRUE, FALSE, TRUE, FALSE))
  classification_analysis_output = single_ligand_activity_score_classification(normalized_ligand_activities,cell_scores_tbl)
  expect_type(classification_analysis_output,"list")

})


