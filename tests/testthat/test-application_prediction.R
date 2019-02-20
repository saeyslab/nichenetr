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

