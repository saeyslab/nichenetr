context("Popularity bias evaluation functions")

test_that("Ligand popularity bias functions are ok for target gene prediction", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(settings)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  performances = dplyr::bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix))

##  ncitations = get_ncitations_genes()
  expect_equal(is.data.frame(ncitations),TRUE)

  performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances,ncitations)
  expect_equal(is.data.frame(performances_ligand_popularity),TRUE)

  slopes_df = get_slope_ligand_popularity("auroc",performances_ligand_popularity)
  expect_equal(is.data.frame(slopes_df),TRUE)
})
test_that("Target gene popularity bias functions are ok for target gene prediction", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(settings)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#  ncitations = get_ncitations_genes()

  performances_target_bins_popularity = evaluate_target_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
  expect_equal(is.data.frame(performances_target_bins_popularity),TRUE)

  slopes_df = get_slope_target_gene_popularity("auroc",performances_target_bins_popularity)
  expect_equal(is.data.frame(slopes_df),TRUE)

  slopes_df2 = get_slope_target_gene_popularity("auroc",performances_target_bins_popularity,method = "all")
  expect_equal(is.data.frame(slopes_df2),TRUE)
})
test_that("Ligand popularity bias functions are ok for ligand activity prediction", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))

  ligand_activity_popularity_bias = lapply(0:3,ligand_activity_performance_top_i_removed, ligand_importances, ncitations) %>% bind_rows()
  slopes_auroc = get_ligand_slope_ligand_prediction_popularity("auroc",ligand_activity_popularity_bias)
  slopes_df = ligand_activity_popularity_bias %>% select(-importance_measure, -popularity_index) %>% colnames() %>% lapply(.,get_ligand_slope_ligand_prediction_popularity,ligand_activity_popularity_bias) %>% bind_rows()
  ##  ncitations = get_ncitations_genes()
  expect_type(ligand_activity_popularity_bias$popularity_index, "double")
  expect_equal(is.data.frame(slopes_df),TRUE)
  expect_equal(is.data.frame(slopes_auroc),TRUE)

})
test_that("Target gene popularity bias functions are ok for ligand activity prediction", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(settings)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  #  ncitations = get_ncitations_genes()

  performances_target_bins_popularity = evaluate_ligand_prediction_per_bin(3,settings,ligand_target_matrix,ncitations,cutoff_method = "quantile") # make quantile-based discrete ligand-target matrix!
  expect_equal(is.data.frame(performances_target_bins_popularity),TRUE)

  slopes_df = get_slope_target_gene_popularity_ligand_prediction("auroc",performances_target_bins_popularity)
  expect_equal(is.data.frame(slopes_df),TRUE)

  slopes_df_all = performances_target_bins_popularity %>% select(-importance_measure,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity_ligand_prediction,performances_target_bins_popularity) %>% bind_rows()
  expect_equal(is.data.frame(slopes_df_all),TRUE)

})
