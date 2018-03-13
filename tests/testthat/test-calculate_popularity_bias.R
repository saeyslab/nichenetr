context("Popularity bias evaluation functions")

test_that("Ligand popularity bias functions are ok", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
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
test_that("Target gene popularity bias functions are ok", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
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
