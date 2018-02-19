context("Model evaluation functions")

test_that("Extract ligands from settings", {
  settings = list(
    list(name = "test1", from = "BMP2"),
    list(name = "test2", from = c("BMP2","TNF"))
  )
  expect_equal(extract_ligands_from_settings(settings),list("BMP2",c("BMP2","TNF"), "TNF"))
})
test_that("Convert expression settings to settings", {
  expect_type(lapply(expression_settings_validation,convert_expression_settings_evaluation),"list")
  expect_type(lapply(expression_settings_validation,convert_expression_settings_evaluation) %>% .[[1]] %>% .$response,"logical")
})
test_that("Evaluate target gene prediction", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = extract_ligands_from_settings(expression_settings_validation)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, algorithm = "PPR", damping_factor = 0.5)
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  performances = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix))
  expect_type(performances,"list")
  performances_discrete = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix %>% make_discrete_ligand_target_matrix))
  expect_type(performances_discrete,"list")

})



