context("Model evaluation functions: ligand prediction")
test_that("Convert expression settings to settings", {
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(expression_settings_validation, combination = FALSE)
  ligands = unlist(ligands)

  settings_ligand_pred = convert_settings_ligand_prediction(settings = settings, all_ligands = ligands, validation = TRUE, single = TRUE)
  expect_type(settings_ligand_pred,"list")

  expect_type(settings_ligand_pred[[1]]$ligand,"character")
  expect_equal(length(settings_ligand_pred[[1]]$ligand),1)


  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")


  expect_equal(length(settings_ligand_pred),length(ligands)*length(settings))


  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = TRUE, single = FALSE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(settings))
  expect_type(settings_ligand_pred[[1]]$ligand,"character")
  expect_equal(length(settings_ligand_pred[[1]]$ligand),1)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = FALSE, single = TRUE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(ligands)*length(settings))
  expect_equal(settings_ligand_pred[[1]]$ligand,NULL)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = FALSE, single = FALSE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(settings))
  expect_equal(settings_ligand_pred[[1]]$ligand,NULL)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

})
test_that("Get ligand importances: single", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = get_single_ligand_importances(settings_ligand_pred,ligand_target_matrix)
  expect_type(ligand_importances,"list")
  ligand_importances2 = get_single_ligand_importances(settings_ligand_pred,ligand_target_matrix %>% make_discrete_ligand_target_matrix())
  expect_type(ligand_importances2,"list")
})
test_that("Get ligand importances: multi", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances_glm = get_multi_ligand_importances(settings_ligand_pred,ligand_target_matrix, algorithm = "glm")
  expect_type(ligand_importances_glm,"list")
  ligand_importances_glm2 = get_multi_ligand_importances(settings_ligand_pred,ligand_target_matrix, algorithm = "glm",filter_genes = TRUE)
  expect_type(ligand_importances_glm2,"list")
})











