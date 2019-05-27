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
test_that("Get ligand importances: single + evaluation", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)

  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances,ligand_target_matrix))
  expect_type(ligand_importances,"list")
  expect_gte(nrow(ligand_importances),1)
  ligand_importances2 = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix %>% make_discrete_ligand_target_matrix()))
  expect_type(ligand_importances2,"list")
  expect_gte(nrow(ligand_importances2),1)

  evaluation_single = evaluate_single_importances_ligand_prediction(ligand_importances,"median")
  expect_type(evaluation_single,"list")
  expect_gte(nrow(evaluation_single),1)
  evaluation_single2= evaluate_single_importances_ligand_prediction(ligand_importances,"mean")
  expect_type(evaluation_single2,"list")
  expect_gte(nrow(evaluation_single2),1)

  evaluation_single3 = ligand_importances$setting %>% unique() %>% lapply(function(x){x}) %>% lapply(wrapper_evaluate_single_importances_ligand_prediction,ligand_importances) %>% bind_rows() %>% inner_join(ligand_importances %>% distinct(setting,ligand))
  expect_type(evaluation_single3,"list")
  expect_gte(nrow(evaluation_single3),1)

  evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")
  expect_type(evaluation,"list")
  expect_gte(nrow(evaluation$performances),1)


})
test_that("Get ligand importances: multi", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances_glm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix, algorithm = "glm"))
  expect_type(ligand_importances_glm,"list")
  expect_gte(nrow(ligand_importances_glm),1)

  ligand_importances_glm2 = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix, algorithm = "glm",filter_genes = TRUE))
  expect_type(ligand_importances_glm2,"list")
  expect_gte(nrow(ligand_importances_glm2),1)

})

test_that("Get ligand importances: random forest", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances_rf = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_rf_importances, ligand_target_matrix, ntrees = 100, mtry = 2))

  expect_type(ligand_importances_rf,"list")
  expect_gte(nrow(ligand_importances_rf),1)

})

test_that("Model-based ligand activity prediction", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix))
  evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")

  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = FALSE, single = TRUE)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix,known = FALSE))
  activity_predictions = model_based_ligand_activity_prediction(ligand_importances, evaluation$model, "median")
  expect_type(activity_predictions,"list")
  expect_gte(nrow(activity_predictions),1)
})

test_that("Get ligand importances regression: single + evaluation", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)

  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances_regression,ligand_target_matrix))
  expect_type(ligand_importances,"list")
  expect_gte(nrow(ligand_importances),1)

  evaluation_single = evaluate_single_importances_ligand_prediction(ligand_importances,"median")
  expect_type(evaluation_single,"list")
  expect_gte(nrow(evaluation_single),1)
  evaluation_single2= evaluate_single_importances_ligand_prediction(ligand_importances,"mean")
  expect_type(evaluation_single2,"list")
  expect_gte(nrow(evaluation_single2),1)

  evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")
  expect_type(evaluation,"list")
  expect_gte(nrow(evaluation$performances),1)

})
test_that("Get ligand importances regression: multi", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances_lm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances_regression, ligand_target_matrix, algorithm = "lm"))
  expect_type(ligand_importances_lm,"list")
  expect_gte(nrow(ligand_importances_lm),1)

  ligand_importances_lm2 = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances_regression, ligand_target_matrix, algorithm = "lm",filter_genes = TRUE))
  expect_type(ligand_importances_lm2,"list")
  expect_gte(nrow(ligand_importances_lm2),1)

})

test_that("Get ligand importances regression: random forest", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances_rf = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_rf_importances_regression, ligand_target_matrix, ntrees = 100, mtry = 2))

  expect_type(ligand_importances_rf,"list")
  expect_gte(nrow(ligand_importances_rf),1)

})

test_that("Expression setting converters: tf upstream analysis", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_tf_pred = convert_settings_tf_prediction(settings, all_tfs = c("SMAD1","STAT1","RELA"), single = TRUE)
  # show how this function can be used to predict activities of TFs
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  tf_target = construct_tf_target_matrix(weighted_networks, tfs_as_cols = TRUE, standalone_output = TRUE)
  tf_importances = dplyr::bind_rows(lapply(settings_tf_pred,get_single_ligand_importances,tf_target,known = FALSE))

  expect_type(tf_importances,"list")

})
test_that("Expression setting converters: top n ligands", {
  settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)

  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
  evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")

  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = FALSE, single = TRUE)
  ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix, known = FALSE))
  settings = lapply(settings,convert_settings_topn_ligand_prediction, importances = ligand_importances, model = evaluation$model, n = 3, normalization = "median" )

  expect_type(settings,"list")
})








