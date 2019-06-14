context("Model evaluation functions: target gene predictions")

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
  expect_type(lapply(expression_settings_validation,convert_expression_settings_evaluation_regression) %>% .[[1]] %>% .$response,"double")

})

test_that("Convert gene list to settings", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  all_genes = unique(c(weighted_networks$gr$from,weighted_networks$gr$to,weighted_networks$lr_sig$from, weighted_networks$lr_sig$to))
  gene_list = c("ID1","ID2","ID3")
  setting = list(convert_gene_list_settings_evaluation(gene_list = c("ID1","ID2","ID3"), name = "test",ligands_oi = "TGFB1", background = all_genes))
  expect_type(setting,"list")
  expect_type(setting %>% .[[1]] %>% .$response,"logical")
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
  performances_target_prediction_medianed = ligands %>% unlist() %>% unique() %>% lapply(wrapper_average_performances, performances,"median") %>% bind_rows() %>% drop_na()
  performances_target_prediction_averaged = ligands %>% unlist() %>% unique() %>% lapply(wrapper_average_performances, performances,"mean") %>% bind_rows() %>% drop_na()
  expect_type(performances_target_prediction_medianed,"list")
  expect_type(performances_target_prediction_averaged,"list")
})

test_that("Evaluate target gene value prediction: regression", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = extract_ligands_from_settings(expression_settings_validation)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, algorithm = "PPR", damping_factor = 0.5)
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation_regression)
  performances = bind_rows(lapply(settings,evaluate_target_prediction_regression,ligand_target_matrix))
  expect_type(performances,"list")


})
test_that("Evaluate target gene prediction: interpretation", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = extract_ligands_from_settings(expression_settings_validation)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, algorithm = "PPR", damping_factor = 0.5)
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  performances = lapply(settings,evaluate_target_prediction_interprete,ligand_target_matrix) %>% .[[1]]
  expect_type(performances,"list")
  performances_discrete = lapply(settings,evaluate_target_prediction_interprete,ligand_target_matrix %>% make_discrete_ligand_target_matrix) %>% .[[1]]
  expect_type(performances_discrete,"list")

  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation_regression)
  performances2 = lapply(settings,evaluate_target_prediction_interprete,ligand_target_matrix) %>% .[[1]]
  expect_type(performances2,"list")

})
test_that("Evaluate target gene prediction multiple ligands", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  setting = convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
  ligands = extract_ligands_from_settings(setting)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  performances = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix,ligands_position = "cols",algorithm = "glm") %>% .[[1]]
  expect_type(performances,"list")
  performances_discrete = lapply(setting,evaluate_multi_ligand_target_prediction,make_discrete_ligand_target_matrix(ligand_target_matrix),ligands_position = "cols",algorithm = "glm" ) %>% .[[1]]
  expect_type(performances_discrete,"list")
  performances_discrete2 = lapply(setting,evaluate_multi_ligand_target_prediction,make_discrete_ligand_target_matrix(ligand_target_matrix),ligands_position = "cols",algorithm = "glm",continuous = FALSE) %>% .[[1]]
  expect_type(performances_discrete2,"list")

  performances = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix %>% t(),algorithm = "glm",ligands_position = "rows") %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix,ligands_position = "cols",algorithm = "glm",ignore_errors = TRUE) %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix,ligands_position = "cols",algorithm = "glm",cv = FALSE) %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix,ligands_position = "cols",algorithm = "glm",var_imps = TRUE) %>% .[[1]]
  expect_type(performances,"list")
})
test_that("Evaluate target gene prediction multiple ligands: regression", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
  setting = convert_expression_settings_evaluation_regression(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
  ligands = extract_ligands_from_settings(setting)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
  performances = lapply(setting,evaluate_multi_ligand_target_prediction_regression,ligand_target_matrix,ligands_position = "cols",algorithm = "lm") %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction_regression,ligand_target_matrix %>% t(),algorithm = "lm",ligands_position = "rows") %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction_regression,ligand_target_matrix,ligands_position = "cols",algorithm = "lm",ignore_errors = TRUE) %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction_regression,ligand_target_matrix,ligands_position = "cols",algorithm = "lm",cv = FALSE) %>% .[[1]]
  expect_type(performances,"list")
  performances = lapply(setting,evaluate_multi_ligand_target_prediction_regression,ligand_target_matrix,ligands_position = "cols",algorithm = "lm",var_imps = TRUE) %>% .[[1]]
  expect_type(performances,"list")
})

