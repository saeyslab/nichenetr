context("Data source and model characterization functions")
#
test_that("Settings leave-one-in can be constructed and hyperparameters added", {
  weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)

  expect_type(weights_settings_loi, "list")
  expect_type(weights_settings_loi[[2]], "list")
  expect_type(weights_settings_loi[[2]]$source_weights, "double")
  expect_type(weights_settings_loi[[2]]$model_name, "character")


  weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(weights_settings_loi[[2]]$lr_sig_hub, "double")
  expect_type(weights_settings_loi[[2]]$gr_hub, "double")
  expect_type(weights_settings_loi[[2]]$damping_factor, "double")
  expect_type(weights_settings_loi[[2]]$ltf_cutoff, "double")

  expect_type(weights_settings_loi[[2]]$algorithm, "character")
  expect_type(weights_settings_loi[[2]]$correct_topology, "logical")

})

test_that("Leave-one-in models can be evaluated and results of this further processed", {

  settings = lapply(expression_settings_validation[1:20], convert_expression_settings_evaluation)
  weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
  weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  output_characterization = lapply(weights_settings_loi[1:2],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, cutoff_method = "quantile") # make discrete ligand-target matrix via the quantile-method!

  target_prediction_performances = process_characterization_target_prediction(output_characterization)
  target_prediction_performances_average = process_characterization_target_prediction_average(output_characterization)
  ligand_prediction_performances_single = process_characterization_ligand_prediction_single_measures(output_characterization)
  popularity_slopes_target_prediction_performances = process_characterization_popularity_slopes_target_prediction(output_characterization)
  popularity_slopes_ligand_prediction_performances = process_characterization_popularity_slopes_ligand_prediction(output_characterization)

  expect_type(target_prediction_performances, "list")
  expect_type(target_prediction_performances_average, "list")
  expect_type(ligand_prediction_performances_single, "list")
  expect_type(popularity_slopes_target_prediction_performances, "list")
  expect_type(popularity_slopes_ligand_prediction_performances, "list")

  expect_gte(nrow(target_prediction_performances),1)
  expect_gte(nrow(target_prediction_performances_average),1)
  expect_gte(nrow(ligand_prediction_performances_single),1)
  expect_gte(nrow(popularity_slopes_target_prediction_performances),1)
  expect_gte(nrow(popularity_slopes_ligand_prediction_performances),1)

})
