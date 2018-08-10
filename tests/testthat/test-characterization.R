context("Data source and model characterization functions")
#
test_that("Settings leave-one-in can be constructed and hyperparameters added", {
  print("Settings leave-one-in can be constructed and hyperparameters added")
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
test_that("Settings leave-one-out can be constructed and hyperparameters added", {
  print("Settings leave-one-out can be constructed and hyperparameters added")
  weights_settings_loo = prepare_settings_leave_one_out_characterization(lr_network,sig_network, gr_network, source_weights_df)

  expect_type(weights_settings_loo, "list")
  expect_type(weights_settings_loo[[2]], "list")
  expect_type(weights_settings_loo[[2]]$source_weights, "double")
  expect_type(weights_settings_loo[[2]]$model_name, "character")

  weights_settings_loo = lapply(weights_settings_loo,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(weights_settings_loo[[2]]$lr_sig_hub, "double")
  expect_type(weights_settings_loo[[2]]$gr_hub, "double")
  expect_type(weights_settings_loo[[2]]$damping_factor, "double")
  expect_type(weights_settings_loo[[2]]$ltf_cutoff, "double")

  expect_type(weights_settings_loo[[2]]$algorithm, "character")
  expect_type(weights_settings_loo[[2]]$correct_topology, "logical")

})
test_that("Settings one-vs-one can be constructed and hyperparameters added", {
  print("Settings one-vs-one can be constructed and hyperparameters added")
  weights_settings_ovo = prepare_settings_one_vs_one_characterization(lr_network,sig_network, gr_network)

  expect_type(weights_settings_ovo, "list")
  expect_type(weights_settings_ovo[[2]], "list")
  expect_type(weights_settings_ovo[[2]]$source_weights, "double")
  expect_type(weights_settings_ovo[[2]]$model_name, "character")
  expect_type(weights_settings_ovo[[2]]$lr_sig_source, "character")
  expect_type(weights_settings_ovo[[2]]$gr_source, "character")

  weights_settings_ovo = lapply(weights_settings_ovo,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(weights_settings_ovo[[2]]$lr_sig_hub, "double")
  expect_type(weights_settings_ovo[[2]]$gr_hub, "double")
  expect_type(weights_settings_ovo[[2]]$damping_factor, "double")
  expect_type(weights_settings_ovo[[2]]$ltf_cutoff, "double")

  expect_type(weights_settings_ovo[[2]]$algorithm, "character")
  expect_type(weights_settings_ovo[[2]]$correct_topology, "logical")

  weights_settings_ovo = prepare_settings_one_vs_one_characterization(lr_network,sig_network, gr_network, lr_network_separate = TRUE)
  expect_type(weights_settings_ovo, "list")
  expect_type(weights_settings_ovo[[2]], "list")
  expect_type(weights_settings_ovo[[2]]$source_weights, "double")
  expect_type(weights_settings_ovo[[2]]$model_name, "character")
  expect_type(weights_settings_ovo[[2]]$lr_source, "character")
  expect_type(weights_settings_ovo[[2]]$sig_source, "character")
  expect_type(weights_settings_ovo[[2]]$gr_source, "character")

  weights_settings_ovo = lapply(weights_settings_ovo,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(weights_settings_ovo[[2]]$lr_sig_hub, "double")
  expect_type(weights_settings_ovo[[2]]$gr_hub, "double")
  expect_type(weights_settings_ovo[[2]]$damping_factor, "double")
  expect_type(weights_settings_ovo[[2]]$ltf_cutoff, "double")

  expect_type(weights_settings_ovo[[2]]$algorithm, "character")
  expect_type(weights_settings_ovo[[2]]$correct_topology, "logical")

})
test_that("Leave-one-in models can be evaluated and results of this further processed", {
  print("Leave-one-in models can be evaluated and results of this further processed")
  settings = lapply(expression_settings_validation[1:5], convert_expression_settings_evaluation)
  weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
  weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  output_characterization = lapply(weights_settings_loi[2:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, cutoff_method = "quantile")
  expect_type(output_characterization, "list")

  output_characterization_partial = lapply(weights_settings_loi[2:3],evaluate_model_cv,lr_network,sig_network, gr_network,settings, cutoff_method = "quantile")
  expect_type(output_characterization_partial, "list")

  output_characterization_random = lapply(weights_settings_loi[2:3],evaluate_random_model, lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = FALSE, calculate_popularity_bias_ligand_prediction = FALSE, ncitations, cutoff_method = "quantile")
  expect_type(output_characterization_random, "list")


  output_characterization_application = lapply(weights_settings_loi[1:2],evaluate_model_application,lr_network,sig_network, gr_network,settings[1:3], cutoff_method = "quantile") # make discrete ligand-target matrix via the quantile-method!
  expect_type(output_characterization_application, "list")

  target_prediction_performances = process_characterization_target_prediction(output_characterization)
  target_prediction_performances_output_characterization_application = process_characterization_target_prediction(output_characterization_application)
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
test_that("Leave-one-in models can be evaluated for multi-ligand applications", {
  print("Leave-one-in models can be evaluated for multi-ligand applications")
  settings = convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
  weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
  weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  output_characterization_application = lapply(weights_settings_loi[1:3],evaluate_model_application_multi_ligand,lr_network,sig_network, gr_network,settings, classification_algorithm = "lda", var_imps = FALSE, cv_number = 5, cv_repeats = 4)
  expect_type(output_characterization_application, "list")
  target_prediction_performances = process_characterization_target_prediction(output_characterization_application)
  expect_gte(nrow(target_prediction_performances),1)

})

test_that("Influence individual data source on ligand-target scores can be assessed", {
  print("Influence individual data source on ligand-target scores can be assessed")
  ligands =  extract_ligands_from_settings(expression_settings_validation[1:4])
  output = assess_influence_source("ontogenet_coarse", lr_network,sig_network, gr_network, source_weights_df, ligands,lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(output, "list")
  expect_type(output[[2]]$targets_higher, "double")
  first = log(output[[2]]$targets_higher[1])
  second = log(output[[2]]$targets_higher[2])
  expect_gt(first - second,0 )

  output = assess_influence_source("ontogenet_coarse", lr_network,sig_network, gr_network, source_weights_df, ligands, rankings = TRUE, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(output, "list")
  expect_type(output[[2]]$targets_higher, "integer")
  first = output[[2]]$targets_higher[1]
  second = output[[2]]$targets_higher[2]
  expect_lt(first - second,0 )

  output = assess_influence_source("ontogenet_coarse", lr_network,sig_network, gr_network, source_weights_df, ligands, rankings = TRUE, matrix_output = TRUE, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(output, "list")
  expect_equal(output[[1]]$model %>% is.matrix(), TRUE)

  output = assess_influence_source("ontogenet_coarse", lr_network,sig_network, gr_network, source_weights_df, ligands, rankings = FALSE,matrix_output = TRUE, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
  expect_type(output, "list")
  expect_equal(output[[1]]$model %>% is.matrix(), TRUE)

})



