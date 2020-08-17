context("Parameter optimization functions")
#
test_that("Objective function to construct and evaluate the model for parameter optimization is ok", {
  nr_datasources = source_weights_df$source %>% unique() %>% length()

  print("PPR-correct topology")
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
  test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "PPR", TRUE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")

  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)

  print("PPR-ltf cutoff")
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5, "ltf_cutoff" = 0.95, "damping_factor" = 0.5)
  test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "PPR", FALSE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")

  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)

  print("direct ligand-tf links")
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5)
  test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "direct", FALSE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")

  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)

  print("SPL")
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5,"ltf_cutoff" = 0.95)
  test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "SPL", FALSE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")

  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)

  print("direct target links")
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5)
  test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "PPR", FALSE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no", damping_factor = 0)

  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)


  print("hyperparameter optimization")
  test_input = list("lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
  source_weights = source_weights_df$weight
  names(source_weights) = source_weights_df$source
  test_evaluation_optimization = model_evaluation_hyperparameter_optimization(test_input, source_weights, "PPR", TRUE, lr_network, sig_network, gr_network, lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")
  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),4)

  nr_datasources = source_weights_df$source %>% unique() %>% length()
  test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
  test_evaluation_optimization = model_evaluation_optimization_application(test_input, source_weights_df$source %>% unique(), algorithm = "PPR", TRUE, lr_network, sig_network, gr_network, list(convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries)), secondary_targets = FALSE, remove_direct_links = "no", classification_algorithm = "lda", var_imps = FALSE, cv_number = 5, cv_repeats = 4)
  expect_type(test_evaluation_optimization, "double")
  expect_equal(length(test_evaluation_optimization),2)
})

test_that("mlrMBO optimization of a multi-objective function can be performed is ok", {
  library(mlrMBO)
  library(parallelMap)
  model_evaluation_optimization_decoy = function(x, source_names, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",...){

    names(x$source_weights) = source_names
    parameters_setting = list(model_name = "query_design", source_weights = x$source_weights)

    if (correct_topology == TRUE){
      parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = 0, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = TRUE)
    } else {
      parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = x$ltf_cutoff, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = FALSE)
    }

    mean_auroc_target_prediction = parameters_setting$source_weights %>% sum()
    mean_aupr_target_prediction = parameters_setting$source_weights %>% max()
    mean_auroc_ligand_prediction =  parameters_setting$source_weights %>% mean()
    mean_aupr_ligand_prediction = parameters_setting$source_weights %>% min()

    return(c(mean_auroc_target_prediction, mean_aupr_target_prediction, mean_auroc_ligand_prediction, mean_aupr_ligand_prediction))
  }

  additional_arguments_topology_correction = list(source_names = source_weights_df$source %>% unique(), algorithm = "PPR", correct_topology = TRUE,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings = lapply(expression_settings_validation,convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no", cutoff_method = "quantile")

  nr_datasources = additional_arguments_topology_correction$source_names %>% length()

  obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",description = "data source weight and hyperparameter optimization: expensive black-box function", fn = model_evaluation_optimization_decoy, par.set = makeParamSet( makeNumericVectorParam("source_weights", len = nr_datasources, lower = 0, upper = 1), makeNumericVectorParam("lr_sig_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("damping_factor", len = 1, lower = 0, upper = 0.99)), has.simple.signature = FALSE,n.objectives = 4, noisy = FALSE,minimize = c(FALSE,FALSE,FALSE,FALSE))


  if(Sys.info()['sysname'] == "Windows"){
    print("windows - skip test on using parallelized mlrmbo optimization")
  } else {
    mlrmbo_optimization_result = lapply(1,mlrmbo_optimization, obj_fun = obj_fun_multi_topology_correction, niter = 2, ncores = 1, nstart = 100, additional_arguments = additional_arguments_topology_correction)

    optimized_parameters = process_mlrmbo_nichenet_optimization(mlrmbo_optimization_result[[1]],additional_arguments_topology_correction$source_names)
    optimized_parameters = process_mlrmbo_nichenet_optimization(mlrmbo_optimization_result,additional_arguments_topology_correction$source_names)
    expect_type(mlrmbo_optimization_result, "list")
    expect_type(mlrmbo_optimization_result[[1]]$pareto.set[[1]]$source_weights, "double")
    expect_equal(length(mlrmbo_optimization_result[[1]]$pareto.set[[1]]$source_weights),nr_datasources)
    expect_type(optimized_parameters, "list")
    expect_equal(is.data.frame(optimized_parameters$source_weight_df), TRUE)
  }
})

test_that("Data source weights can be estimated from leave-one-in and leave-one-out performances", {

  # run characterization loi
  settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
  weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df)
  weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR", damping_factor = 0.2, correct_topology = TRUE)

  job_characterization_loi = lapply(weights_settings_loi[1:4], evaluate_model,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings,calculate_popularity_bias_target_prediction = FALSE, calculate_popularity_bias_ligand_prediction = FALSE, ncitations)
  loi_performances = process_characterization_target_prediction_average(job_characterization_loi)

  # run characterization loo
  weights_settings_loo = prepare_settings_leave_one_out_characterization(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df)
  weights_settings_loo = lapply(weights_settings_loo,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR", damping_factor = 0.2, correct_topology = TRUE)

  job_characterization_loo = lapply(weights_settings_loo[1:4], evaluate_model,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings,calculate_popularity_bias_target_prediction = FALSE, calculate_popularity_bias_ligand_prediction = FALSE,ncitations)
  loo_performances = process_characterization_target_prediction_average(job_characterization_loo)

  # run the regression
  sources_oi = c("kegg_cytokines")
  output = estimate_source_weights_characterization(loi_performances,loo_performances,source_weights_df %>% filter(source != "kegg_cytokines"), sources_oi, random_forest =FALSE)
})
