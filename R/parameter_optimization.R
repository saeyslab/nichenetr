#' @title Construct and evaluate a ligand-target model given input parameters with the purpose of parameter optimization.
#'
#' @description \code{model_evaluation_optimization} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input validation settings (average performance for both target gene prediction and ligand activity prediction, as measured via the auroc and aupr).
#'
#' @usage
#' model_evaluation_optimization(x, source_names, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",...)
#'
#' @inheritParams evaluate_model
#' @inheritParams construct_ligand_target_matrix
#' @param x A list containing the following elements. $source_weights: numeric vector representing the weight for each data source; $lr_sig_hub: hub correction factor for the ligand-signaling network; $gr_hub: hub correction factor for the gene regulatory network; $damping_factor: damping factor in the PPR algorithm if using PPR and optionally $ltf_cutoff: the cutoff on the ligand-tf matrix. For more information about these parameters: see \code{construct_ligand_target_matrix} and \code{apply_hub_correction}.
#' @param source_names Character vector containing the names of the data sources. The order of data source names accords to the order of weights in x$source_weights.
#' @param correct_topology This parameter indicates whether the PPR-constructed ligand-target matrix will be subtracted by a PR-constructed target matrix. TRUE or FALSE.
#' @param ... Additional arguments to \code{make_discrete_ligand_target_matrix}.
#'
#' @return A numeric vector of length 4 containing the average auroc for target gene prediction, average aupr (corrected for TP fraction) for target gene prediction, average auroc for ligand activity prediction and average aupr for ligand activity prediction.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' nr_datasources = source_weights_df$source %>% unique() %>% length()
#' test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
# test_evaluation_optimization = model_evaluation_optimization(test_input, source_weights_df$source %>% unique(), "PPR", TRUE, lr_network, sig_network, gr_network, lapply(expression_settings_validation,convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")
#' }
#'
#' @export
#'
model_evaluation_optimization = function(x, source_names, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",...){

  requireNamespace("dplyr")

  #input check
  if (!is.list(x))
    stop("x should be a list!")
  if (!is.numeric(x$source_weights))
    stop("x$source_weights should be a numeric vector")
  if (x$lr_sig_hub < 0 | x$lr_sig_hub > 1)
    stop("x$lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (x$gr_hub < 0 | x$gr_hub > 1)
    stop("x$gr_hub must be a number between 0 and 1 (0 and 1 included)")
  if(is.null(x$ltf_cutoff)){
    if(correct_topology == FALSE)
      warning("Did you not forget to give a value to x$ltf_cutoff?")
  } else {
    if (x$ltf_cutoff < 0 | x$ltf_cutoff > 1)
      stop("x$ltf_cutoff must be a number between 0 and 1 (0 and 1 included) or NULL")
  }
  if(algorithm == "PPR"){
    if (x$damping_factor < 0 | x$damping_factor >= 1)
      stop("x$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }

  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if (correct_topology != TRUE & correct_topology != FALSE)
    stop("correct_topology must be TRUE or FALSE")
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(settings))
    stop("settings should be a list!")
  if(!is.character(settings[[1]]$from) | !is.character(settings[[1]]$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(settings[[1]]$response) | is.null(names(settings[[1]]$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")
  if(!is.character(source_names))
    stop("source_names should be a character vector")
  if(length(source_names) != length(x$source_weights))
    stop("Length of source_names should be the same as length of x$source_weights")
  if(correct_topology == TRUE && !is.null(x$ltf_cutoff))
    warning("Because PPR-ligand-target matrix will be corrected for topology, the proposed cutoff on the ligand-tf matrix will be ignored (x$ltf_cutoff")
  if(correct_topology == TRUE && algorithm != "PPR")
    warning("Topology correction is PPR-specific and makes no sense when the algorithm is not PPR")

  names(x$source_weights) = source_names
  parameters_setting = list(model_name = "query_design", source_weights = x$source_weights)

  if (correct_topology == TRUE){
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = 0, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = TRUE)
  } else {
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = x$ltf_cutoff, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = FALSE)
  }

  output_evaluation = evaluate_model(parameters_setting, lr_network, sig_network, gr_network, settings,calculate_popularity_bias_target_prediction = FALSE,calculate_popularity_bias_ligand_prediction=FALSE,ncitations = ncitations, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links, n_target_bins = 3, ...)

  mean_auroc_target_prediction = output_evaluation$performances_target_prediction$auroc %>% mean() %>% unique()
  mean_aupr_target_prediction = output_evaluation$performances_target_prediction$aupr_corrected %>% mean() %>% unique()

  mean_auroc_ligand_prediction = output_evaluation$performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% .$auroc %>% unique() # unique necessary because possible that two different ligand importance measures result in same maximal performance
  mean_aupr_ligand_prediction = output_evaluation$performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% .$aupr_corrected %>% max() # get aupr corresponding to importance measure resulting in best auroc; take max if two measures have maximal auroc with different aupr.

  return(c(mean_auroc_target_prediction, mean_aupr_target_prediction, mean_auroc_ligand_prediction, mean_aupr_ligand_prediction))
}
#' @title Optimization of objective functions via model-based optimization.
#'
#' @description \code{mlrmbo_optimization} will execute multi-objective model-based optimization of an objective function. The defined surrogate learner here is "kriging".
#'
#' @usage
#' mlrmbo_optimization(run_id,obj_fun,niter,ncores,nstart,additional_arguments)
#'
#' @param run_id Indicate the id of the optimization run.
#' @param obj_fun An objective function as created by the function \code{mlrMBO::makeMultiObjectiveFunction}.
#' @param niter The number of iterations during the optimization process.
#' @param ncores The number of cores on which several parameter settings will be evaluated in parallel.
#' @param nstart The number of different parameter settings used in the begin design.
#' @param additional_arguments A list of named additional arguments that will be passed on the objective function.
#'
#' @return A result object from the function \code{mlrMBO::mbo}. Among other things, this contains the optimal parameter settings, the output corresponding to every input etc.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(mlrMBO)
#' library(parallelMap)
#' additional_arguments_topology_correction = list(source_names = source_weights_df$source %>% unique(), algorithm = "PPR", correct_topology = TRUE,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings = lapply(expression_settings_validation,convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no", cutoff_method = "quantile")
#' nr_datasources = additional_arguments_topology_correction$source_names %>% length()
#'
#' obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",description = "data source weight and hyperparameter optimization: expensive black-box function", fn = model_evaluation_optimization, par.set = makeParamSet( makeNumericVectorParam("source_weights", len = nr_datasources, lower = 0, upper = 1), makeNumericVectorParam("lr_sig_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("damping_factor", len = 1, lower = 0, upper = 0.99)), has.simple.signature = FALSE,n.objectives = 4, noisy = FALSE,minimize = c(FALSE,FALSE,FALSE,FALSE))
#'
#' mlrmbo_optimization = lapply(1,mlrmbo_optimization, obj_fun = obj_fun_multi_topology_correction, niter = 3, ncores = 8, nstart = 100, additional_arguments = additional_arguments_topology_correction)
#'
#' }
#'
#' @export
#'
mlrmbo_optimization = function(run_id,obj_fun,niter,ncores,nstart,additional_arguments){

  requireNamespace("mlrMBO")
  requireNamespace("parallelMap")
  requireNamespace("dplyr")

  # input check

  if (length(run_id) != 1)
    stop("run_id should be a vector of length 1")
  if(!is.function(obj_fun) | !is.list(attributes(obj_fun)$par.set$pars))
    stop("obj_fun should be a function (and generated by mlrMBO::makeMultiObjectiveFunction)")
   if(niter <= 0)
    stop("niter should be a number higher than 0")
  if(ncores <= 0)
    stop("ncores should be a number higher than 0")
  nparams = attributes(obj_fun)$par.set$pars %>% lapply(function(x){x$len}) %>% unlist() %>% sum()
  if(nstart < nparams)
    stop("nstart should be equal or larger than the number of parameters")
  if (!is.list(additional_arguments))
    stop("additional_arguments should be a list!")


  ctrl = makeMBOControl(n.objectives = attributes(obj_fun) %>% .$n.objectives, propose.points = ncores)
  ctrl = setMBOControlMultiObj(ctrl, method = "dib",dib.indicator = "sms")
  ctrl = setMBOControlInfill(ctrl, crit = makeMBOInfillCritDIB(cb.lambda = 2L))
  ctrl = setMBOControlMultiPoint(ctrl, method = "cb")
  ctrl = setMBOControlTermination(ctrl, iters = niter)

  design = generateDesign(n = nstart, par.set = getParamSet(obj_fun))

  configureMlr(on.learner.warning = "quiet", show.learner.output = FALSE)
  parallelStartMulticore(cpus = ncores, show.info = TRUE)

  surr.rf = makeLearner("regr.km", predict.type = "se")

  print(design)
  print(ctrl)

  res = mbo(obj_fun, design = design, learner = surr.rf ,control = ctrl, show.info = TRUE, more.args = additional_arguments)

  parallelStop()
  return(res)
}
#' @title Construct and evaluate a ligand-target model given input parameters with the purpose of hyperparameter optimization.
#'
#' @description \code{model_evaluation_hyperparameter_optimization} will take as input a setting of parameters (hyperparameters), data source weights and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input validation settings (average performance for both target gene prediction and ligand activity prediction, as measured via the auroc and aupr).
#'
#' @usage
#' model_evaluation_hyperparameter_optimization(x, source_weights, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",...)
#'
#' @inheritParams model_evaluation_optimization
#' @param x A list containing the following elements. $lr_sig_hub: hub correction factor for the ligand-signaling network; $gr_hub: hub correction factor for the gene regulatory network; $damping_factor: damping factor in the PPR algorithm if using PPR and optionally $ltf_cutoff: the cutoff on the ligand-tf matrix. For more information about these parameters: see \code{construct_ligand_target_matrix} and \code{apply_hub_correction}.
#' @param source_weights A named numeric vector indicating the weight for every data source.
#' @param ... Additional arguments to \code{make_discrete_ligand_target_matrix}.
#'
#' @return A numeric vector of length 4 containing the average auroc for target gene prediction, average aupr (corrected for TP fraction) for target gene prediction, average auroc for ligand activity prediction and average aupr for ligand activity prediction.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' nr_datasources = source_weights_df$source %>% unique() %>% length()
#' test_input = list("lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
#' source_weights = source_weights_df$weight
#' names(source_weights) = source_weights_df$source
# test_evaluation_optimization = model_evaluation_hyperparameter_optimization(test_input, source_weights, "PPR", TRUE, lr_network, sig_network, gr_network, lapply(expression_settings_validation,convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no")
#' }
#'
#' @export
#'
model_evaluation_hyperparameter_optimization = function(x, source_weights, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",...){

  requireNamespace("dplyr")

  #input check
  if (!is.list(x))
    stop("x should be a list!")
  if (x$lr_sig_hub < 0 | x$lr_sig_hub > 1)
    stop("x$lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (x$gr_hub < 0 | x$gr_hub > 1)
    stop("x$gr_hub must be a number between 0 and 1 (0 and 1 included)")
  if(is.null(x$ltf_cutoff)){
    if(correct_topology == FALSE)
      warning("Did you not forget to give a value to x$ltf_cutoff?")
  } else {
    if (x$ltf_cutoff < 0 | x$ltf_cutoff > 1)
      stop("x$ltf_cutoff must be a number between 0 and 1 (0 and 1 included) or NULL")
  }
  if (!is.numeric(source_weights) | is.null(names(source_weights)))
    stop("source_weights should be a named numeric vector")
  if(algorithm == "PPR"){
    if (x$damping_factor < 0 | x$damping_factor >= 1)
      stop("x$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if (correct_topology != TRUE & correct_topology != FALSE)
    stop("correct_topology must be TRUE or FALSE")
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(settings))
    stop("settings should be a list!")
  if(!is.character(settings[[1]]$from) | !is.character(settings[[1]]$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(settings[[1]]$response) | is.null(names(settings[[1]]$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")
  if(correct_topology == TRUE && !is.null(x$ltf_cutoff))
    warning("Because PPR-ligand-target matrix will be corrected for topology, the proposed cutoff on the ligand-tf matrix will be ignored (x$ltf_cutoff")
  if(correct_topology == TRUE && algorithm != "PPR")
    warning("Topology correction is PPR-specific and makes no sense when the algorithm is not PPR")

  parameters_setting = list(model_name = "query_design", source_weights = source_weights)

  if (correct_topology == TRUE){
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = 0, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = TRUE)
  } else {
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = x$ltf_cutoff, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = FALSE)
  }

  output_evaluation = evaluate_model(parameters_setting, lr_network, sig_network, gr_network, settings,calculate_popularity_bias_target_prediction = FALSE,calculate_popularity_bias_ligand_prediction=FALSE,ncitations = ncitations, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links, n_target_bins = 3, ...)

  mean_auroc_target_prediction = output_evaluation$performances_target_prediction$auroc %>% mean()
  mean_aupr_target_prediction = output_evaluation$performances_target_prediction$aupr_corrected %>% mean()

  mean_auroc_ligand_prediction = output_evaluation$performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% .$auroc %>% unique() # unique necessary because possible that two different ligand importance measures result in same maximal performance
  mean_aupr_ligand_prediction = output_evaluation$performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% .$aupr_corrected %>% unique()

  return(c(mean_auroc_target_prediction, mean_aupr_target_prediction, mean_auroc_ligand_prediction, mean_aupr_ligand_prediction))
}
#' @title Process the output of mlrmbo multi-objective optimization to extract optimal parameter values.
#'
#' @description \code{process_mlrmbo_nichenet_optimization} will process the output of multi-objective mlrmbo optimization. As a result, a list containing the optimal parameter values for model construction will be returned.
#'
#' @usage
#' process_mlrmbo_nichenet_optimization(optimization_results,source_names,parameter_set_index = NULL)
#'
#' @param optimization_results A list generated as output from multi-objective optimization by mlrMBO. Should contain the elements $pareto.front, $pareto.set See \code{mlrmbo_optimization}.
#' @param source_names Character vector containing the names of the data sources. The order of data source names accords to the order of weights in x$source_weights.
#' @param parameter_set_index Number indicating which of the proposed solutions must be selected to extract optimal parameters. If NULL: the solution of which that average z-score of the objective functions outputs is the highest will be selected. Default: NULL.
#'
#' @return A list containing the parameter values leading to maximal performance and thus with the following elements: $source_weight_df, $lr_sig_hub, $gr_hub, $ltf_cutoff, $damping_factor
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(mlrMBO)
#' library(parallelMap)
#' additional_arguments_topology_correction = list(source_names = source_weights_df$source %>% unique(), algorithm = "PPR", correct_topology = TRUE,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings = lapply(expression_settings_validation,convert_expression_settings_evaluation), secondary_targets = FALSE, remove_direct_links = "no", cutoff_method = "quantile")
#' nr_datasources = additional_arguments_topology_correction$source_names %>% length()
#'
#' obj_fun_multi_topology_correction = makeMultiObjectiveFunction(name = "nichenet_optimization",description = "data source weight and hyperparameter optimization: expensive black-box function", fn = model_evaluation_optimization, par.set = makeParamSet( makeNumericVectorParam("source_weights", len = nr_datasources, lower = 0, upper = 1), makeNumericVectorParam("lr_sig_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("gr_hub", len = 1, lower = 0, upper = 1),  makeNumericVectorParam("damping_factor", len = 1, lower = 0, upper = 0.99)), has.simple.signature = FALSE,n.objectives = 4, noisy = FALSE,minimize = c(FALSE,FALSE,FALSE,FALSE))
#'
#' mlrmbo_optimization_result = lapply(1,mlrmbo_optimization, obj_fun = obj_fun_multi_topology_correction, niter = 3, ncores = 8, nstart = 100, additional_arguments = additional_arguments_topology_correction)
#' optimized_parameters = process_mlrmbo_nichenet_optimization(mlrmbo_optimization_result[[1]],additional_arguments_topology_correction$source_names)
#'
#' }
#'
#' @export
#'
process_mlrmbo_nichenet_optimization = function(optimization_results,source_names,parameter_set_index = NULL){

  requireNamespace("dplyr")
  requireNamespace("tibble")

  if(length(optimization_results) == 1){
    optimization_results = optimization_results[[1]]
  }
  # input check
  if (!is.list(optimization_results))
    stop("optimization_results should be a list!")
  if (!is.list(optimization_results$pareto.set))
    stop("optimization_results$pareto.set should be a list! Are you sure you provided the output of mlrMBO::mbo (multi-objective)?")
  if (!is.matrix(optimization_results$pareto.front))
    stop("optimization_results$pareto.front should be a matrix! Are you sure you provided the output of mlrMBO::mbo (multi-objective?")
  if (!is.character(source_names))
    stop("source_names should be a character vector")
  if(!is.numeric(parameter_set_index)  & !is.null(parameter_set_index))
    stop("parameter_set_index should be a number or NULL")

  # winning parameter set
  if(is.null(parameter_set_index)){
    # parameter_set_index = optimization_results$pareto.front %>% tbl_df() %>% mutate(average = apply(.,1,mean), index = seq(nrow(.))) %>% filter(average == max(average)) %>% .$index
    parameter_set_index = optimization_results$pareto.front %>% tbl_df() %>% apply(2,function(x){(x-mean(x))/sd(x)}) %>% tbl_df() %>% mutate(average = apply(.,1,mean), index = seq(nrow(.))) %>% filter(average == max(average)) %>% .$index # take the best parameter setting considering the average of z-scores for each objective function result
  }
  if(parameter_set_index > nrow(optimization_results$pareto.front))
    stop("parameter_set_index may not be a number higher than the total number of proposed solutions")
  parameter_set = optimization_results$pareto.set[[parameter_set_index]]

  # data source weight model parameter
  source_weights = parameter_set$source_weights
  names(source_weights) = source_names

  # "hyperparameters"
  lr_sig_hub = parameter_set$lr_sig_hub
  gr_hub =  parameter_set$gr_hub
  ltf_cutoff =  parameter_set$ltf_cutoff
  damping_factor = parameter_set$damping_factor

  source_weight_df = tibble(source = names(source_weights), weight = source_weights)
  output_optimization = list(source_weight_df = source_weight_df, lr_sig_hub = lr_sig_hub, gr_hub = gr_hub,ltf_cutoff = ltf_cutoff, damping_factor = damping_factor)
  return(output_optimization)

}
#' @title Construct and evaluate a ligand-target model given input parameters with the purpose of parameter optimization for multi-ligand application.
#'
#' @description \code{model_evaluation_optimization_application} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input application settings (average performance for target gene prediction, as measured via the auroc and aupr).
#'
#' @usage
#' model_evaluation_optimization_application(x, source_names, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",classification_algorithm = "lda",...)
#'
#' @inheritParams model_evaluation_optimization
#' @param classification_algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost", "pls" (load "pls" package before!); without: "lda","naive_bayes", "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param ... Additional arguments to \code{evaluate_multi_ligand_target_prediction}.
#'
#' @return A numeric vector of length 2 containing the average auroc and aupr for target gene prediction.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' nr_datasources = source_weights_df$source %>% unique() %>% length()
#' test_input = list("source_weights" = rep(0.5, times = nr_datasources), "lr_sig_hub" = 0.5, "gr_hub" = 0.5, "damping_factor" = 0.5)
# test_evaluation_optimization = model_evaluation_optimization_application(test_input, source_weights_df$source %>% unique(), algorithm = "PPR", TRUE, lr_network, sig_network, gr_network, list(convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries)), secondary_targets = FALSE, remove_direct_links = "no", classification_algorithm = "lda", var_imps = FALSE, cv_number = 5, cv_repeats = 4)
#' }
#'
#' @export
#'
model_evaluation_optimization_application = function(x, source_names, algorithm, correct_topology, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",classification_algorithm = "lda",...){

  requireNamespace("dplyr")

  #input check
  if (!is.list(x))
    stop("x should be a list!")
  if (!is.numeric(x$source_weights))
    stop("x$source_weights should be a numeric vector")
  if (x$lr_sig_hub < 0 | x$lr_sig_hub > 1)
    stop("x$lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (x$gr_hub < 0 | x$gr_hub > 1)
    stop("x$gr_hub must be a number between 0 and 1 (0 and 1 included)")
  if(is.null(x$ltf_cutoff)){
    if(correct_topology == FALSE)
      warning("Did you not forget to give a value to x$ltf_cutoff?")
  } else {
    if (x$ltf_cutoff < 0 | x$ltf_cutoff > 1)
      stop("x$ltf_cutoff must be a number between 0 and 1 (0 and 1 included) or NULL")
  }
  if(algorithm == "PPR"){
    if (x$damping_factor < 0 | x$damping_factor >= 1)
      stop("x$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }

  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if (correct_topology != TRUE & correct_topology != FALSE)
    stop("correct_topology must be TRUE or FALSE")
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(settings))
    stop("settings should be a list!")
  if(!is.character(settings[[1]]$from) | !is.character(settings[[1]]$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(settings[[1]]$response) | is.null(names(settings[[1]]$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")
  if(!is.character(source_names))
    stop("source_names should be a character vector")
  if(length(source_names) != length(x$source_weights))
    stop("Length of source_names should be the same as length of x$source_weights")
  if(correct_topology == TRUE && !is.null(x$ltf_cutoff))
    warning("Because PPR-ligand-target matrix will be corrected for topology, the proposed cutoff on the ligand-tf matrix will be ignored (x$ltf_cutoff")
  if(correct_topology == TRUE && algorithm != "PPR")
    warning("Topology correction is PPR-specific and makes no sense when the algorithm is not PPR")
  if(!is.character(classification_algorithm))
    stop("classification_algorithm should be a character vector of length 1")
  names(x$source_weights) = source_names
  parameters_setting = list(model_name = "query_design", source_weights = x$source_weights)

  if (correct_topology == TRUE){
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = 0, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = TRUE)
  } else {
    parameters_setting = add_hyperparameters_parameter_settings(parameters_setting, lr_sig_hub = x$lr_sig_hub, gr_hub = x$gr_hub, ltf_cutoff = x$ltf_cutoff, algorithm = algorithm,damping_factor = x$damping_factor,correct_topology = FALSE)
  }

  output_evaluation = evaluate_model_application_multi_ligand(parameters_setting, lr_network, sig_network, gr_network, settings, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links, classification_algorithm = classification_algorithm,...)

  mean_auroc_target_prediction = output_evaluation$performances_target_prediction$auroc %>% mean()
  mean_aupr_target_prediction = output_evaluation$performances_target_prediction$aupr_corrected %>% mean()

  return(c(mean_auroc_target_prediction, mean_aupr_target_prediction))
}
#' @title Estimate data source weights of data sources of interest based on leave-one-in and leave-one-out characterization performances.
#'
#' @description \code{estimate_source_weights_characterization} will estimate data source weights of data sources of interest based on a model that was trained to predict weights of data sources based on leave-one-in and leave-one-out characterization performances.
#'
#' @usage
#' estimate_source_weights_characterization(loi_performances,loo_performances,source_weights_df, sources_oi, random_forest =FALSE)
#'
#' @param loi_performances Performances of models in which a particular data source of interest was the only data source in or the ligand-signaling or the gene regulatory network.
#' @param loo_performances Performances of models in which a particular data source of interest was removed from the ligand-signaling or the gene regulatory network before model construction.
#' @param source_weights_df A data frame / tibble containing the weights associated to each individual data source. Sources with higher weights will contribute more to the final model performance (required columns: source, weight). Note that only interactions described by sources included here, will be retained during model construction.
#' @param sources_oi The names of the data sources of which data source weights should be estimated based on leave-one-in and leave-one-out performances.
#' @param random_forest Indicate whether for the regression between leave-one-in + leave-one-out performances and data source weights a random forest model should be trained (TRUE) or a linear model (FALSE). Default: FALSE
#'
#' @return A list containing two elements. $source_weights_df (the input source_weights_df extended by the estimated source_weighs for data sources of interest) and $model (model object of the regression between leave-one-in, leave-one-out performances and data source weights).
#'
#' @importFrom purrr reduce
#' @importFrom randomForest randomForest
#'
#' @examples
#' \dontrun{
#' library(dplyr)
# run characterization loi
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR", damping_factor = 0.2, correct_topology = TRUE)

#' doMC::registerDoMC(cores = 4)
#' job_characterization_loi = parallel::mclapply(weights_settings_loi[1:4], evaluate_model,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings,calculate_popularity_bias_target_prediction = FALSE, calculate_popularity_bias_ligand_prediction = FALSE, ncitations, mc.cores = 4)
#' loi_performances = process_characterization_target_prediction_average(job_characterization_loi)

# run characterization loo
#' weights_settings_loo = prepare_settings_leave_one_out_characterization(lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, source_weights_df)
#' weights_settings_loo = lapply(weights_settings_loo,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR", damping_factor = 0.2, correct_topology = TRUE)

#' doMC::registerDoMC(cores = 4)
#' job_characterization_loo = parallel::mclapply(weights_settings_loo[1:4], evaluate_model,lr_network = lr_network, sig_network = sig_network, gr_network = gr_network, settings,calculate_popularity_bias_target_prediction = FALSE, calculate_popularity_bias_ligand_prediction = FALSE,ncitations,mc.cores = 4)
#' loo_performances = process_characterization_target_prediction_average(job_characterization_loo)

# run the regression
#' sources_oi = c("kegg_cytokines")
#' output = estimate_source_weights_characterization(loi_performances,loo_performances,source_weights_df %>% filter(source != "kegg_cytokines"), sources_oi, random_forest =FALSE)
#' }
#'
#' @export
#'
estimate_source_weights_characterization = function(loi_performances,loo_performances,source_weights_df, sources_oi, random_forest =FALSE){

  requireNamespace("dplyr")
  requireNamespace("tibble")

  #input check
  if(!is.data.frame(loi_performances))
    stop("loi_performances should be a data frame")
  if(!is.character(loi_performances$model_name))
    stop("loi_performances$model_name should be a character vector")
  if(!is.data.frame(loo_performances))
    stop("loo_performances should be a data frame")
  if(!is.character(loo_performances$model_name))
    stop("loo_performances$model_name should be a character vector")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")
  if(!is.character(sources_oi))
    stop("sources_oi should be a character vector")
  if(random_forest != TRUE & random_forest != FALSE)
    stop("random_forest should be TRUE or FALSE")

  loi_performances_train = loi_performances %>% filter((model_name %in% sources_oi) == FALSE)
  loo_performances_train = loo_performances %>% filter((model_name %in% sources_oi) == FALSE)

  loi_performances_test = loi_performances %>% filter(model_name == "complete_model" | (model_name %in% sources_oi))
  loo_performances_test = loo_performances %>% filter(model_name == "complete_model" | (model_name %in% sources_oi))

  output_regression_model = regression_characterization_optimization(loi_performances_train, loo_performances_train, source_weights_df, random_forest = random_forest)

  new_source_weight_df = assign_new_weight(loi_performances_test, loo_performances_test,output_regression_model,source_weights_df)
  return(list(source_weights_df = new_source_weight_df, model = output_regression_model))
}
