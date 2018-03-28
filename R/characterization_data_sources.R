#' @title Prepare settings for leave-one-in characterization
#'
#' @description \code{prepare_settings_leave_one_in_characterization} will generate a list of lists containing the data source weights that need to be used for model construction. Every sub-list will contain the data source weights needed to make so called leave-one-in models in which only one ligand-signaling data source is used and all gene regulatory data sources (or vice versa).
#'
#' @usage
#' prepare_settings_leave_one_in_characterization(lr_network, sig_network, gr_network, source_weights_df)
#'
#' @inheritParams construct_weighted_networks
#' @return A list of lists. Every sub-list contains 2 elements: $source: the name of the left-in data source; $source_weights: named numeric vector containing the data source weights that will be used for the construction of leave-one-in models.
#'
#' @examples
#'
#'
#' @export
#'
#'
prepare_settings_leave_one_in_characterization = function(lr_network, sig_network, gr_network, source_weights_df){
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")

  requireNamespace("dplyr")

  lr_sig_sources = source_weights_df %>% filter(source %in% c(lr_network$source %>% unique(), sig_network$source %>% unique())) %>% .$source
  gr_sources = source_weights_df %>% filter(source %in% unique(gr_network$source)) %>% .$source

  weights_lr_sig = source_weights_df %>% filter(source %in% lr_sig_sources) %>% .$weight
  names(weights_lr_sig) = lr_sig_sources

  weights_gr = source_weights_df %>% filter(source %in% gr_sources) %>% .$weight
  names(weights_gr) = gr_sources

  weights_settings = list()
  weights_settings[[1]] = list("complete_model",c(weights_lr_sig, weights_gr))
  names(weights_settings[[1]]) = c("model_name","source_weights")


  weights_settings_lr_sig  = list()
  for(i in 1:length(weights_lr_sig)){
    novel_weights = rep(0,times=length(weights_lr_sig))
    names(novel_weights) = names(weights_lr_sig)
    novel_weights[i]=1
    weights_settings_lr_sig[[i]]= list(names(weights_lr_sig[i]),c(novel_weights,weights_gr))
    names(weights_settings_lr_sig[[i]]) = c("model_name","source_weights")
  }

  weights_settings_gr  = list()
  for(i in 1:length(weights_gr)){
    novel_weights = rep(0,times=length(weights_gr))
    names(novel_weights) = names(weights_gr)
    novel_weights[i]= 1
    weights_settings_gr[[i]]= list(names(weights_gr[i]),c(weights_lr_sig,novel_weights))
    names(weights_settings_gr[[i]]) = c("model_name","source_weights")

  }
  weights_settings_loi = c(weights_settings,weights_settings_lr_sig, weights_settings_gr)
  return(weights_settings_loi)
}
#' @title Prepare settings for leave-one-out characterization
#'
#' @description \code{prepare_settings_leave_one_out_characterization} will generate a list of lists containing the data source weights that need to be used for model construction. Every sub-list will contain the data source weights needed to make so called leave-one-out models in one ligand-signaling data source is or gene regulatory data source is left out.
#'
#' @usage
#' prepare_settings_leave_one_out_characterization(lr_network, sig_network, gr_network, source_weights_df)
#'
#' @inheritParams construct_weighted_networks
#' @return A list of lists. Every sub-list contains 2 elements: $source: the name of the left-in data source; $source_weights: named numeric vector containing the data source weights that will be used for the construction of leave-one-in models.
#'
#' @examples
#'
# weights_settings_loo = prepare_settings_leave_one_out_characterization(lr_network,sig_network, gr_network, source_weights_df)
#'
#' @export
#'
#'
prepare_settings_leave_one_out_characterization = function(lr_network, sig_network, gr_network, source_weights_df){
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")

  requireNamespace("dplyr")

  lr_sig_sources = source_weights_df %>% filter(source %in% c(lr_network$source %>% unique(), sig_network$source %>% unique())) %>% .$source
  gr_sources = source_weights_df %>% filter(source %in% unique(gr_network$source)) %>% .$source

  weights_lr_sig = source_weights_df %>% filter(source %in% lr_sig_sources) %>% .$weight
  names(weights_lr_sig) = lr_sig_sources

  weights_gr = source_weights_df %>% filter(source %in% gr_sources) %>% .$weight
  names(weights_gr) = gr_sources

  weights_settings = list()
  weights_settings[[1]] = list("complete_model",c(weights_lr_sig, weights_gr))
  names(weights_settings[[1]]) = c("model_name","source_weights")


  weights_settings_lr_sig  = list()
  for(i in 1:length(weights_lr_sig)){
    novel_weights = weights_lr_sig
    novel_weights[i]=0
    weights_settings_lr_sig[[i]]= list(names(weights_lr_sig[i]),c(novel_weights,weights_gr))
    names(weights_settings_lr_sig[[i]]) = c("model_name","source_weights")
  }

  weights_settings_gr  = list()
  for(i in 1:length(weights_gr)){
    novel_weights = weights_gr
    novel_weights[i]= 0
    weights_settings_gr[[i]]= list(names(weights_gr[i]),c(weights_lr_sig,novel_weights))
    names(weights_settings_gr[[i]]) = c("model_name","source_weights")

  }
  weights_settings_loi = c(weights_settings,weights_settings_lr_sig, weights_settings_gr)
  return(weights_settings_loi)
}
#' @title Add hyperparameters to existing parameter settings
#'
#' @description \code{add_hyperparameters_parameter_settings} will generate a list of lists containing the parameter values that need to be used for model construction.
#'
#' @usage
#' add_hyperparameters_parameter_settings(parameter_setting, lr_sig_hub, gr_hub, ltf_cutoff, algorithm, damping_factor, correct_topology)
#'
#' @param parameter_setting A list of lists. Sublists contains parameters (like data source weights) and novel sublists will be created by this function by adding additional parameters.
#' @param correct_topology This parameter indicates whether the PPR-constructed ligand-target matrix will be subtracted by a PR-constructed target matrix. TRUE or FALSE.
#' @inheritParams construct_ligand_target_matrix
#' @inheritParams apply_hub_corrections
#'
#' @return A list of lists. Every sub-list contains parameter values for a different parameter.
#'
#' @examples
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#'
#' @export
#'
#'
add_hyperparameters_parameter_settings = function(parameter_setting,lr_sig_hub,gr_hub,ltf_cutoff,algorithm,damping_factor,correct_topology){

  # input check
  if (!is.list(parameter_setting))
    stop("parameter_setting must be a list")

  requireNamespace("dplyr")

  parameter_setting$lr_sig_hub = lr_sig_hub
  parameter_setting$gr_hub = gr_hub
  parameter_setting$ltf_cutoff = ltf_cutoff
  parameter_setting$algorithm = algorithm
  parameter_setting$damping_factor = damping_factor
  parameter_setting$correct_topology = correct_topology
  return(parameter_setting)
}
#' @title Construct and evaluate a ligand-target model given input parameters.
#'
#' @description \code{evaluate_model} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input validation settings (both target gene prediction and ligand activity prediction).
#'
#' @usage
#' evaluate_model(parameters_setting, lr_network, sig_network, gr_network, settings,...)
#'
#' @inheritParams construct_weighted_networks
#' @param parameters_setting A list containing following elements: $model_name, $source_weights, $lr_sig_hub, $gr_hub, $ltf_cutoff, $algorithm, $damping_factor, $correct_topology. See \code{prepare_settings_leave_one_in_characterization} and \code{add_hyperparameters_parameter_settings}.
#' @param settings A list of lists for which each sub-list contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active ligand(s) or a FALSE.
#' @param ... Additional arguments to \code{construct_ligand_target_matrix} such as "secondary_targets" and "remove_direct_links".
#'
#' @return A list containing following elements: $model_name, $performances_target_prediction, $performances_ligand_prediction, $performances_ligand_prediction_single
#'
#' @import tibble tibble
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' }
#'
#' @export
#'
evaluate_model = function(parameters_setting, lr_network, sig_network, gr_network, settings,...){
  # input check
  if (!is.list(parameters_setting))
    stop("parameters_setting should be a list!")
  if (!is.character(parameters_setting$model_name))
    stop("parameters_setting$model_name should be a character vector")
  if (!is.numeric(parameters_setting$source_weights) | is.null(names(parameters_setting$source_weights)))
    stop("parameters_setting$source_weights should be a named numeric vector")
  if (parameters_setting$lr_sig_hub < 0 | parameters_setting$lr_sig_hub > 1)
    stop("parameters_setting$lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (parameters_setting$gr_hub < 0 | parameters_setting$gr_hub > 1)
    stop("parameters_setting$gr_hub must be a number between 0 and 1 (0 and 1 included)")
  if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
    stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")
  if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
    stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  if (parameters_setting$correct_topology != TRUE & parameters_setting$correct_topology != FALSE)
    stop("parameters_setting$correct_topology must be TRUE or FALSE")

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

  requireNamespace("dplyr")

  # read in parameters
  model_name = parameters_setting$model_name

  source_weights = parameters_setting$source_weights
  source_weights_df = tibble(source = names(source_weights), weight = source_weights)

  lr_sig_hub = parameters_setting$lr_sig_hub
  gr_hub = parameters_setting$gr_hub
  ltf_cutoff = parameters_setting$ltf_cutoff
  algorithm = parameters_setting$algorithm
  damping_factor = parameters_setting$damping_factor
  correct_topology = parameters_setting$correct_topology

  # construct weighted networks
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df) %>% apply_hub_corrections(lr_sig_hub, gr_hub)

  # extract ligands and construct ligand-target matrix
  ligands = extract_ligands_from_settings(settings)
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks,
                                                        ligands = ligands,
                                                        ltf_cutoff = ltf_cutoff,
                                                        algorithm = algorithm,
                                                        damping_factor =  damping_factor,
                                                        ...)
  if (correct_topology == TRUE & algorithm == "PPR"){
    ligand_target_matrix = correct_topology_ppr(ligand_target_matrix, weighted_networks)
  }

  # transcriptional response evaluation
  performances_target_prediction = bind_rows(lapply(settings,evaluate_target_prediction, ligand_target_matrix))
  # performances_target_prediction_discrete = bind_rows(lapply(settings,evaluate_target_prediction, ligand_target_matrix %>% make_discrete_ligand_target_matrix(cutoff_method = "quantile")))

  # ligand activity state prediction
  all_ligands = unlist(extract_ligands_from_settings(settings, combination = FALSE))
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = TRUE)
  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix[, all_ligands]))
  # ligand_importances_discrete = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix[, all_ligands] %>% make_discrete_ligand_target_matrix(cutoff_method = "quantile"))

  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = FALSE)
  ligand_importances_glm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix[,all_ligands], algorithm = "glm", cv = FALSE)) %>% rename(glm_imp = importance)

  all_importances = inner_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand"))
  # all_importances = inner_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand")) %>% inner_join(ligand_importances_discrete, by = c("setting","test_ligand", "ligand"))

  evaluation = evaluate_importances_ligand_prediction(all_importances, "median","lda",cv_number = 3, cv_repeats = 20)
  performances_ligand_prediction = evaluation$performances

  performances_ligand_prediction_single = evaluate_single_importances_ligand_prediction(all_importances, "median")

    return(list(model_name = model_name, performances_target_prediction = performances_target_prediction,performances_ligand_prediction = performances_ligand_prediction, performances_ligand_prediction_single = performances_ligand_prediction_single ))

}
#' @title Process the output of model evaluation for data source characterization purposes on the target prediction performance
#'
#' @description \code{process_characterization_target_prediction} will process output formed by model evaluation to get a data frame containing performance measures in target gene prediction.
#'
#' @usage
#' process_characterization_target_prediction(output_characterization)
#'
#' @param output_characterization a list of lists containing the results of evaluation of different models (e.g. after execution of \code{evaluate_model}. Every sublist should contain at least the following elements: $model_name and $performances_target_prediction.
#' @return A data frame containing the target gene prediction performances on every validation dataset for all the models that were evaluated.
#'
#' @examples
#'
#' \dontrun{
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' target_prediction_performances = process_characterization_target_prediction(output_characterization)
#' }
#'
#' @export
#'
#'
process_characterization_target_prediction = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$performances_target_prediction))
    stop("output_characterization[[1]]$performances_target_prediction should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  performances_target_prediction = output_characterization %>% lapply(function(x){
    performances = x$performances_target_prediction %>% mutate(model_name = x$model_name)
    return(performances)
  }) %>% bind_rows()
  return(performances_target_prediction)
}
#' @title Process the output of model evaluation for data source characterization purposes on the target prediction performance (average)
#'
#' @description \code{output_characterization} will process output formed by model evaluation to get a data frame containing performance measures in target gene prediction (averaged over all validation datasets).
#'
#' @usage
#' process_characterization_target_prediction_average(output_characterization)
#'
#' @inheritParams  process_characterization_target_prediction
#' @return A data frame containing the target gene prediction performances averaged over all validation datasets for all the models that were evaluated.
#'
#' @examples
#'
#' \dontrun{
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' target_prediction_performances = process_characterization_target_prediction_average(output_characterization)
#' }
#'
#' @export
#'
process_characterization_target_prediction_average = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$performances_target_prediction))
    stop("output_characterization[[1]]$performances_target_prediction should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  performances_target_prediction = process_characterization_target_prediction(output_characterization)
  performances_target_prediction_averages = performances_target_prediction %>% group_by(model_name) %>% select(-ligand, -setting) %>% mutate_all(mean) %>% distinct()
  return(performances_target_prediction_averages)
}
#' @title Process the output of model evaluation for data source characterization purposes on the ligand prediction performance
#'
#' @description \code{process_characterization_ligand_prediction} will process output formed by model evaluation to get a data frame containing performance measures in ligand prediction.
#'
#' @usage
#' process_characterization_ligand_prediction(output_characterization)
#'
#' @inheritParams  process_characterization_target_prediction
#' @return A data frame containing the ligand activity prediction performance for all the models that were evaluated.
#'
#' @examples
#'
#' \dontrun{
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' ligand_prediction_performances = process_characterization_ligand_prediction(output_characterization)
#' }
#'
#' @export
#'
process_characterization_ligand_prediction = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$performances_ligand_prediction))
    stop("output_characterization[[1]]$performances_ligand_prediction should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  performances_ligand_prediction = output_characterization %>% lapply(function(x){
    performances = x$performances_ligand_prediction %>% mutate(model_name = x$model_name)
    return(performances)
  }) %>% bind_rows()
  performances_ligand_prediction = performances_ligand_prediction %>% select(-Resample) %>% group_by(model_name) %>% mutate_all(mean) %>% distinct()
  return(performances_ligand_prediction)
}
#' @title Process the output of model evaluation for data source characterization purposes on the ligand prediction performance (for every importance score individually)
#'
#' @description \code{process_characterization_ligand_prediction_single_measures} will process output formed by model evaluation to get a data frame containing performance measures in ligand prediction for every importance score individually.
#'
#' @usage
#' process_characterization_ligand_prediction_single_measures(output_characterization)
#'
#' @inheritParams  process_characterization_target_prediction
#' @return A data frame containing the ligand activity prediction performance of every ligand importance measure for all the models that were evaluated.
#'
#' @examples
#'
#' \dontrun{
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' ligand_prediction_performances = process_characterization_ligand_prediction_single_measures(output_characterization)
#' }
#'
#' @export
process_characterization_ligand_prediction_single_measures = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$performances_ligand_prediction_single))
    stop("output_characterization[[1]]$performances_ligand_prediction_single should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  performances_ligand_prediction_single = output_characterization %>% lapply(function(x){
    performances = x$performances_ligand_prediction_single %>% mutate(model_name = x$model_name)
    return(performances)
  }) %>% bind_rows()
  return(performances_ligand_prediction_single)
}
