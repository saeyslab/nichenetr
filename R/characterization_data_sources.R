#' @title Prepare settings for leave-one-in characterization
#'
#' @description \code{prepare_settings_leave_one_in_characterization} will generate a list of lists containing the data source weights that need to be used for model construction. Every sub-list will contain the data source weights needed to make so called leave-one-in models in which only one ligand-signaling data source is used and all gene regulatory data sources (or vice versa).
#'
#' @usage
#' prepare_settings_leave_one_in_characterization(lr_network, sig_network, gr_network, source_weights_df)
#'
#' @inheritParams construct_weighted_networks
#' @return A list of lists. Every sub-list contains 2 elements: $model_name: the name of the left-in data source; $source_weights: named numeric vector containing the data source weights that will be used for the construction of leave-one-in models.
#'
#' @examples
#' \dontrun{
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' }
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
#' @return A list of lists. Every sub-list contains 2 elements: $model_name: the name of the left-out data source; $source_weights: named numeric vector containing the data source weights that will be used for the construction of leave-one-out models.
#'
#' @examples
#' \dontrun{
# weights_settings_loo = prepare_settings_leave_one_out_characterization(lr_network,sig_network, gr_network, source_weights_df)
#'}
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
  weights_settings_loo = c(weights_settings,weights_settings_lr_sig, weights_settings_gr)
  return(weights_settings_loo)
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
#' \dontrun{
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' }
#' @export
#'
#'
add_hyperparameters_parameter_settings = function(parameter_setting,lr_sig_hub,gr_hub,ltf_cutoff,algorithm,damping_factor,correct_topology){

  # input check
  if (!is.list(parameter_setting))
    stop("parameter_setting must be a list")
  if (lr_sig_hub < 0 | lr_sig_hub > 1)
    stop("lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (gr_hub < 0 | gr_hub > 1)
    stop("gr_hub must be a number between 0 and 1 (0 and 1 included)")

  if(is.null(ltf_cutoff)){
    if( algorithm == "PPR" | algorithm == "SPL" )
      warning("Did you not forget to give a value to ltf_cutoff?")
  } else {
    if (ltf_cutoff < 0 | ltf_cutoff > 1)
      stop("ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(algorithm == "PPR"){
    if (damping_factor < 0 | damping_factor >= 1)
      stop("damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
  if (correct_topology != TRUE & correct_topology != FALSE)
    stop("correct_topology must be TRUE or FALSE")
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
#' evaluate_model(parameters_setting, lr_network, sig_network, gr_network, settings, calculate_popularity_bias_target_prediction,calculate_popularity_bias_ligand_prediction, ncitations = ncitations, secondary_targets = FALSE, remove_direct_links = "no", n_target_bins = 3,...)
#'
#' @inheritParams construct_weighted_networks
#' @inheritParams construct_ligand_target_matrix
#' @param parameters_setting A list containing following elements: $model_name, $source_weights, $lr_sig_hub, $gr_hub, $ltf_cutoff, $algorithm, $damping_factor, $correct_topology. See \code{prepare_settings_leave_one_in_characterization} and \code{add_hyperparameters_parameter_settings}.
#' @param settings A list of lists for which each sub-list contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active ligand(s) or a FALSE.
#' @param calculate_popularity_bias_target_prediction Indicate whether popularity bias in target gene prediction performance should be calculated (TRUE or FALSE).
#' @param calculate_popularity_bias_ligand_prediction Indicate whether popularity bias in ligand activity prediction performance should be calculated (TRUE or FALSE).
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'. Default: ncitations (variable contained in this package). See function \code{get_ncitations_genes} for a function that makes this data frame from current Pubmed information.
#' @param n_target_bins Indicate the number of bins the target genes will be divided in according to popularity. Only relevant when calculate_popularity_bias_target_prediction and/or calculate_popularity_bias_ligand_prediction is/are TRUE.  Default = 3.
#' @param ... Additional arguments to \code{make_discrete_ligand_target_matrix}.
#'
#' @return A list containing following elements: $model_name, $performances_target_prediction, $performances_ligand_prediction, $performances_ligand_prediction_single
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
#' }
#'
#' @export
#'
evaluate_model = function(parameters_setting, lr_network, sig_network, gr_network, settings,calculate_popularity_bias_target_prediction,calculate_popularity_bias_ligand_prediction ,ncitations = ncitations, secondary_targets = FALSE, remove_direct_links = "no", n_target_bins = 3, ...){

  requireNamespace("dplyr")

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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }

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

  if (calculate_popularity_bias_target_prediction != TRUE & calculate_popularity_bias_target_prediction != FALSE)
    stop("calculate_popularity_bias_target_prediction must be TRUE or FALSE")
  if (calculate_popularity_bias_ligand_prediction != TRUE & calculate_popularity_bias_ligand_prediction != FALSE)
    stop("calculate_popularity_bias_ligand_prediction must be TRUE or FALSE")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")
  if(n_target_bins < 0 | n_target_bins > 10)
    stop("n_target_bins should be a number higher than 0. In addition we recommend to keep this number not higher than 10 in order to have a meaningful analysis")
  # construct model
  ligands =  extract_ligands_from_settings(settings)
  output_model_construction = construct_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links)
  model_name = output_model_construction$model_name
  ligand_target_matrix = output_model_construction$model
  # ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)

  ## if in ligand-target matrix: all targets are zero for some ligands
  ligands_zero = ligand_target_matrix %>% colnames() %>% sapply(function(ligand){sum(ligand_target_matrix[,ligand]) == 0}) %>% .[. == TRUE]
  if (length(ligands_zero > 0)){
    noisy_target_scores = runif(nrow(ligand_target_matrix), min = 0, max = min(ligand_target_matrix[ligand_target_matrix>0])) # give ligands not in model a very low noisy random score; why not all 0 --> ties --> problem aupr calculation
    ligand_target_matrix[,names(ligands_zero)] = noisy_target_scores
  }
  # transcriptional response evaluation
  performances_target_prediction = bind_rows(lapply(settings,evaluate_target_prediction, ligand_target_matrix))
  # performances_target_prediction_discrete = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix_discrete))
  # performances_target_prediction = performances_target_prediction %>% full_join(performances_target_prediction_discrete, by = c("setting", "ligand"))
  if (calculate_popularity_bias_target_prediction == TRUE){
    ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)
    performances_target_prediction = performances_target_prediction %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})

    # print(performances_target_prediction) ########################## print

    # ligand-level
    performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances_target_prediction, ncitations)

    # print(performances_ligand_popularity) ########################## print


    ligand_slopes_df = performances_ligand_popularity %>% select(-setting,-ligand,-ncitations) %>% colnames() %>% lapply(.,get_slope_ligand_popularity,performances_ligand_popularity) %>% bind_rows()

    # print(ligand_slopes_df) ########################## print

     # target-level
    performances_target_bins_popularity = evaluate_target_prediction_per_bin(n_target_bins,settings,ligand_target_matrix, ncitations)

    # print(performances_target_bins_popularity) ########################## print


    target_slopes_df = performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) ,method = "all") %>% bind_rows()

    # print(target_slopes_df) ########################## print

    performances_target_bins_popularity = evaluate_target_prediction_per_bin(n_target_bins,settings,ligand_target_matrix_discrete, ncitations)

    # print(performances_target_bins_popularity) ########################## print

    target_slopes_df_discrete = performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}), method = "all") %>% bind_rows()

    # print(target_slopes_df_discrete) ########################## print


    target_slopes_df = bind_rows(target_slopes_df, target_slopes_df_discrete)

    popularity_slopes_target_prediction = inner_join(ligand_slopes_df, target_slopes_df, by = "metric")
  } else {
    popularity_slopes_target_prediction = NULL
  }

  # ligand activity state prediction
  all_ligands = unlist(extract_ligands_from_settings(settings, combination = FALSE))

  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = TRUE)
  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix[, all_ligands]))
  # ligand_importances_discrete = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix_discrete[, all_ligands]))
  # ligand_importances_discrete = ligand_importances_discrete %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})
  # if(sum(is.na(ligand_importances_discrete$fisher_odds)) > 0){
  #   ligand_importances_discrete = ligand_importances_discrete %>% select(-fisher_odds) %>% select(-fisher_pval_log) # because contains too much NA sometimes in leave one in models
  # }

  # settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = FALSE)
  # ligand_importances_glm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix[,all_ligands], algorithm = "glm", cv = FALSE)) %>% rename(glm_imp = importance)

  ligand_importances$pearson[is.na(ligand_importances$pearson)] = 0
  ligand_importances$spearman[is.na(ligand_importances$spearman)] = 0
  ligand_importances$pearson_log_pval[is.na(ligand_importances$pearson_log_pval)] = 0
  ligand_importances$spearman_log_pval[is.na(ligand_importances$spearman_log_pval)] = 0
  ligand_importances$mean_rank_GST_log_pval[is.na(ligand_importances$mean_rank_GST_log_pval)] = 0
  ligand_importances$pearson_log_pval[is.infinite(ligand_importances$pearson_log_pval)] = 10000
  ligand_importances$spearman_log_pval[is.infinite(ligand_importances$spearman_log_pval)] = 10000
  ligand_importances$mean_rank_GST_log_pval[is.infinite(ligand_importances$mean_rank_GST_log_pval)] = 10000

  all_importances = ligand_importances %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})
  # all_importances = full_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand")) %>% full_join(ligand_importances_discrete, by = c("setting","test_ligand", "ligand"))

  # evaluation = suppressWarnings(evaluate_importances_ligand_prediction(all_importances, "median","lda",cv_number = 3, cv_repeats = 20))
  # warning lda here: variables are collinear --> not problematic but logical here
  # performances_ligand_prediction = evaluation$performances
  # all_importances = all_importances %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})
  performances_ligand_prediction_single = all_importances$setting %>% unique() %>% lapply(function(x){x}) %>%
    lapply(wrapper_evaluate_single_importances_ligand_prediction,all_importances) %>%
    bind_rows() %>% inner_join(all_importances %>% distinct(setting,ligand))

  # performances_ligand_prediction_single = evaluate_single_importances_ligand_prediction(all_importances, "median")

  if (calculate_popularity_bias_ligand_prediction == TRUE){
    ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)
    # print(all_importances) ########################## print
    # print(all_importances$test_ligand %>% unique()) ########################## print
    # ligand level
    i_max = round(0.75*length(all_ligands))
    ligand_activity_popularity_bias = lapply(0:i_max,ligand_activity_performance_top_i_removed, all_importances, ncitations) %>% bind_rows()
    slopes_df_ligand = ligand_activity_popularity_bias %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-importance_measure, -popularity_index) %>% colnames() %>% lapply(.,get_ligand_slope_ligand_prediction_popularity ,ligand_activity_popularity_bias %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})) %>% bind_rows()

    # # target level
    performances_target_bins_popularity = evaluate_ligand_prediction_per_bin(3,settings,ligand_target_matrix,ncitations)
    slopes_df_target = performances_target_bins_popularity  %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-importance_measure,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity_ligand_prediction,performances_target_bins_popularity  %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})) %>% bind_rows()
    popularity_slopes_ligand_prediction = inner_join(slopes_df_ligand, slopes_df_target, by = "metric")

    }
  else {
    popularity_slopes_ligand_prediction = NULL
  }
  return(list(model_name = model_name, performances_target_prediction = performances_target_prediction, performances_ligand_prediction_single = performances_ligand_prediction_single, popularity_slopes_target_prediction = popularity_slopes_target_prediction,popularity_slopes_ligand_prediction = popularity_slopes_ligand_prediction,all_importances = all_importances))
    # return(list(model_name = model_name, performances_target_prediction = performances_target_prediction,performances_ligand_prediction = performances_ligand_prediction, performances_ligand_prediction_single = performances_ligand_prediction_single, popularity_slopes_target_prediction = popularity_slopes_target_prediction,popularity_slopes_ligand_prediction = popularity_slopes_ligand_prediction))

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
#' library(dplyr)
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
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
#' library(dplyr)
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
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
  performances_target_prediction_averages = performances_target_prediction %>% group_by(model_name) %>% select(-ligand, -setting) %>% mutate_all(mean) %>% distinct() %>% ungroup()
  return(performances_target_prediction_averages)
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
#' library(dplyr)
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
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
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
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
#' @title Process the output of model evaluation for data source characterization purposes on the popularity bias assessment of target prediction performance
#'
#' @description \code{process_characterization_popularity_slopes_target_prediction} will process output formed by model evaluation to get a data frame containing popularity bias measures in performance of target gene prediction.
#'
#' @usage
#' process_characterization_popularity_slopes_target_prediction(output_characterization)
#'
#' @param output_characterization a list of lists containing the results of evaluation of different models (e.g. after execution of \code{evaluate_model}. Every sublist should contain at least the following elements: $model_name and $performances_target_prediction.
#' @return A data frame containing the popularity bias slopes in target gene prediction performances on every validation dataset for all the models that were evaluated.
#'
#' @examples
#'
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)

#' popularity_slopes_target_prediction_performances = process_characterization_popularity_slopes_target_prediction(output_characterization)
#' }
#'
#' @export
#'
#'
process_characterization_popularity_slopes_target_prediction = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$popularity_slopes_target_prediction))
    stop("output_characterization[[1]]$popularity_slopes_target_prediction should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  popularity_slopes_target_prediction = output_characterization %>% lapply(function(x){
    performances = x$popularity_slopes_target_prediction %>% mutate(model_name = x$model_name)
    return(performances)
  }) %>% bind_rows()
  return(popularity_slopes_target_prediction)
}
#' @title Process the output of model evaluation for data source characterization purposes on the popularity bias assessment of ligand activity performance
#'
#' @description \code{process_characterization_popularity_slopes_ligand_prediction} will process output formed by model evaluation to get a data frame containing popularity bias measures in performance of ligand activity prediction.
#'
#' @usage
#' process_characterization_popularity_slopes_ligand_prediction(output_characterization)
#'
#' @param output_characterization a list of lists containing the results of evaluation of different models (e.g. after execution of \code{evaluate_model}. Every sublist should contain at least the following elements: $model_name and $performances_target_prediction.
#' @return A data frame containing the popularity bias slopes in ligand activity prediction performances.
#'
#' @examples
#'
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation, convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)

#' popularity_slopes_ligand_prediction_performances = process_characterization_popularity_slopes_ligand_prediction(output_characterization)
#' }
#'
#' @export
#'
#'
process_characterization_popularity_slopes_ligand_prediction = function(output_characterization){
  # input check
  if (!is.list(output_characterization))
    stop("output_characterization should be a list!")
  if (!is.data.frame(output_characterization[[1]]$popularity_slopes_ligand_prediction))
    stop("output_characterization[[1]]$popularity_slopes_ligand_prediction should be a data frame")
  if (!is.character(output_characterization[[1]]$model_name))
    stop("output_characterization[[1]]$model_name should be a character vector")

  requireNamespace("dplyr")

  popularity_slopes_ligand_prediction = output_characterization %>% lapply(function(x){
    performances = x$popularity_slopes_ligand_prediction %>% mutate(model_name = x$model_name)
    return(performances)
  }) %>% bind_rows()
  return(popularity_slopes_ligand_prediction)
}
#' @title Prepare settings for one-vs-one characterization
#'
#' @description \code{prepare_settings_one_vs_one_characterization} will generate a list of lists containing the data source weights that need to be used for model construction. Every sub-list will contain the data source weights needed to make so called one-vs-one models in which only one ligand-signaling data source and only one gene regulatory data source is used. It is also possible to construct one-vs-one-vs-one models by keeping 1 ligand-receptor, 1 signaling and 1 gene regulatory data source.
#'
#' @usage
#' prepare_settings_one_vs_one_characterization(lr_network, sig_network, gr_network, lr_network_separate = FALSE)
#'
#' @inheritParams construct_weighted_networks
#' @param lr_network_separate Indicate whether the one-vs-one models should contain 1 ligand-receptor, 1 signaling and 1 gene regulatory network source (TRUE) or just 1 ligand-signaling combined and 1 gene regulatory source (FALSE). Default: FALSE.
#' @return A list of lists. Every sub-list contains following elements: $model_name; $source_lr_sig: the name of the left-in ligand-signaling data source (or $source_lr and $source_sig if lr_network_separate is TRUE); $source_gr: the name of the left-in gene regulatory data source and $source_weights: named numeric vector containing the data source weights that will be used for the construction of one-vs-one models.
#'
#' @examples
#' \dontrun{
#' weights_settings_ovo = prepare_settings_one_vs_one_characterization(lr_network,sig_network, gr_network)
#' weights_settings_ovo_lr_separate = prepare_settings_one_vs_one_characterization(lr_network,sig_network, gr_network, lr_network_separate = TRUE)
#' }
#' @export
#'
#'
prepare_settings_one_vs_one_characterization = function(lr_network, sig_network, gr_network, lr_network_separate = FALSE){
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
   if(lr_network_separate != TRUE & lr_network_separate != FALSE)
    stop("lr_network_separate should be TRUE or FALSE")

  requireNamespace("dplyr")
  if (lr_network_separate == FALSE){
    lr_sig_sources = c(lr_network$source, sig_network$source) %>% unique()
    gr_sources = unique(gr_network$source)

    source_possibilities_df = expand.grid(lr_sig_sources, gr_sources) %>% as_tibble() %>% rename(lr_sig_source = Var1, gr_source = Var2) %>% mutate(lr_sig_source = as.character(lr_sig_source), gr_source = as.character(gr_source))

    source_weight_settings_list = lapply(seq(nrow(source_possibilities_df)), function(i, source_possibilities_df){
      row_oi = source_possibilities_df[i,]
      novel_weights = rep(0,times=length(c(lr_sig_sources,gr_sources)))
      names(novel_weights) = c(lr_sig_sources,gr_sources)
      novel_weights[row_oi$lr_sig_source]=1
      novel_weights[row_oi$gr_source]=1
      lr_sig_source = row_oi$lr_sig_source
      gr_source = row_oi$gr_source
      model_name = paste(lr_sig_source, gr_source,sep = "-")
      return(list(model_name = model_name, lr_sig_source = lr_sig_source, gr_source = gr_source, source_weights = novel_weights))
    }, source_possibilities_df)
  } else {
    lr_sources = unique(lr_network$source)
    sig_sources = unique(sig_network$source)
    gr_sources = unique(gr_network$source)

    source_possibilities_df = expand.grid(lr_sources,sig_sources, gr_sources) %>% as_tibble() %>% rename(lr_source = Var1, sig_source = Var2, gr_source = Var3) %>% mutate(lr_source = as.character(lr_source), sig_source = as.character(sig_source), gr_source = as.character(gr_source))

    source_weight_settings_list = lapply(seq(nrow(source_possibilities_df)), function(i, source_possibilities_df){
      row_oi = source_possibilities_df[i,]
      novel_weights = rep(0,times=length(c(lr_sources, sig_sources, gr_sources)))
      names(novel_weights) = c(lr_sources, sig_sources, gr_sources)
      novel_weights[row_oi$lr_source]=1
      novel_weights[row_oi$sig_source]=1
      novel_weights[row_oi$gr_source]=1
      lr_source = row_oi$lr_source
      sig_source = row_oi$sig_source
      gr_source = row_oi$gr_source
      model_name = paste(lr_source, sig_source, gr_source,sep = "-")
      return(list(model_name = model_name, lr_source = lr_source, sig_source = sig_source, gr_source = gr_source, source_weights = novel_weights))
    }, source_possibilities_df)
  }
   return(source_weight_settings_list)
}
#' @title Construct and evaluate a ligand-target model given input parameters (for application purposes).
#'
#' @description \code{evaluate_model_application} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input application settings (only target gene prediction).
#'
#' @usage
#' evaluate_model_application(parameters_setting, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no", ...)
#'
#' @inheritParams evaluate_model
#'
#' @return A list containing following elements: $model_name, $performances_target_prediction.
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model_application,lr_network,sig_network, gr_network,settings, mc.cores = 3)
#' }
#'
#' @export
#'
evaluate_model_application = function(parameters_setting, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no", ...){

  requireNamespace("dplyr")

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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
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

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")

  # construct model
  ligands =  extract_ligands_from_settings(settings)
  output_model_construction = construct_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links)
  model_name = output_model_construction$model_name
  ligand_target_matrix = output_model_construction$model
  ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)

  # transcriptional response evaluation
  performances_target_prediction = bind_rows(lapply(settings,evaluate_target_prediction, ligand_target_matrix))
  performances_target_prediction_discrete = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix_discrete))
  performances_target_prediction = performances_target_prediction %>% full_join(performances_target_prediction_discrete, by = c("setting", "ligand"))

  return(list(model_name = model_name, performances_target_prediction = performances_target_prediction))
}

#' @title Construct a ligand-target model given input parameters.
#'
#' @description \code{construct_model} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix.
#'
#' @usage
#' construct_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = FALSE, remove_direct_links = "no")
#'
#' @inheritParams evaluate_model
#' @param ligands List of ligands for which the model should be constructed
#'
#' @return A list containing following elements: $model_name and $model.
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' ligands =  extract_ligands_from_settings(settings)
#' models_characterization = parallel::mclapply(weights_settings_loi[1:3],construct_model,lr_network,sig_network, gr_network,ligands, mc.cores = 3)
#' }
#'
#' @export
#'
construct_model = function(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = FALSE, remove_direct_links = "no"){

  requireNamespace("dplyr")


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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")

  if (parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
  if (parameters_setting$correct_topology != TRUE & parameters_setting$correct_topology != FALSE)
    stop("parameters_setting$correct_topology must be TRUE or FALSE")
  if(parameters_setting$correct_topology == TRUE && parameters_setting$algorithm != "PPR")
    warning("Topology correction is PPR-specific and makes no sense when the algorithm is not PPR")


  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(ligands))
    stop("ligands should be a list!")

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")

  # read in parameters
  model_name = parameters_setting$model_name

  source_weights = parameters_setting$source_weights
  source_weights_df = tibble::tibble(source = names(source_weights), weight = source_weights)

  lr_sig_hub = parameters_setting$lr_sig_hub
  gr_hub = parameters_setting$gr_hub
  ltf_cutoff = parameters_setting$ltf_cutoff
  algorithm = parameters_setting$algorithm
  damping_factor = parameters_setting$damping_factor
  correct_topology = parameters_setting$correct_topology

  # construct weighted networks
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df) %>% apply_hub_corrections(lr_sig_hub, gr_hub)

  # extract ligands and construct ligand-target matrix
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks,
                                                        ligands = ligands,
                                                        ltf_cutoff = ltf_cutoff,
                                                        algorithm = algorithm,
                                                        damping_factor =  damping_factor,
                                                        secondary_targets = secondary_targets,
                                                        remove_direct_links = remove_direct_links)
  if (correct_topology == TRUE & algorithm == "PPR"){
    ligand_target_matrix = correct_topology_ppr(ligand_target_matrix, weighted_networks)
  }
  return(list(model_name = model_name, model = ligand_target_matrix))
}
#' @title Assess the influence of an individual data source on ligand-target probability scores
#'
#' @description \code{assess_influence_source} will assess the influence of an individual data source on ligand-target probability scores (or rankings of these). Possible output: the ligand-target matrices of the complete model vs the leave-one-out model in which the data source of interest was left out; or a list indicating which target genes for every ligand of interest are affected the most by leaving out the data source of interest.
#'
#' @usage
#' assess_influence_source(source, lr_network, sig_network, gr_network, source_weights_df, ligands, rankings = FALSE, matrix_output = FALSE,  secondary_targets = FALSE, remove_direct_links = "no", ...)
#'
#' @inheritParams construct_weighted_networks
#' @inheritParams construct_model
#' @param source Name of the data source that will be left out to assess its influence.
#' @param ... Argumentes for the function \code{add_hyperparameters_parameter_settings}
#' @param  rankings Indicate whether the output of the models should be the ranking of target gene probability scores (TRUE; top target gene rank = 1) or the scores themselves (FALSE). Default: FALSE.
#' @param matrix_output Indicate whether the output should be the 2 ligand-target matrices (complete model and leave-one-out model) (TRUE) or a listing of genes of which the ligand-target scores/rankings were influenced the most (FALSE). Default: FALSE.
#' @return If matrix_output == TRUE: A list of sublists; every sublist contains the elements $model_name and $model: the constructed ligand-target matrix. If matrix_output == FALSE: A list of sublist: every sublist contains; $ligand: name of the ligand tested; $targets_higher: sorted vector of ligand-target scores or rankings of target that score higher in the complete model compared to the leave-one-out model; targets_lower: sorted vector of ligand-target scores or rankings of target that score lower in the complete model compared to the leave-one-out model.
#'
#' @examples
#' \dontrun{
#' ligands =  extract_ligands_from_settings(expression_settings_validation[1:4])
#' output = assess_influence_source("ontogenet", lr_network,sig_network, gr_network, source_weights_df, ligands,lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' }
#' @export
#'
#'
assess_influence_source = function(source, lr_network, sig_network, gr_network, source_weights_df, ligands, rankings = FALSE, matrix_output = FALSE,  secondary_targets = FALSE, remove_direct_links = "no", ...){
  # input check
  all_sources = unique(c(lr_network$source,sig_network$source,gr_network$source))
  if (!is.character(source) | (source %in% all_sources) == FALSE)
    stop("source must be a character vector and exist in one of the networks")
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")
  if (rankings != TRUE & rankings != FALSE)
    stop("rankings must be TRUE or FALSE")
  if (matrix_output != TRUE & matrix_output != FALSE)
    stop("matrix_output must be TRUE or FALSE")

  if(!is.list(ligands))
    stop("ligands should be a list")
  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")

  requireNamespace("dplyr")

  ## make weight settings for the complete and leave-one-out model

  all_sources = source_weights_df$source %>% unique()
  all_weights = source_weights_df %>% filter(source %in% all_sources) %>% .$weight
  names(all_weights) = all_sources

  weights_settings = list()

  # complete
  weights_settings[[1]] = list("complete_model",all_weights)
  names(weights_settings[[1]]) = c("model_name","source_weights")

  # leave one out
  all_novel_weights = all_weights
  all_novel_weights[source] = 0

  weights_settings[[2]] = list(source,all_novel_weights)
  names(weights_settings[[2]]) = c("model_name","source_weights")

  ## Add other model construction (hyper)parameters to the model settings

  weights_settings = lapply(weights_settings,add_hyperparameters_parameter_settings, ...)

  ## Construct models
  models_output = lapply(weights_settings,construct_model,lr_network,sig_network, gr_network,ligands, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links)

  ## because different data sources used, target genes of the models can be different
  ## therefore we will keep only target genes present in both models

  intersecting_rownames = intersect(rownames(models_output[[1]]$model), rownames(models_output[[2]]$model))
  models_output[[1]]$model = models_output[[1]]$model %>% .[rownames(.) %in% intersecting_rownames,]
  models_output[[2]]$model = models_output[[2]]$model %>% .[rownames(.) %in% intersecting_rownames,]

  if(rankings == TRUE){
    models_output[[1]]$model = models_output[[1]]$model %>% apply(2,rank_desc)
    models_output[[2]]$model = models_output[[2]]$model %>% apply(2,rank_desc)
  }
  if (matrix_output == TRUE){
    return(models_output)
  } else {
      lt_diff_source_complete = models_output[[1]]$model - models_output[[2]]$model
      affected_targets_output = get_affected_targets_output(lt_diff_source_complete, rankings)
    }
    return(affected_targets_output)
}
#' @title Construct and evaluate a ligand-target model given input parameters (for application purposes + multi-ligand predictive model).
#'
#' @description \code{evaluate_model_application_multi_ligand} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input application settings (only target gene prediction; multi-ligand classification).
#'
#' @usage
#' evaluate_model_application_multi_ligand(parameters_setting, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",classification_algorithm = "lda", ...)
#'
#' @inheritParams evaluate_model_application
#' @param classification_algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost", "pls" (load "pls" package before!); without: "lda","naive_bayes", "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param ... Optional arguments to \code{evaluate_multi_ligand_target_prediction}.
#'
#' @return A list containing following elements: $model_name, $performances_target_prediction.
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_model_application_multi_ligand,lr_network,sig_network, gr_network,settings, classification_algorithm = "lda", var_imps = FALSE, cv_number = 5, cv_repeats = 4, parallel = TRUE, mc.cores = 3)
#'
#' }
#'
#' @export
#'
evaluate_model_application_multi_ligand = function(parameters_setting, lr_network, sig_network, gr_network, settings, secondary_targets = FALSE, remove_direct_links = "no",classification_algorithm = "lda", ...){

  requireNamespace("dplyr")

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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
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

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")

  # construct model
  ligands =  extract_ligands_from_settings(settings)
  output_model_construction = construct_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links)
  model_name = output_model_construction$model_name
  ligand_target_matrix = output_model_construction$model

  # transcriptional response evaluation
  performances_target_prediction = bind_rows(lapply(settings,wrapper_evaluate_multi_ligand_target_prediction, ligand_target_matrix,algorithm = classification_algorithm, ...))

  return(list(model_name = model_name, performances_target_prediction = performances_target_prediction))
}
#' @title Construct and evaluate a randomised ligand-target model given input parameters.
#'
#' @description \code{evaluate_random_model} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix and evaluate its performance on input validation settings (both target gene prediction and ligand activity prediction) after randomisation of the networks by edge swapping.
#'
#' @usage
#' evaluate_random_model(parameters_setting, lr_network, sig_network, gr_network, settings, calculate_popularity_bias_target_prediction,calculate_popularity_bias_ligand_prediction, ncitations = ncitations, secondary_targets = FALSE, remove_direct_links = "no", n_target_bins = 3,...)
#'
#' @inheritParams evaluate_model
#'
#' @return A list containing following elements: $model_name, $performances_target_prediction, $performances_ligand_prediction, $performances_ligand_prediction_single
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' output_characterization = parallel::mclapply(weights_settings_loi[1:3],evaluate_random_model,lr_network,sig_network, gr_network,settings,calculate_popularity_bias_target_prediction = TRUE, calculate_popularity_bias_ligand_prediction = TRUE, ncitations, mc.cores = 3)
#' }
#'
#' @export
#'
evaluate_random_model = function(parameters_setting, lr_network, sig_network, gr_network, settings,calculate_popularity_bias_target_prediction,calculate_popularity_bias_ligand_prediction ,ncitations = ncitations, secondary_targets = FALSE, remove_direct_links = "no", n_target_bins = 3, ...){

  requireNamespace("dplyr")

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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }

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

  if (calculate_popularity_bias_target_prediction != TRUE & calculate_popularity_bias_target_prediction != FALSE)
    stop("calculate_popularity_bias_target_prediction must be TRUE or FALSE")
  if (calculate_popularity_bias_ligand_prediction != TRUE & calculate_popularity_bias_ligand_prediction != FALSE)
    stop("calculate_popularity_bias_ligand_prediction must be TRUE or FALSE")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")
  if(n_target_bins < 0 | n_target_bins > 10)
    stop("n_target_bins should be a number higher than 0. In addition we recommend to keep this number not higher than 10 in order to have a meaningful analysis")
  # construct model
  ligands =  extract_ligands_from_settings(settings)
  output_model_construction = construct_random_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = secondary_targets, remove_direct_links = remove_direct_links)
  model_name = output_model_construction$model_name
  ligand_target_matrix = output_model_construction$model
  ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)

  # transcriptional response evaluation
  performances_target_prediction = bind_rows(lapply(settings,evaluate_target_prediction, ligand_target_matrix))
  performances_target_prediction_discrete = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix_discrete))
  performances_target_prediction = performances_target_prediction %>% full_join(performances_target_prediction_discrete, by = c("setting", "ligand"))
  if (calculate_popularity_bias_target_prediction == TRUE){
    performances_target_prediction = performances_target_prediction %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})

    # print(performances_target_prediction) ########################## print

    # ligand-level
    performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances_target_prediction, ncitations)

    # print(performances_ligand_popularity) ########################## print


    ligand_slopes_df = performances_ligand_popularity %>% select(-setting,-ligand,-ncitations) %>% colnames() %>% lapply(.,get_slope_ligand_popularity,performances_ligand_popularity) %>% bind_rows()

    # print(ligand_slopes_df) ########################## print

    # target-level
    performances_target_bins_popularity = evaluate_target_prediction_per_bin(n_target_bins,settings,ligand_target_matrix, ncitations)

    # print(performances_target_bins_popularity) ########################## print


    target_slopes_df = performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) ,method = "all") %>% bind_rows()

    # print(target_slopes_df) ########################## print


    performances_target_bins_popularity = evaluate_target_prediction_per_bin(n_target_bins,settings,ligand_target_matrix_discrete, ncitations)

    # print(performances_target_bins_popularity) ########################## print

    target_slopes_df_discrete = performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}), method = "all") %>% bind_rows()

    # print(target_slopes_df_discrete) ########################## print


    target_slopes_df = bind_rows(target_slopes_df, target_slopes_df_discrete)

    popularity_slopes_target_prediction = inner_join(ligand_slopes_df, target_slopes_df, by = "metric")
  } else {
    popularity_slopes_target_prediction = NULL
  }

  # ligand activity state prediction
  all_ligands = unlist(extract_ligands_from_settings(settings, combination = FALSE))

  # print(all_ligands) ########################## print

  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = TRUE)
  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix[, all_ligands]))

  # print(ligand_importances) ########################## print

  ligand_importances_discrete = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix_discrete[, all_ligands]))

  # print(ligand_importances_discrete) ########################## print


  # ligand_importances_discrete = ligand_importances_discrete %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})
  # if(sum(is.na(ligand_importances_discrete$fisher_odds)) > 0){
  #   ligand_importances_discrete = ligand_importances_discrete %>% select(-fisher_odds) %>% select(-fisher_pval_log) # because contains too much NA sometimes in leave one in models
  # }

  # print(ligand_importances_discrete) ########################## print


  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = FALSE)
  ligand_importances_glm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix[,all_ligands], algorithm = "glm", cv = FALSE)) %>% rename(glm_imp = importance)

  # print(ligand_importances_glm) ########################## print

  all_importances = full_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand")) %>% full_join(ligand_importances_discrete, by = c("setting","test_ligand", "ligand"))
  # all_importances = inner_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand"))

  # evaluation = suppressWarnings(evaluate_importances_ligand_prediction(all_importances, "median","lda",cv_number = 3, cv_repeats = 20))
  # warning lda here: variables are collinear --> not problematic but logical here
  # performances_ligand_prediction = evaluation$performances
  all_importances = all_importances %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})

  performances_ligand_prediction_single = evaluate_single_importances_ligand_prediction(all_importances, "median")

  if (calculate_popularity_bias_ligand_prediction == TRUE){
    # print(all_importances) ########################## print
    # print(all_importances$test_ligand %>% unique()) ########################## print
    # ligand level
    i_max = round(0.75*length(all_ligands))
    ligand_activity_popularity_bias = lapply(0:i_max,ligand_activity_performance_top_i_removed, all_importances, ncitations) %>% bind_rows()
    slopes_df_ligand = ligand_activity_popularity_bias %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-importance_measure, -popularity_index) %>% colnames() %>% lapply(.,get_ligand_slope_ligand_prediction_popularity ,ligand_activity_popularity_bias %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})) %>% bind_rows()

    # # target level
    performances_target_bins_popularity = evaluate_ligand_prediction_per_bin(3,settings,ligand_target_matrix,ncitations)
    slopes_df_target = performances_target_bins_popularity  %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) %>% select(-importance_measure,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity_ligand_prediction,performances_target_bins_popularity  %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})) %>% bind_rows()
    popularity_slopes_ligand_prediction = inner_join(slopes_df_ligand, slopes_df_target, by = "metric")

  }
  else {
    popularity_slopes_ligand_prediction = NULL
  }
  return(list(model_name = model_name, performances_target_prediction = performances_target_prediction, performances_ligand_prediction_single = performances_ligand_prediction_single, popularity_slopes_target_prediction = popularity_slopes_target_prediction,popularity_slopes_ligand_prediction = popularity_slopes_ligand_prediction))
  # return(list(model_name = model_name, performances_target_prediction = performances_target_prediction,performances_ligand_prediction = performances_ligand_prediction, performances_ligand_prediction_single = performances_ligand_prediction_single, popularity_slopes_target_prediction = popularity_slopes_target_prediction,popularity_slopes_ligand_prediction = popularity_slopes_ligand_prediction))

}
#' @title Construct a randomised ligand-target model given input parameters.
#'
#' @description \code{construct_random_model} will take as input a setting of parameters (data source weights and hyperparameters) and layer-specific networks to construct a ligand-target matrix after randomization by edge swapping.
#'
#' @usage
#' construct_random_model(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = FALSE, remove_direct_links = "no")
#'
#' @inheritParams evaluate_model
#' @param ligands List of ligands for which the model should be constructed
#'
#' @return A list containing following elements: $model_name and $model.
#'
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:4], convert_expression_settings_evaluation)
#' weights_settings_loi = prepare_settings_leave_one_in_characterization(lr_network,sig_network, gr_network, source_weights_df)
#' weights_settings_loi = lapply(weights_settings_loi,add_hyperparameters_parameter_settings, lr_sig_hub = 0.25,gr_hub = 0.5,ltf_cutoff = 0,algorithm = "PPR",damping_factor = 0.8,correct_topology = TRUE)
#' doMC::registerDoMC(cores = 8)
#' ligands =  extract_ligands_from_settings(settings)
#' models_characterization = parallel::mclapply(weights_settings_loi[1:3],construct_random_model,lr_network,sig_network, gr_network,ligands, mc.cores = 3)
#' }
#'
#' @export
#'
construct_random_model = function(parameters_setting, lr_network, sig_network, gr_network, ligands, secondary_targets = FALSE, remove_direct_links = "no"){

  requireNamespace("dplyr")


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
  if(is.null(parameters_setting$ltf_cutoff)){
    if( parameters_setting$algorithm == "PPR" | parameters_setting$algorithm == "SPL" )
      warning("Did you not forget to give a value to parameters_setting$ltf_cutoff?")
  } else {
    if (parameters_setting$ltf_cutoff < 0 | parameters_setting$ltf_cutoff > 1)
      stop("parameters_setting$ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }
  if (parameters_setting$algorithm != "PPR" & parameters_setting$algorithm != "SPL" & parameters_setting$algorithm != "direct")
    stop("parameters_setting$algorithm must be 'PPR' or 'SPL' or 'direct'")

  if (parameters_setting$algorithm == "PPR"){
    if (parameters_setting$damping_factor < 0 | parameters_setting$damping_factor >= 1)
      stop("parameters_setting$damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
  if (parameters_setting$correct_topology != TRUE & parameters_setting$correct_topology != FALSE)
    stop("parameters_setting$correct_topology must be TRUE or FALSE")
  if(parameters_setting$correct_topology == TRUE && parameters_setting$algorithm != "PPR")
    warning("Topology correction is PPR-specific and makes no sense when the algorithm is not PPR")


  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.list(ligands))
    stop("ligands should be a list!")

  if (secondary_targets != TRUE & secondary_targets != FALSE)
    stop("secondary_targets must be TRUE or FALSE")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be  'no' or 'ligand' or 'ligand-receptor'")

  # read in parameters
  model_name = parameters_setting$model_name

  source_weights = parameters_setting$source_weights
  source_weights_df = tibble::tibble(source = names(source_weights), weight = source_weights)

  lr_sig_hub = parameters_setting$lr_sig_hub
  gr_hub = parameters_setting$gr_hub
  ltf_cutoff = parameters_setting$ltf_cutoff
  algorithm = parameters_setting$algorithm
  damping_factor = parameters_setting$damping_factor
  correct_topology = parameters_setting$correct_topology

  # construct weighted networks
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df) %>% apply_hub_corrections(lr_sig_hub, gr_hub)
  weighted_networks$lr_sig = weighted_networks$lr_sig %>% randomize_network(output_weighted = TRUE)
  weighted_networks$gr = weighted_networks$gr %>% randomize_network(output_weighted = TRUE)
  # extract ligands and construct ligand-target matrix
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks = weighted_networks,
                                                        ligands = ligands,
                                                        ltf_cutoff = ltf_cutoff,
                                                        algorithm = algorithm,
                                                        damping_factor =  damping_factor,
                                                        secondary_targets = secondary_targets,
                                                        remove_direct_links = remove_direct_links)
  if (correct_topology == TRUE & algorithm == "PPR"){
    ligand_target_matrix = correct_topology_ppr(ligand_target_matrix, weighted_networks)
  }
  return(list(model_name = model_name, model = ligand_target_matrix))
}




















