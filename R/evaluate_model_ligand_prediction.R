#' @title Convert settings to correct settings format for ligand prediction.
#'
#' @description \code{convert_settings_ligand_prediction} Converts settings to correct settings format for ligand activity prediction. In this prediction problem, ligands (out of a set of possibly active ligands) will be ranked based on feature importance scores. The format can be made suited for: 1) validation of ligand activity state prediction by calculating individual feature importane scores or 2) feature importance based on models with embedded feature importance determination; applications in which ligands need to be scores based on their possible upstream activity: 3) by calculating individual feature importane scores or 4) feature importance based on models with embedded feature importance determination.
#'
#' @usage
#' convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = TRUE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest.
#' @param all_ligands A character vector of possible ligands that will be considered for the ligand activity state prediction.
#' @param validation TRUE if seetings need to be prepared for validation of ligand activity state predictions (this implies that the true active ligand of a setting is known); FALSE for application purposes when the true active ligand(s) is/are not known.
#' @param single TRUE if feature importance scores for ligands will be calculated by looking at ligans individually. FALSE if the goal is to calculate the feature importance scores via sophisticated classification algorithms like random forest.

#' @return A list with following elements: $name, $ligand: name of active ligand(s) (only if validation is TRUE), $from (ligand(s) that will be tested for activity prediction), $response
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
#' ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE))
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = TRUE, single = TRUE)
#' }
#' @export
#'
#'
convert_settings_ligand_prediction = function(settings,all_ligands,validation = TRUE, single = TRUE){

  # input check
  if(!is.list(settings))
    stop("settings should be a list")
  if(!is.character(all_ligands))
    stop("all_ligands should be a character vector")
  if(!is.logical(validation) | length(validation) != 1)
    stop("validation should be TRUE or FALSE")
  if(!is.logical(single) | length(single) != 1)
    stop("single should be TRUE or FALSE")

  requireNamespace("dplyr")

  new_settings = list()
  if (validation == TRUE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_ligands)){
        test_ligand = all_ligands[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_single_validation(setting,test_ligand))
      }
    }
  } else if (validation == TRUE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_multi_validation(setting,all_ligands))
    }
  } else if (validation == FALSE && single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_ligands)){
        test_ligand = all_ligands[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_single_application(setting,test_ligand))
      }
    }
  } else if (validation == FALSE && single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_multi_application(setting,all_ligands))
    }
  }
  return(new_settings %>% unlist(recursive = FALSE))
}
#' @title Get ligand importances based on target gene prediction performance of single ligands.
#'
#' @description \code{get_single_ligand_importances} Get ligand importance measures for ligands based on how well a single, individual, ligand can predict an observed response. Assess how well every ligand of interest is able to predict the observed transcriptional response in a particular dataset, according to the ligand-target model. It can be assumed that the ligand that best predicts the observed response, is more likely to be the true ligand.
#'
#' @usage
#' get_single_ligand_importances(setting,ligand_target_matrix, ligands_position = "cols", known = TRUE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @inheritParams evaluate_target_prediction
#'
#' @return A data.frame with for each ligand - data set combination, classification evaluation metrics indicating how well the query ligand predicts the response in the particular dataset. Evaluation metrics are the same as in \code{\link{evaluate_target_prediction}}. In addition to the metrics, the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' print(head(ligand_importances))
#' }
#' @export
#'
get_single_ligand_importances = function(setting,ligand_target_matrix, ligands_position = "cols", known = TRUE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  metrics = evaluate_target_prediction(setting, ligand_target_matrix, ligands_position)
  metrics = metrics %>% rename(test_ligand = ligand)
  if (known == TRUE){
    true_ligand = setting$ligand
    metrics_meta = metrics %>% select(setting,test_ligand) %>% bind_cols(tibble(ligand = true_ligand))
    metrics = inner_join(metrics_meta, metrics, by = c("setting","test_ligand"))
  }
  return(metrics)
}
#' @title Get ligand importances from a multi-ligand classfication model.
#'
#' @description \code{get_multi_ligand_importances} A classificiation algorithm chosen by the user is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features) in order to predict the observed response in a particular dataset. Variable importance scores that indicate for each ligand the importance for response prediction, are extracted. It can be assumed that ligands with higher variable importance scores are more likely to be a true active ligand.
#'
#' @usage
#' get_multi_ligand_importances(setting,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost"; without: "lda","naive_bayes","pls"(because bug in current version of pls package), "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#' @param continuous Indicate whether during training of the model, model training and evaluation should be done on class probabilities or discrete class labels. For huge class imbalance, we recommend setting this value to TRUE. Default: TRUE.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @param filter_genes Indicate whether 50 per cent of the genes that are the least variable in ligand-target scores should be removed in order to reduce the training of the model. Default: FALSE.
#'
#' @return A data.frame with for each ligand - data set combination, feature importance scores indicating how important the query ligand is for the prediction of the response in the particular dataset, when prediction is done via a trained classification model with all possible ligands as input. In addition to the importance score(s), the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances_glm = dplyr::bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances,ligand_target_matrix, algorithm = "glm"))
#' print(head(ligand_importances_glm))
#' }
#' @export
#'
get_multi_ligand_importances = function(setting,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    ligand_target_matrix = filter_genes_ligand_target_matrix(ligand_target_matrix,ligands_position)
  }

  setting_name = setting$name
  output = evaluate_multi_ligand_target_prediction(setting, ligand_target_matrix, ligands_position,algorithm, var_imps = TRUE, cv, cv_number, cv_repeats, parallel, n_cores, ignore_errors, continuous)
  metrics = output$var_imps
  metrics = metrics %>% mutate(setting = setting_name) %>% rename(test_ligand = feature)

  if (known == TRUE){
    true_ligand = setting$ligand
    metrics = metrics %>% mutate(ligand = true_ligand)
    metrics = metrics %>% select(setting, test_ligand, ligand, importance)
    return(metrics)
  }
  metrics = metrics %>% select(setting, test_ligand, importance)
  return(metrics)

}
#' @title Evaluation of ligand activity prediction based on ligand importance scores.
#'
#' @description \code{evaluate_importances_ligand_prediction} Evaluate how well a trained model of ligand importance scores is able to predict the true activity state of a ligand. For this it is assumed, that ligand importance measures for truely active ligands will be higher than for non-active ligands. A classificiation algorithm chosen by the user is trained to construct one model based on the ligand importance scores of all ligands of interest (ligands importance scores are considered as features). Several classification evaluation metrics for the prediction are calculated and variable importance scores can be extracted to rank the different importance measures in order of importance for ligand activity state prediction.
#'
#' @usage
#' evaluate_importances_ligand_prediction(importances, normalization, algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4,ignore_errors = FALSE)
#'
#' @param importances A data frame containing at least folowing variables: $setting, $test_ligand, $ligand and one or more feature importance scores. $test_ligand denotes the name of a possibly active ligand, $ligand the name of the truely active ligand.
#' @param normalization Way of normalization of the importance measures: "mean" (classifcal z-score) or "median" (modified z-score)
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost"; without: "lda","naive_bayes","pls"(because bug in current version of pls package), "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param var_imps Indicate whether in addition to classification evaluation performances, variable importances should be calculated. Default: TRUE.
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#'
#' @return A list with the following elements. $performances: data frame containing classification evaluation measure for classification on the test folds during training via cross-validation; $performances_training: data frame containing classification evaluation measures for classification of the final model (discrete class assignments) on the complete data set (performance can be severly optimistic due to overfitting!); $performance_training_continuous: data frame containing classification evaluation measures for classification of the final model (class probability scores) on the complete data set (performance can be severly optimistic due to overfitting!) $var_imps: data frame containing the variable importances of the different ligands (embbed importance score for some classification algorithms, otherwise just the auroc); $prediction_response_df: data frame containing for each ligand-setting combination the ligand importance scores for the individual importance scores, the complete model of importance scores and the ligand activity as well (TRUE or FALSE); $model: the caret model object that can be used on new importance scores to predict the ligand activity state.
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @importFrom limma wilcoxGST
#' @import caret
#' @importFrom purrr safely
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")
#' print(head(evaluation))
#' }
#' @export
#'
evaluate_importances_ligand_prediction = function(importances, normalization, algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE){
  if (!is.data.frame(importances))
    stop("importances must be a data frame")
  if(!is.character(importances$setting) | !is.character(importances$test_ligand) | !is.character(importances$ligand))
    stop("importances$setting, importances$test_ligand and importances$ligand should be character vectors")
  if(normalization != "mean" & normalization != "median")
    stop("normalization should be 'mean' or 'median'")
   if(!is.character(algorithm))
    stop("algorithm should be a character vector")
  if(!is.logical(var_imps) | length(var_imps) > 1)
    stop("var_imps should be a logical vector: TRUE or FALSE")
  if(!is.logical(cv) | length(cv) > 1)
    stop("cv should be a logical vector: TRUE or FALSE")
  if(!is.numeric(cv_number) | length(cv_number) > 1)
    stop("cv_number should be a numeric vector of length 1")
  if(!is.numeric(cv_repeats) | length(cv_repeats) > 1)
    stop("cv_repeats should be a numeric vector of length 1")
  if(!is.logical(parallel) | length(parallel) > 1)
    stop("parallel should be a logical vector: TRUE or FALSE")
  if(!is.numeric(n_cores) | length(n_cores) > 1)
    stop("n_cores should be a numeric vector of length 1")
  if(!is.logical(ignore_errors) | length(ignore_errors) > 1)
    stop("ignore_errors should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

#  importances = importances %>% tidyr::drop_na()
  added = is_ligand_active(importances)
  importances = importances %>% mutate(class = added)

  if (normalization == "mean"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-ligand,-test_ligand,-class) %>% mutate_all(funs(scaling_zscore)) %>% ungroup() %>% select(-setting)
  } else if (normalization == "median"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-ligand,-test_ligand,-class) %>% mutate_all(funs(scaling_modified_zscore)) %>% ungroup() %>% select(-setting)
  }

  response_vector = importances$class %>% make.names() %>% as.factor()
  train_data = normalized_importances %>% mutate(obs = response_vector) %>% data.frame()

  output = wrapper_caret_classification(train_data,algorithm,TRUE,var_imps,cv,cv_number,cv_repeats,parallel,n_cores,prediction_response_df = bind_cols(importances %>% select(setting,ligand,test_ligand,class), normalized_importances),ignore_errors,return_model = TRUE)
  return(output)
}
#' @title Evaluation of ligand activity prediction performance of single ligand importance scores: aggregate all datasets.
#'
#' @description \code{evaluate_single_importances_ligand_prediction} Evaluate how well a single ligand importance score is able to predict the true activity state of a ligand. For this it is assumed, that ligand importance measures for truely active ligands will be higher than for non-active ligands. Several classification evaluation metrics for the prediction are calculated and variable importance scores can be extracted to rank the different importance measures in order of importance for ligand activity state prediction.
#'
#' @usage
#' evaluate_single_importances_ligand_prediction(importances,normalization)
#'
#' @param importances A data frame containing at least folowing variables: $setting, $test_ligand, $ligand and one or more feature importance scores. $test_ligand denotes the name of a possibly active ligand, $ligand the name of the truely active ligand.
#' @param normalization Way of normalization of the importance measures: "mean" (classifcal z-score) or "median" (modified z-score) or "no" (use unnormalized feature importance scores - only recommended when evaluating ligand activity prediction on individual datasets)
#'
#' @return A data frame containing classification evaluation measures for the ligand activity state prediction single, individual feature importance measures.
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @importFrom limma wilcoxGST
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' evaluation = evaluate_single_importances_ligand_prediction(ligand_importances,normalization = "median")
#' print(head(evaluation))
#' }
#' @export
#'
evaluate_single_importances_ligand_prediction = function(importances,normalization){
  if (!is.data.frame(importances))
    stop("importances must be a data frame")
  if(!is.character(importances$setting) | !is.character(importances$test_ligand) | !is.character(importances$ligand))
    stop("importances$setting, importances$test_ligand and importances$ligand should be character vectors")
  if(normalization != "mean" & normalization != "median" & normalization != "no")
    stop("normalization should be 'mean' or 'median' or 'no'")

  requireNamespace("dplyr")
  importances0 = importances %>% select(-setting,-ligand,-test_ligand)
#  importances = importances %>% tidyr::drop_na()
  added = is_ligand_active(importances)

  if (nrow(importances) == 0){
    performances = lapply(importances, classification_evaluation_continuous_pred, added, iregulon = FALSE)
    output = tibble(importance_measure = names(performances))
    performances = bind_rows(performances)
    return(bind_cols(output,performances))
  }

  importances = importances %>% select_if(.predicate = function(x) {
    sum(is.na(x)) == 0
  })

  if (normalization == "mean"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-ligand,-test_ligand) %>% mutate_all(funs(scaling_zscore)) %>% ungroup() %>% select(-setting)
  } else if (normalization == "median"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-ligand,-test_ligand) %>% mutate_all(funs(scaling_modified_zscore)) %>% ungroup() %>% select(-setting)
  } else if (normalization == "no") {
    normalized_importances = importances %>% select(-c(setting,test_ligand,ligand))
  }

  performances = lapply(normalized_importances, classification_evaluation_continuous_pred, added, iregulon = FALSE)
  output = tibble(importance_measure = names(performances))
  performances = bind_rows(performances)
  return(bind_cols(output,performances))
}
#' @title Prediction of ligand activity prediction by a model trained on ligand importance scores.
#'
#' @description \code{model_based_ligand_activity_prediction} Predict the activity state of a ligand based on a classification model that was trained to predict ligand activity state based on ligand importance scores.
#'
#' @usage
#' model_based_ligand_activity_prediction(importances, model, normalization)
#'
#' @param model A model object of a classification object as e.g. generated via caret.
#' @param importances A data frame containing at least folowing variables: $setting, $test_ligand, $ligand and one or more feature importance scores. $test_ligand denotes the name of a possibly active ligand, $ligand the name of the truely active ligand.
#' @param normalization Way of normalization of the importance measures: "mean" (classifcal z-score) or "median" (modified z-score)
#'
#' @return A data frame containing the ligand importance scores and the probabilities that according to the trained model, the ligands are active based on their importance scores.
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")
#'
#' settings = lapply(expression_settings_validation[5:10],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = FALSE, single = TRUE)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix, known = FALSE))
#' activity_predictions = model_based_ligand_activity_prediction(ligand_importances, evaluation$model,"median")
#' print(head(activity_predictions))
#' }
#'
#' @export
#'
model_based_ligand_activity_prediction = function(importances, model, normalization){
  if (!is.list(model))
    stop("model must be a list, derived as model object from model training (e.g. via the caret package)")
  if(model$finalModel$problemType != "Classification" & model$finalModel$problemType != "Regression")
    stop("model should be model object (derived from model training)")
  if (!is.data.frame(importances))
    stop("importances must be a data frame")
  if(!is.character(importances$setting) | !is.character(importances$test_ligand))
    stop("importances$setting and importances$test_ligand should be character vectors")
  if(normalization != "mean" & normalization != "median")
    stop("normalization should be 'mean' or 'median'")

  requireNamespace("dplyr")

#  importances = importances %>% tidyr::drop_na()

  if (normalization == "mean"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-test_ligand) %>% mutate_all(funs(scaling_zscore)) %>% ungroup() %>% select(-setting)
  } else if (normalization == "median"){
    normalized_importances = importances %>% group_by(setting) %>% dplyr::select(-test_ligand) %>% mutate_all(funs(scaling_modified_zscore)) %>% ungroup() %>% select(-setting)
  }

  final_model_predictions = predict(model,newdata = normalized_importances, type = "prob")
  final_model_predictions = final_model_predictions %>% as_tibble() %>% mutate(active = TRUE. > FALSE.) %>% select(-FALSE.) %>% rename(model = TRUE.)
  return(bind_cols(importances,final_model_predictions) %>% as_tibble())

}
#' @title Get ligand importances from a multi-ligand trained random forest model.
#'
#' @description \code{get_multi_ligand_rf_importances} A random forest is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features) in order to predict the observed response in a particular dataset. Variable importance scores that indicate for each ligand the importance for response prediction, are extracted. It can be assumed that ligands with higher variable importance scores are more likely to be a true active ligand.
#'
#' @usage
#' get_multi_ligand_rf_importances(setting,ligand_target_matrix, ligands_position = "cols", ntrees = 1000, mtry = 2, continuous = TRUE, known = TRUE, filter_genes = FALSE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param ntrees Indicate the number of trees used in the random forest algorithm. The more trees, the longer model training takes, but the more robust the extraced importance scores will be. Default: 1000. Recommended for robustness to have till 10000 trees.
#' @param mtry n**(1/mtry) features of the n features will be sampled at each split during the training of the random forest algorithm. Default: 2 (square root).
#' @param continuous Indicate whether during training of the model, model training and evaluation should be done on class probabilities or discrete class labels. For huge class imbalance, we recommend setting this value to TRUE. Default: TRUE.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @param filter_genes Indicate whether 50 per cent of the genes that are the least variable in ligand-target scores should be removed in order to reduce the training of the model. Default: FALSE.
#'
#' @return A data.frame with for each ligand - data set combination, feature importance scores indicating how important the query ligand is for the prediction of the response in the particular dataset, when prediction is done via a trained classification model with all possible ligands as input. In addition to the importance score(s), the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @importFrom randomForest randomForest importance
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances_rf = dplyr::bind_rows(lapply(settings_ligand_pred, get_multi_ligand_rf_importances,ligand_target_matrix, ntrees = 100, mtry = 2))
#' print(head(ligand_importances_rf))
#' }
#'
#' @export
#'
get_multi_ligand_rf_importances = function(setting,ligand_target_matrix, ligands_position = "cols", ntrees = 1000, mtry = 2,  continuous = TRUE, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")
  if(ntrees <= 1)
    stop("ntrees should be higher than 1")
  if(mtry <= 1)
    stop("mtry should be higher than 1")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    ligand_target_matrix = filter_genes_ligand_target_matrix(ligand_target_matrix,ligands_position)
  }

  setting_name = setting$name
  ligands_oi = setting$from

  if (ligands_position == "cols"){
    if(sum((ligands_oi %in% colnames(ligand_target_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_target_matrix")
    prediction_matrix = ligand_target_matrix[,ligands_oi]
    target_genes = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    if(sum((ligands_oi %in% rownames(ligand_target_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_target_matrix")
    prediction_matrix = ligand_target_matrix[ligands_oi,] %>% t()
    target_genes = colnames(ligand_target_matrix)
  }

  response_vector = setting$response
  response_df = tibble(gene = names(response_vector), response = response_vector %>% make.names() %>% as.factor())

  prediction_df = prediction_matrix %>% data.frame() %>% as_tibble()

  if(is.double(prediction_matrix) == FALSE){
    convert_categorical_factor = function(x){
      x = x %>% make.names() %>% as.factor()
    }
    prediction_df = prediction_df %>% mutate_all(funs(convert_categorical_factor))
  }

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  train_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()

  rf_model = randomForest::randomForest(y = train_data$obs,
                                        x = train_data[,-(which(colnames(train_data) == "obs"))],
                                        ntree = ntrees,
                                        mtry = ncol(train_data[,-(which(colnames(train_data) == "obs"))])**(1/mtry) %>% ceiling(),
                                        importance = TRUE
  )

  metrics = randomForest::importance(rf_model) %>% data.frame() %>% tibble::rownames_to_column("test_ligand") %>% as_tibble() %>% mutate(setting = setting_name)

  if (known == TRUE){
    true_ligand = setting$ligand
    metrics = metrics %>% mutate(ligand = true_ligand)
    metrics = metrics %>% select(setting, test_ligand, ligand, MeanDecreaseAccuracy, MeanDecreaseGini)
    return(metrics)
  }
  metrics = metrics %>% select(setting, test_ligand, MeanDecreaseAccuracy, MeanDecreaseGini)
  return(metrics)


}
#' @title Get ligand importances based on target gene value prediction performance of single ligands (regression).
#'
#' @description \code{get_single_ligand_importances_regression} Get ligand importance measures for ligands based on how well a single, individual, ligand can predict an observed response. Assess how well every ligand of interest is able to predict the observed transcriptional response in a particular dataset, according to the ligand-target model. It can be assumed that the ligand that best predicts the observed response, is more likely to be the true ligand. Response: continuous values associated to a gene, e.g. a log fold change value.
#'
#' @usage
#' get_single_ligand_importances_regression(setting,ligand_target_matrix, ligands_position = "cols", known = TRUE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @inheritParams evaluate_target_prediction_regression
#'
#' @return A data.frame with for each ligand - data set combination, regression model fit metrics indicating how well the query ligand predicts the response in the particular dataset. Evaluation metrics are the same as in \code{\link{evaluate_target_prediction_regression}}. In addition to the metrics, the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances_regression,ligand_target_matrix))
#' print(head(ligand_importances))
#' }
#' @export
#'
get_single_ligand_importances_regression = function(setting,ligand_target_matrix, ligands_position = "cols", known = TRUE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  metrics = evaluate_target_prediction_regression(setting, ligand_target_matrix, ligands_position)
  metrics = metrics %>% rename(test_ligand = ligand)
  if (known == TRUE){
    true_ligand = setting$ligand
    metrics_meta = metrics %>% select(setting,test_ligand) %>% bind_cols(tibble(ligand = true_ligand))
    metrics = inner_join(metrics_meta, metrics, by = c("setting","test_ligand"))
  }
  return(metrics)
}
#' @title Get ligand importances from a multi-ligand regression model.
#'
#' @description \code{get_multi_ligand_importances_regression} A regression algorithm chosen by the user is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features) in order to predict the observed response in a particular dataset (respone: e.g. absolute value of log fold change). Variable importance scores that indicate for each ligand the importance for response prediction, are extracted. It can be assumed that ligands with higher variable importance scores are more likely to be a true active ligand.
#'
#' @usage
#' get_multi_ligand_importances_regression(setting,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, known = TRUE, filter_genes = FALSE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost"; without: "lda","naive_bayes","pls"(because bug in current version of pls package), "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @param filter_genes Indicate whether 50 per cent of the genes that are the least variable in ligand-target scores should be removed in order to reduce the training of the model. Default: FALSE.
#'
#' @return A data.frame with for each ligand - data set combination, feature importance scores indicating how important the query ligand is for the prediction of the response in the particular dataset, when prediction is done via a trained regression model with all possible ligands as input. In addition to the importance score(s), the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances_lm = dplyr::bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances_regression,ligand_target_matrix, algorithm = "lm"))
#' print(head(ligand_importances_lm))
#' }
#' @export
#'
get_multi_ligand_importances_regression = function(setting,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    ligand_target_matrix = filter_genes_ligand_target_matrix(ligand_target_matrix,ligands_position)
  }

  setting_name = setting$name
  output = evaluate_multi_ligand_target_prediction_regression(setting, ligand_target_matrix, ligands_position,algorithm, var_imps = TRUE, cv, cv_number, cv_repeats, parallel, n_cores, ignore_errors)
  metrics = output$var_imps
  metrics = metrics %>% mutate(setting = setting_name) %>% rename(test_ligand = feature)

  if (known == TRUE){
    true_ligand = setting$ligand
    metrics = metrics %>% mutate(ligand = true_ligand)
    metrics = metrics %>% select(setting, test_ligand, ligand, importance)
    return(metrics)
  }
  metrics = metrics %>% select(setting, test_ligand, importance)
  return(metrics)

}
#' @title Get ligand importances from a multi-ligand trained random forest regression model.
#'
#' @description \code{get_multi_ligand_rf_importances_regression} A random forest is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features) in order to predict the observed response in a particular dataset (response: e.g. absolute values of log fold change). Variable importance scores that indicate for each ligand the importance for response prediction, are extracted. It can be assumed that ligands with higher variable importance scores are more likely to be a true active ligand.
#'
#' @usage
#' get_multi_ligand_rf_importances_regression(setting,ligand_target_matrix, ligands_position = "cols", ntrees = 1000, mtry = 2, known = TRUE, filter_genes = FALSE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param ntrees Indicate the number of trees used in the random forest algorithm. The more trees, the longer model training takes, but the more robust the extraced importance scores will be. Default: 1000. Recommended for robustness to have till 10000 trees.
#' @param mtry n**(1/mtry) features of the n features will be sampled at each split during the training of the random forest algorithm. Default: 2 (square root).
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @param filter_genes Indicate whether 50 per cent of the genes that are the least variable in ligand-target scores should be removed in order to reduce the training of the model. Default: FALSE.
#'
#' @return A data.frame with for each ligand - data set combination, feature importance scores indicating how important the query ligand is for the prediction of the response in the particular dataset, when prediction is done via a trained regression model with all possible ligands as input. In addition to the importance score(s), the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @importFrom randomForest randomForest importance
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation_regression)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances_rf = dplyr::bind_rows(lapply(settings_ligand_pred, get_multi_ligand_rf_importances_regression,ligand_target_matrix, ntrees = 100, mtry = 2))
#' print(head(ligand_importances_rf))
#' }
#'
#' @export
#'
get_multi_ligand_rf_importances_regression = function(setting,ligand_target_matrix, ligands_position = "cols", ntrees = 1000, mtry = 2, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")
  if(ntrees <= 1)
    stop("ntrees should be higher than 1")
  if(mtry <= 1)
    stop("mtry should be higher than 1")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    ligand_target_matrix = filter_genes_ligand_target_matrix(ligand_target_matrix,ligands_position)
  }

  setting_name = setting$name
  ligands_oi = setting$from

  if (ligands_position == "cols"){
    if(sum((ligands_oi %in% colnames(ligand_target_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_target_matrix")
    prediction_matrix = ligand_target_matrix[,ligands_oi]
    target_genes = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    if(sum((ligands_oi %in% rownames(ligand_target_matrix)) == FALSE) > 0)
      stop("ligands should be in ligand_target_matrix")
    prediction_matrix = ligand_target_matrix[ligands_oi,] %>% t()
    target_genes = colnames(ligand_target_matrix)
  }

  response_vector = setting$response
  response_df = tibble(gene = names(response_vector), response = response_vector)

  prediction_df = prediction_matrix %>% data.frame() %>% as_tibble()

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  train_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()

  rf_model = randomForest::randomForest(y = train_data$obs,
                                        x = train_data[,-(which(colnames(train_data) == "obs"))],
                                        ntree = ntrees,
                                        mtry = ncol(train_data[,-(which(colnames(train_data) == "obs"))])**(1/mtry) %>% ceiling(),
                                        importance = TRUE
  )

  metrics = randomForest::importance(rf_model) %>% data.frame() %>% tibble::rownames_to_column("test_ligand") %>% as_tibble() %>% mutate(setting = setting_name)

  if (known == TRUE){
    true_ligand = setting$ligand
    metrics = metrics %>% mutate(ligand = true_ligand)
    metrics = metrics %>% select(setting, test_ligand, ligand, X.IncMSE, IncNodePurity) %>% rename(IncMSE = X.IncMSE)
    return(metrics)
  }
  metrics = metrics %>% select(setting, test_ligand, X.IncMSE, IncNodePurity) %>% rename(IncMSE = X.IncMSE)
  return(metrics)
}
#' @title Convert settings to correct settings format for TF prediction.
#'
#' @description \code{convert_settings_tf_prediction} Converts settings to correct settings format for TF activity prediction. In this prediction problem, TFs (out of a set of possibly active TFs) will be ranked based on feature importance scores. The format can be made suited for applications in which TFs need to be scored based on their possible upstream activity: 3) by calculating individual feature importane scores or 4) feature importance based on models with embedded feature importance determination. Remark that upstream regulator analysis for TFs here is experimental and was not thoroughly validated in the study accompanying this package.
#'
#' @usage
#' convert_settings_tf_prediction(settings, all_tfs, single = TRUE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the tf(s) active in the setting of interest; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest.
#' @param all_tfs A character vector of possible tfs that will be considered for the tf activity state prediction.
#' @param single TRUE if feature importance scores for tfs will be calculated by looking at ligans individually. FALSE if the goal is to calculate the feature importance scores via sophisticated classification algorithms like random forest.

#' @return A list with following elements: $name, $tf: name of active tf(s) (only if validation is TRUE), $from (tf(s) that will be tested for activity prediction), $response
#'
#' @examples
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_tf_pred = convert_settings_tf_prediction(settings, all_tfs = c("SMAD1","STAT1","RELA"), single = TRUE)
#' # show how this function can be used to predict activities of TFs
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' tf_target = construct_tf_target_matrix(weighted_networks, tfs_as_cols = TRUE, standalone_output = TRUE)
#' tf_importances = dplyr::bind_rows(lapply(settings_tf_pred,get_single_ligand_importances,tf_target,known = FALSE))
#' print(head(tf_importances))
#'
#' @export
#'
#'
convert_settings_tf_prediction = function(settings,all_tfs, single = TRUE){

  # input check
  if(!is.list(settings))
    stop("settings should be a list")
  if(!is.character(all_tfs))
    stop("all_tfs should be a character vector")
  if(!is.logical(single) | length(single) != 1)
    stop("single should be TRUE or FALSE")

  requireNamespace("dplyr")

  new_settings = list()
  if (single == TRUE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      for (k in 1:length(all_tfs)){
        test_tf = all_tfs[[k]]
        new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_single_application(setting,test_tf))
      }
    }
  } else if (single == FALSE){
    for (i in 1:length(settings)){
      setting = settings[[i]]
      new_settings[[length(new_settings) + 1]] = list(make_new_setting_ligand_prediction_multi_application(setting,all_tfs))
    }
  }
  return(new_settings %>% unlist(recursive = FALSE))
}
#' @title Converts expression settings to format in which the total number of potential ligands is reduced up to n top-predicted active ligands.
#'
#' @description \code{convert_expression_settings_evaluation} Converts expression settings to format in which the total number of potential ligands is reduced up to n top-predicted active ligands.(useful for applications when a lot of ligands are potentially active, a lot of settings need to be predicted and a multi-ligand model is trained).
#'
#' @usage
#' convert_settings_topn_ligand_prediction(setting, importances, model, n, normalization)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$diffexp: data frame or tibble containing at least 3 variables= $gene, $lfc (log fold change treated vs untreated) and $qval (fdr-corrected p-value)
#' @param n The top n number of ligands according to the ligand activity state prediction model will be considered as potential ligand for the generation of a new setting.
#' @inheritParams model_based_ligand_activity_prediction

#' @return A list with following elements: $name, $from, $response. $response will be a gene-named logical vector indicating whether the gene's transcription was influenced by the active ligand(s) in the setting of interest.
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' evaluation = evaluate_importances_ligand_prediction(ligand_importances,"median","lda")
#'
#' settings = lapply(expression_settings_validation[5:10],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = FALSE, single = TRUE)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix, known = FALSE))
#' settings = lapply(settings,convert_settings_topn_ligand_prediction, importances = ligand_importances, model = evaluation$model, n = 3, normalization = "median" )
#' }
#'
#' @export
#'
convert_settings_topn_ligand_prediction = function(setting, importances, model, n, normalization){
  # input check
  if(!is.list(setting))
    stop("setting should be a list")
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(setting$response))
    stop("setting$response should be a logical vector")

  requireNamespace("dplyr")
  setting_name = setting$name
  importances_oi = importances %>% filter(setting == setting_name)
  output = model_based_ligand_activity_prediction(importances_oi, model,normalization)
  top_ligands = output %>% top_n(n,model) %>% .$test_ligand

  new_setting = list()
  new_setting$name = setting$name
  new_setting$from = top_ligands
  new_setting$response = setting$response

  return(new_setting)
}

#' @title Evaluation of ligand activity prediction performance of single ligand importance scores: each dataset individually.
#'
#' @description \code{wrapper_evaluate_single_importances_ligand_prediction} Evaluate how well a single ligand importance score is able to predict the true activity state of a ligand. For this it is assumed, that ligand importance measures for truely active ligands will be higher than for non-active ligands. Several classification evaluation metrics for the prediction are calculated and variable importance scores can be extracted to rank the different importance measures in order of importance for ligand activity state prediction.
#'
#' @usage
#' wrapper_evaluate_single_importances_ligand_prediction(group,ligand_importances)
#'
#' @param group Name of the dataset (setting) you want to calculate ligand activity performance for.
#' @param ligand_importances A data frame containing at least folowing variables: $setting, $test_ligand, $ligand and one or more feature importance scores. $test_ligand denotes the name of a possibly active ligand, $ligand the name of the truely active ligand.
#'
#' @return A data frame containing classification evaluation measures for the ligand activity state prediction single, individual feature importance measures.
#'
#' @examples
#' \dontrun{
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#' evaluation = ligand_importances$setting %>% unique() %>% lapply(function(x){x}) %>% lapply(wrapper_evaluate_single_importances_ligand_prediction,ligand_importances) %>% bind_rows() %>% inner_join(ligand_importances %>% distinct(setting,ligand))
#' print(head(evaluation))
#' }
#' @export
#'
wrapper_evaluate_single_importances_ligand_prediction = function(group,ligand_importances){
  if (!is.data.frame(ligand_importances))
    stop("ligand_importances must be a data frame")
  if (!is.character(group))
    stop("group must be a character")
  ligand_importances %>% filter(setting %in% group) %>% evaluate_single_importances_ligand_prediction(normalization = "no") %>% mutate(setting = group)
}

