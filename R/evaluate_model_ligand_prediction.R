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
#' settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
#' ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE))
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = TRUE, single = TRUE)
#'
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
#' get_single_ligand_importances(settings,ligand_target_matrix, ligands_position = "cols", known = TRUE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
#' @param known Indicate whether the true active ligand for a particular dataset is known or not. Default: TRUE. The true ligand will be extracted from the $ligand slot of the setting.
#' @inheritParams evaluate_target_prediction
#'
#' @return A data.frame with for each ligand - data set combination, classification evaluation metrics indicating how well the query ligand predicts the response in the particular dataset. Evaluation metrics are the same as in \code{\link{evaluate_target_prediction}}. In addition to the metrics, the name of the particular setting ($setting), the name of the query ligand($test_ligand), the name of the true active ligand (if known: $ligand).
#'
#' @examples
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = get_single_ligand_importances(settings_ligand_pred,ligand_target_matrix)
#'
#' @export
#'
get_single_ligand_importances = function(settings,ligand_target_matrix, ligands_position = "cols", known = TRUE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")

  importances = bind_rows(lapply(settings,wrapper_evaluate_target_prediction_ligand_prediction, ligand_target_matrix, ligands_position,known))

  return(importances)
}
#' @title Get ligand importances from a multi-ligand classfication model.
#'
#' @description \code{get_multi_ligand_importances} A classificiation algorithm chosen by the user is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features) in order to predict the observed response in a particular dataset. Variable importance scores that indicate for each ligand the importance for response prediction, are extracted. It can be assumed that ligands with higher variable importance scores are more likely to be a true active ligand.
#'
#' @usage
#' get_multi_ligand_importances(settings,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) of which the predictve performance need to be assessed; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest. $ligand: NULL or the name of the ligand(s) that are known to be active in the setting of interest.
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
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = FALSE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances_glm = get_multi_ligand_importances(settings_ligand_pred,ligand_target_matrix, algorithm = "glm")
#'
#' @export
#'
get_multi_ligand_importances = function(settings,ligand_target_matrix, ligands_position = "cols", algorithm, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE, known = TRUE, filter_genes = FALSE){

  if(!is.logical(known) | length(known) > 1)
    stop("known should be a logical vector: TRUE or FALSE")
  if(!is.logical(filter_genes) | length(filter_genes) > 1)
    stop("filter_genes should be a logical vector: TRUE or FALSE")

  requireNamespace("dplyr")

  if (filter_genes == TRUE){
    ligand_target_matrix = filter_genes_ligand_target_matrix(ligand_target_matrix,ligands_position)
  }

  importances = bind_rows(lapply(settings,wrapper_evaluate_target_prediction_multi_ligand_prediction, ligand_target_matrix, ligands_position, algorithm, cv, cv_number, cv_repeats, parallel, n_cores, ignore_errors, continuous, known))

  return(importances)

}






