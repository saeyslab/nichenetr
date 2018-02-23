#' @title Convert settings to correct settings format for ligand prediction.
#'
#' @description \code{convert_settings_ligand_prediction} Converts settings to correct settings format for ligand activity prediction. In this prediction problem, ligands (out of a set of possibly active ligands) will be ranked based on feature importance scores. The format can be made suited for: 1) validation of ligand activity state prediction by calculating individual feature importane scores or 2) feature importance based on models with embedded feature importance determination; applications in which ligands need to be scores based on their possible upstream activity: 3) by calculating individual feature importane scores or 4) feature importance based on models with embedded feature importance determination.
#'
#' @usage
#' convert_settings_ligand_prediction(setting, all_ligands, validation = TRUE, single = TRUE)
#'
#' @param settings A list of lists. Eeach sublist contains the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: the observed target response: indicate for a gene whether it was a target or not in the setting of interest.
#' @param all_ligands A character vector of possible ligands that will be considered for the ligand activity state prediction.
#' @param validation TRUE if seetings need to be prepared for validation of ligand activity state predictions (this implies that the true active ligand of a setting is known); FALSE for application purposes when the true active ligand(s) is/are not known.
#' @param single TRUE if feature importance scores for ligands will be calculated by looking at ligans individually. FALSE if the goal is to calculate the feature importance scores via sophisticated classification algorithms like random forest.

#' @return A list with following elements: $name, $ligand: name of active ligand(s) (only if validation is TRUE), $from (ligand(s) that will be tested for activity prediction), $response
#'
#' @examples
#' settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings,combination = FALSE) %>% unlist()
#' settings_ligand_pred =  = convert_settings_ligand_prediction(settings, ligands, validation = TRUE, single = TRUE)
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
        test_ligand = ligands[[k]]
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
        test_ligand = ligands[[k]]
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
