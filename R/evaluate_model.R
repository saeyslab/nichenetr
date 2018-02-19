#' @title Extract ligands of interest from settings
#'
#' @description \code{extract_ligands_from_settings} Extract ligands of interest from (expression) settings in correct to construct the ligand-target matrix.
#'
#' @usage
#' extract_ligands_from_settings(settings)
#'
#' @param settings A list of lists for which each sub-list contains the information about (expression) datasets; with minimally the following elements: name of the setting ($name), ligands (possibly) active in the setting of interest ($from).
#'
#' @return A list containing the ligands and ligands combinations for which a ligand-target matrix should be constructed. When for a particular dataset multiple ligands are possibly active (i.e. more than ligand in .$from slot of sublist of settings), then both the combination of these multiple ligands and each of these multiple ligands individually will be select for model construction.
#'
#' @examples
#' extract_ligands_from_settings(expression_settings_validation)
#'
#' @export
#'
extract_ligands_from_settings = function(settings){

  # input check
  if (!is.list(settings))
    stop("settings must be a list")
  if(sum(sapply(settings,function(x){is.character(x$from)})) != length(settings))
    stop("settings$.$from must be a character vector containing ligands")

  ligands_oi = list()
  for (i in 1:length(settings)){
    setting = settings[[i]]
    ligand = setting$from
    if (length(ligand) == 1) {
      ligands_oi[length(ligands_oi) + 1] = ligand
    } else {# if multiple ligands added
      ligands_oi[[length(ligands_oi) + 1]] = ligand # ligands together
      for (l in ligand) {
        ligands_oi[length(ligands_oi) + 1] = l # ligands separate
      }
    }
  }
  ligands_oi = unique(ligands_oi)
  return(ligands_oi)
}
#' @title Convert expression settings to correct settings format for evaluation of target gene prediction.
#'
#' @description \code{convert_expression_settings_evaluation} Converts expression settings to correct settings format for evaluation of target gene prediction.
#'
#' @usage
#' convert_expression_settings_evaluation(setting)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$diffexp: data frame or tibble containing at least 3 variables= $gene, $lfc (log fold change treated vs untreated) and $qval (fdr-corrected p-value)

#' @return A list with following elements: $name, $from, $response
#'
#' @examples
#' settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
#'
#' @export
#'
#'
convert_expression_settings_evaluation = function(setting) {
  # input check
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.data.frame(setting$diffexp))
    stop("setting$diffexp should be data frame")
  if(is.null(setting$diffexp$lfc) | is.null(setting$diffexp$gene) | is.null(setting$diffexp$qval))
    stop("setting$diffexp should contain the variables 'lfc', 'qval' and 'diffexp'")

  requireNamespace("dplyr")

  diffexp_df = setting$diffexp %>% mutate(diffexp = (abs(lfc) >= 1) & (qval <= 0.1))
  diffexp_vector = diffexp_df$diffexp
  names(diffexp_vector) = diffexp_df$gene
  if((diffexp_vector %>% sum) == 0) {
    print(setting$name)
    warning("No differentially expressed genes, remove this expression dataset")
  }
  return(list(name = setting$name, from = setting$from, response = diffexp_vector))
}
#' @title Evaluation of target gene prediction.
#'
#' @description \code{evaluate_target_prediction} Evaluate how well the model (i.e. the inferred ligand-target probability scores) is able to predict the observed response to a ligand (e.g. the set of DE genes after treatment of cells by a ligand). It shows several classification evaluation metrics for the prediction. Different classification metrics are calculated depending on whether the input ligand-target matrix contains probability scores for targets or discrete target assignments.
#'
#' @usage
#' evaluate_target_prediction(setting,ligand_target_matrix, ligands_position = "cols")
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active ligand(s) or a FALSE.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"

#' @return A data.frame with following variables: setting, ligand and for probabilistic predictions: auroc, aupr, aupr_corrected (aupr - aupr for random prediction), sensitivity_roc (proxy measure, inferred from ROC), specificity_roc (proxy measure, inferred from ROC), mean_rank_GST_log_pval (-log10 of p-value of mean-rank gene set test), pearson (correlation coefficient), spearman (correlation coefficient); whereas for categorical predictions: accuracy, recall, specificity, precision, F1, F0.5, F2, mcc, informedness, markedness, fisher_pval_log (which is -log10 of p-value fisher exact test), fisher odds.
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @importFrom data.table data.table
#' @importFrom limma wilcoxGST
#'
#' @examples
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = lapply(expression_settings_validation[1],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(setting)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' perf1 =lapply(setting,evaluate_target_prediction,ligand_target_matrix)
#' perf2 = lapply(setting,evaluate_target_prediction,make_discrete_ligand_target_matrix(ligand_target_matrix))
#'
#' @export
#'
evaluate_target_prediction = function(setting,ligand_target_matrix, ligands_position = "cols"){
  ## still make evaluation multiple ligands possible
  # input check
  if (!is.list(setting))
    stop("setting must be a list")
  if(!is.character(setting$from) | !is.character(setting$name))
    stop("setting$from and setting$name should be character vectors")
  if(!is.logical(setting$response) | is.null(names(setting$response)))
    stop("setting$response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")

  requireNamespace("dplyr")

  if (length(setting$from) == 1){
    ligand_oi = setting$from
  } else {
    ligand_oi = paste0(setting$from,collapse = "-")
  }
  if (ligands_position == "cols"){
    if((ligand_oi %in% colnames(ligand_target_matrix)) == FALSE)
      stop("ligand should be in ligand_target_matrix")
    prediction_vector = ligand_target_matrix[,ligand_oi]
    names(prediction_vector) = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    if((ligand_oi %in% rownames(ligand_target_matrix)) == FALSE)
      stop("ligand should be in ligand_target_matrix")
    prediction_vector = ligand_target_matrix[ligand_oi,]
    names(prediction_vector) = colnames(ligand_target_matrix)
  }
  response_vector = setting$response

  performance = evaluate_target_prediction_strict(response_vector,prediction_vector,is.double(prediction_vector))
  performance = performance %>% mutate(setting = setting$name, ligand = ligand_oi) # ligand_oi should be the multi thing
}
# new_settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
# performances = evaluate_target_prediction(new_settings$rat_Tnfa,ligand_target_matrix)
