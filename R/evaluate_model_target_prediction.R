#' @title Extract ligands of interest from settings
#'
#' @description \code{extract_ligands_from_settings} Extract ligands of interest from (expression) settings in correct to construct the ligand-target matrix.
#'
#' @usage
#' extract_ligands_from_settings(settings,combination = TRUE)
#'
#' @param settings A list of lists for which each sub-list contains the information about (expression) datasets; with minimally the following elements: name of the setting ($name), ligands (possibly) active in the setting of interest ($from).
#' @param combination Indicate whether in case multiple ligands are possibly active ligand combinations should be extracted or only individual ligands. Default: TRUE.
#'
#' @return A list containing the ligands and ligands combinations for which a ligand-target matrix should be constructed. When for a particular dataset multiple ligands are possibly active (i.e. more than ligand in .$from slot of sublist of settings), then both the combination of these multiple ligands and each of these multiple ligands individually will be select for model construction.
#'
#' @examples
#' ligands = extract_ligands_from_settings(expression_settings_validation)
#'
#' @export
#'
extract_ligands_from_settings = function(settings, combination = TRUE){

  # input check
  if (!is.list(settings))
    stop("settings must be a list")
  if(sum(sapply(settings,function(x){is.character(x$from)})) != length(settings))
    stop("settings$.$from must be a character vector containing ligands")

  ligands_oi = list()
  if (combination == TRUE){
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
  } else {
    for (i in 1:length(settings)){
      setting = settings[[i]]
      ligand = setting$from
      if (length(ligand) == 1) {
        ligands_oi[length(ligands_oi) + 1] = ligand
      } else {# if multiple ligands added
          for (l in ligand) {
          ligands_oi[length(ligands_oi) + 1] = l # ligands separate
        }
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
  diffexp_vector = diffexp_vector[unique(names(diffexp_vector))]
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
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (or discrete target assignments).
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
#' perf1 = lapply(setting,evaluate_target_prediction,ligand_target_matrix)
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
  output = bind_cols(tibble(setting = setting$name, ligand = ligand_oi), performance)

  return(output)
}
#' @title Evaluation of target gene prediction for multiple ligands.
#'
#' @description \code{evaluate_multi_ligand_target_prediction} Evaluate how well a trained model is able to predict the observed response to a combination of ligands (e.g. the set of DE genes after treatment of cells by multiple ligands). A classificiation algorithm chosen by the user is trained to construct one model based on the target gene predictions of all ligands of interest (ligands are considered as features). Several classification evaluation metrics for the prediction are calculated depending on whether the input ligand-target matrix contains probability scores for targets or discrete target assignments. In addition, variable importance scores can be extracted to rank the possible active ligands in order of importance for response prediction.
#'
#' @usage
#' evaluate_multi_ligand_target_prediction(setting,ligand_target_matrix, ligands_position = "cols", algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4,ignore_errors = FALSE,continuous = TRUE)
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active ligand(s) or a FALSE.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (recommended) or discrete target assignments (not-recommended).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#' @param algorithm The name of the classification algorithm to be applied. Should be supported by the caret package. Examples of algorithms we recommend: with embedded feature selection: "rf","glm","fda","glmnet","sdwd","gam","glmboost"; without: "lda","naive_bayes","pls"(because bug in current version of pls package), "pcaNNet". Please notice that not all these algorithms work when the features (i.e. ligand vectors) are categorical (i.e. discrete class assignments).
#' @param var_imps Indicate whether in addition to classification evaluation performances, variable importances should be calculated. Default: TRUE.
#' @param cv Indicate whether model training and hyperparameter optimization should be done via cross-validation. Default: TRUE. FALSE might be useful for applications only requiring variable importance, or when final model is not expected to be extremely overfit.
#' @param cv_number The number of folds for the cross-validation scheme: Default: 4; only relevant when cv == TRUE.
#' @param cv_repeats The number of repeats during cross-validation. Default: 2; only relevant when cv == TRUE.
#' @param parallel Indiciate whether the model training will occur parallelized. Default: FALSE. TRUE only possible for non-windows OS.
#' @param n_cores The number of cores used for parallelized model training via cross-validation. Default: 4. Only relevant on non-windows OS.
#' @param ignore_errors Indiciate whether errors during model training by caret should be ignored such that another model training try will be initiated until model is trained without raising errors. Default: FALSE.
#' @param continuous Indicate whether during training of the model, model training and evaluation should be done on class probabilities or discrete class labels. For huge class imbalance, we recommend setting this value to TRUE. Default: TRUE.
#'
#' @return A list with the following elements. $performances: data frame containing classification evaluation measure for classification on the test folds during training via cross-validation; $performances_training: data frame containing classification evaluation measures for classification of the final model (discrete class assignments) on the complete data set (performance can be severly optimistic due to overfitting!); $performance_training_continuous: data frame containing classification evaluation measures for classification of the final model (class probability scores) on the complete data set (performance can be severly optimistic due to overfitting!) $var_imps: data frame containing the variable importances of the different ligands (embbed importance score for some classification algorithms, otherwise just the auroc); $prediction_response_df: data frame containing for each gene the ligand-target predictions of the individual ligands, the complete model and the response as well; $setting: name of the specific setting that needed to be evaluated; $ligands: ligands of interest.
#'
#' @importFrom ROCR prediction performance
#' @importFrom caTools trapz
#' @importFrom limma wilcoxGST
#' @import caret
#' @importFrom purrr safely
#'
#' @examples
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = convert_expression_settings_evaluation(expression_settings_validation$TGFB_IL6_timeseries) %>% list()
#' ligands = extract_ligands_from_settings(setting)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' output = lapply(setting,evaluate_multi_ligand_target_prediction,ligand_target_matrix,ligands_position = "cols",algorithm = "glm")
#' output = lapply(setting,evaluate_multi_ligand_target_prediction,make_discrete_ligand_target_matrix(ligand_target_matrix),ligands_position = "cols",algorithm = "glm" )
#'
#' @export
#'
evaluate_multi_ligand_target_prediction = function(setting,ligand_target_matrix, ligands_position = "cols", algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4, ignore_errors = FALSE, continuous = TRUE){
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
  if(!is.logical(continuous) | length(continuous) > 1)
    stop("continuous should be a logical vector: TRUE or FALSE")
  requireNamespace("dplyr")

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

  prediction_df = prediction_matrix %>% data.frame() %>% tbl_df()

  if(is.double(prediction_matrix) == FALSE){
    convert_categorical_factor = function(x){
      x = x %>% make.names() %>% as.factor()
    }
    prediction_df = prediction_df %>% mutate_all(funs(convert_categorical_factor))
  }

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  train_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()

  output = wrapper_caret_classification(train_data,algorithm,continuous = continuous,var_imps,cv,cv_number,cv_repeats,parallel,n_cores,prediction_response_df = combined,ignore_errors)
  output$setting = setting$name
  output$ligands = ligands_oi
  return(output)
}
#' @title Evaluation of target gene prediction.
#'
#' @description \code{evaluate_target_prediction_interprete} Evaluate how well the model (i.e. the inferred ligand-target probability scores) is able to predict the observed response to a ligand (e.g. the set of DE genes after treatment of cells by a ligand). It shows several classification evaluation metrics for the prediction. Different classification metrics are calculated depending on whether the input ligand-target matrix contains probability scores for targets or discrete target assignments.
#'
#' @usage
#' evaluate_target_prediction_interprete(setting,ligand_target_matrix, ligands_position = "cols")
#'
#' @param setting A list containing the following elements: .$name: name of the setting; .$from: name(s) of the ligand(s) active in the setting of interest; .$response: named logical vector indicating whether a target is a TRUE target of the possibly active ligand(s) or a FALSE.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (or discrete target assignments).
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"

#' @return A list with the elements $performances and $prediction_response_df. $performance is a data.frame with following variables: setting, ligand and for probabilistic predictions: auroc, aupr, aupr_corrected (aupr - aupr for random prediction), sensitivity_roc (proxy measure, inferred from ROC), specificity_roc (proxy measure, inferred from ROC), mean_rank_GST_log_pval (-log10 of p-value of mean-rank gene set test), pearson (correlation coefficient), spearman (correlation coefficient); whereas for categorical predictions: accuracy, recall, specificity, precision, F1, F0.5, F2, mcc, informedness, markedness, fisher_pval_log (which is -log10 of p-value fisher exact test), fisher odds. $prediction_response_df shows for each gene, the model prediction (ligand-target probability score or discrete target assignment TRUE/FALSE) and whether the gene the gene is a target or not according to the observed response.
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
#' perf1 = lapply(setting,evaluate_target_prediction_interprete,ligand_target_matrix)
#'
#' @export
#'
evaluate_target_prediction_interprete = function(setting,ligand_target_matrix, ligands_position = "cols"){
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

  output = evaluate_target_prediction_strict(response_vector,prediction_vector,is.double(prediction_vector), prediction_response_df = TRUE)
  output$setting = setting$name
  output$ligand = ligand_oi
  colnames(output$prediction_response_df) = c("gene","response",ligand_oi)

  return(output)
}
#' @title Convert gene list to correct settings format for evaluation of target gene prediction.
#'
#' @description \code{convert_gene_list_settings_evaluation} Converts a gene list to correct settings format for evaluation of target gene prediction.
#'
#' @usage
#' convert_gene_list_settings_evaluation(gene_list, name, ligands_oi, background)
#'
#' @param gene_list A character vector of target gene names
#' @param name The name that will be given to the setting
#' @param ligands_oi The possibly active ligands
#' @param background A character vector of names of all genes (contains thus genes belonging to the gene_list of interest as well).

#' @return A list with following elements: $name, $from, $response
#'
#' @examples
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' all_genes = unique(c(weighted_networks$gr$from,weighted_networks$gr$to,weighted_networks$lr_sig$from, weighted_networks$lr_sig$to))
#' gene_list = c("ID1","ID2","ID3")
#' setting = list(convert_gene_list_settings_evaluation(gene_list = c("ID1","ID2","ID3"), name = "test",ligands_oi = "TGFB1", background = all_genes))
#'
#' @export
#'
#'
convert_gene_list_settings_evaluation = function(gene_list, name, ligands_oi, background) {
  # input check
  if(!is.character(gene_list))
    stop("gene_list should be character vector")
  if(!is.character(name) | length(name) > 1)
    stop("name should be character vector of length 1")
  if(!is.character(ligands_oi))
    stop("ligands_oi should be character vector")
  if(!is.character(background))
    stop("background should be character vector")

  requireNamespace("dplyr")

  excluded_genes = gene_list[(gene_list %in% background) == FALSE]
  if(length(excluded_genes) > 0)
    warning(paste("some genes are excluded because not present in the background: ",excluded_genes))
  gene_list_vector = (background %in% gene_list)
  names(gene_list_vector) = background

  return(list(name = name, from = ligands_oi, response = gene_list_vector))
}



