#' @title Get the number of times of gene is mentioned in the pubmed literature
#'
#' @description \code{get_ncitations_genes}: Get the number of times of gene is mentioned in the pubmed literature
#'
#' @usage
#' get_ncitations_genes(file = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz")
#'
#' @param file A file containing a data frame denoting the pubmed ids in which a particular gene entrez of a particulat species is mentioned (variables: taxid, entrez, pubmedid). Default = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"
#'
#' @return A data.frame with following variables: entrez, ncitations, symbol,entrez_mouse, symbol_mouse
#'
#' @import readr
#'
#' @examples
#' \dontrun{
#' ncitations = get_ncitations_genes(file = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz")
#' }
#'
#' @export
#'
get_ncitations_genes = function(file = "ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz"){

  if (!is.character(file))
    stop("file must be a character vector")

  requireNamespace("dplyr")
  requireNamespace("readr")
  read_tsv(file, col_names = c("taxid", "entrez", "pubmedid"), skip = 1, col_types = cols(
    taxid = col_integer(),
    entrez = col_character(),
    pubmedid = col_integer())) %>% filter(taxid==9606) %>% group_by(entrez) %>% summarize(ncitations=n()) %>% left_join(geneinfo_human, "entrez") %>% arrange(-ncitations)
}
#' @title Merge target gene prediction performances with popularity measures of ligands
#'
#' @description \code{add_ligand_popularity_measures_to_perfs}: Get a data.frame in which the performance measures for target gene prediction of a ligand are merged with the number of times the ligand is mentioned in the Pubmed literature. Serves to investigate popularity bias (i.e. it can be expected that frequenly cited ligands will have better predictive performance because they are better studied).
#'
#' @usage
#' add_ligand_popularity_measures_to_perfs(performances, ncitations)
#'
#' @param performances A data.frame in which the performance measures for target gene predictions of ligands are denoted
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'. Default: ncitations (variable contained in this package). See function \code{get_ncitations_genes} for a function that makes this data frame from current Pubmed information.
#'
#' @return A data.frame in which the performance measures for target gene prediction of a ligand are merged with the number of times the ligand is mentioned in the Pubmed literature.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' performances = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix))
#' # ncitations = get_ncitations_genes()
#' performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances,ncitations)
#' }
#'
#' @export
#'
add_ligand_popularity_measures_to_perfs = function(performances,ncitations){
  if (!is.data.frame(performances))
    stop("performances must be a data frame")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")
  if(!is.character(performances$ligand))
    stop("performances$ligand should be a character vector")

  requireNamespace("dplyr")
  performances_ligand_pop = performances %>% inner_join(ncitations %>% rename(ligand = symbol) %>% select(ligand,ncitations), by = "ligand")
  return(performances_ligand_pop)
}
#' @title Regression analysis between ligand popularity and target gene predictive performance
#'
#' @description \code{get_slope_ligand_popularity}: Performs regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of the ligand.
#'
#' @usage
#' get_slope_ligand_popularity(metric, performances)
#'
#' @param metric The name of the performance metric of which the trend with the popularity of the ligand should be calculated.
#' @param performances A data.frame in which the performance measures for target gene predictions of ligands are denoted together with the popularity of the ligand. (should contain at least following variables: ligand, ncitations and the metric of interest)
#'
#' @return A data.frame in which the regression coefficient estimate, p-value and corresponding R-squared value are shown for the regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of the ligand.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' performances = bind_rows(lapply(settings,evaluate_target_prediction,ligand_target_matrix))
#' # ncitations = get_ncitations_genes()
#' performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances,ncitations)
#' slopes_auroc = get_slope_ligand_popularity("auroc",performances_ligand_popularity)
#' slopes_df = performances_ligand_popularity %>% select(-setting,-ligand,-ncitations) %>% colnames() %>% lapply(.,get_slope_ligand_popularity,performances_ligand_popularity) %>% bind_rows()
#' }
#'
#'
#' @export
#'
get_slope_ligand_popularity = function(metric,performances){
  if (!is.data.frame(performances))
    stop("performances must be a data frame")
  if (!is.character(metric))
    stop("metric must be a character vector")
  if((metric %in% colnames(performances)) == FALSE)
    stop("metric must be in colnames of performances")
  if(!is.numeric(performances$ncitations))
    stop("performances$ncitations should be a numeric vector")

  requireNamespace("dplyr")
  performances = performances %>% select(paste(metric),ncitations)
  colnames(performances) = c("metric","ncitations")



  if (nrow(performances) == 0) {
    output = tibble(metric = metric, ligand_slope = NA, ligand_slope_pval =  NA, ligand_slope_rsquared = NA)
    colnames(output) = variable
    return(output)
    }
  ligand_pop_regression = lm(metric ~ log(ncitations),performances)
  ligand_slope = summary(ligand_pop_regression) %>% .$coefficients %>% .[2,1]
  ligand_slope_pval = summary(ligand_pop_regression) %>% .$coefficients %>% .[2,4]
  ligand_slope_rsquared = summary(ligand_pop_regression) %>% .$r.squared
  output = tibble(metric = metric, ligand_slope = ligand_slope, ligand_slope_pval =  ligand_slope_pval, ligand_slope_rsquared = ligand_slope_rsquared)
  return(output)
}
#' @title Evaluate target gene predictions for different bins/groups of targets genes
#'
#' @description \code{evaluate_target_prediction_per_bin}: Evaluate target gene predictions for different bins/groups of targets genes. Bins are constructed such that genes that are similarly frequently cited are grouped together and the different bins have similar size.
#'
#' @usage
#' evaluate_target_prediction_per_bin(nbins,settings,ligand_target_matrix,ncitations,ligands_position = "cols")
#'
#' @param nbins The number of bins the target genes should be divided in based on popularity.
#' @param settings  list of lists for which each sub-list contains the information about (expression) datasets; with minimally the following elements: name of the setting ($name), ligands (possibly) active in the setting of interest ($from).
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores (or discrete target assignments).
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'. Default: ncitations (variable contained in this package). See function \code{get_ncitations_genes} for a function that makes this data frame from current Pubmed information.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A data.frame containing several classification evaluation metrics for target gene prediction. Predictions were evaluated for n different bins of target genes. The specific bin is indicated in the variable target_bin_id. target_bin_id = 1: target genes that are least mentioned in the Pubmed literature.
#'
#' @importFrom Hmisc cut2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' # ncitations = get_ncitations_genes()
#' performances_target_bins_popularity = evaluate_target_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
#' }
#'
#' @export
#'
evaluate_target_prediction_per_bin = function(nbins,settings,ligand_target_matrix,ncitations,ligands_position = "cols"){
  if (!is.numeric(nbins))
    stop("nbins must be a numeric vector")
  if (!is.list(settings))
    stop("settings must be a list")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")

  if (ligands_position == "cols"){
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    targets = colnames(ligand_target_matrix)
  }
  ncitations = ncitations %>% filter(symbol %in% targets) %>% mutate(bin = Hmisc::cut2(log(ncitations), m = 1000, g = nbins) %>% as.numeric())
  performances = bind_rows(lapply(seq(nbins),evaluate_target_prediction_bins_direct,settings,ligand_target_matrix,ligands_position,ncitations))
}
#' @title Regression analysis between target gene popularity and target gene predictive performance
#'
#' @description \code{get_slope_target_gene_popularity}: Performs regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of target genes.
#'
#' @usage
#' get_slope_target_gene_popularity(metric,performances,method = "individual")
#'
#' @param metric The name of the performance metric of which the trend with the popularity of the target genes should be calculated.
#' @param performances A data.frame in which the performance measures for target gene predictions of ligands are denoted together with the popularity bin of the target genes for which predictions were evaluated (should contain at least following variables: target_bin_id and the metric of interest)
#' @param method 'All': calculate slope by considering all datasets in settings. 'Individual': calculate slope for every dataset in settings separately to investigate dataset-specific popularity bias. Default: 'individual'.
#'
#' @return A data.frame in which the regression coefficient estimate, p-value and corresponding R-squared value are shown for the regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of the target genes.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' # ncitations = get_ncitations_genes()
#' performances_target_bins_popularity = evaluate_target_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
#' slopes_auroc = get_slope_target_gene_popularity("auroc",performances_target_bins_popularity)
#' slopes_df = performances_target_bins_popularity %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity,method = "individual") %>% bind_rows()
#' slopes_df2 = performances_target_bins_popularity %>% select(-setting,-ligand,-target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity,performances_target_bins_popularity,method = "all") %>% bind_rows()
#' }
#'
#' @export
#'
get_slope_target_gene_popularity = function(metric,performances,method = "individual"){
  if (!is.data.frame(performances))
    stop("performances must be a data frame")
  if (!is.character(metric))
    stop("metric must be a character vector")
  if((metric %in% colnames(performances)) == FALSE)
    stop("metric must be in colnames of performances")
  if(!is.numeric(performances$target_bin_id))
    stop("performances$target_bin_id should be a numeric vector")
  if(method != "individual" & method != "all")
    stop("method should be 'indvidual' or 'all'")

  requireNamespace("dplyr")
  # performance per bin!
  performances = performances %>% select(setting,target_bin_id, paste(metric))
  colnames(performances) = c("setting","bin_id", "metric")



  if (method == "all"){
    target_pop_regression = lm(metric ~ bin_id,performances)
    target_slope =  summary(target_pop_regression) %>% .$coefficients %>% .[2,1]
    target_slope_pval = summary(target_pop_regression) %>% .$coefficients %>% .[2,4]
    target_slope_rsquared = summary(target_pop_regression) %>% .$r.squared
    output = tibble(metric = metric, target_slope = target_slope, target_slope_pval = target_slope_pval, target_slope_rsquared = target_slope_rsquared)
    return(output)
  } else if (method == "individual"){
    settings_names = performances$setting %>% unique()
    slope_df = bind_rows(lapply(settings_names,function(x){
      performances = performances %>% filter(setting == x)
      target_pop_regression = lm(metric ~ bin_id,performances)
      target_slope =  summary(target_pop_regression) %>% .$coefficients %>% .[2,1]
      target_slope_pval = summary(target_pop_regression) %>% .$coefficients %>% .[2,4]
      target_slope_rsquared = summary(target_pop_regression) %>% .$r.squared

      output = tibble(setting = x,  metric = metric, target_slope = target_slope, target_slope_pval = target_slope_pval, target_slope_rsquared = target_slope_rsquared)

      return(output)
    }))
  }
}
#' @title Calculate ligand activity performance without considering evaluation datasets belonging to the top i most frequently cited ligands
#'
#' @description \code{ligand_activity_performance_top_i_removed}: calculates ligand activity performance (given ligand importances measures) without considering evaluation datasets belonging to the top i most frequently cited ligands.
#'
#' @usage
#' ligand_activity_performance_top_i_removed(i, importances, ncitations)
#'
#' @param i Indicate the number of most popular ligands for which the corresponding datasets will be removed for performance calculation.
#' @param importances Data frame containing the ligand importance measures for every ligand-dataset combination. See \code{get_single_ligand_importances} and \code{get_multi_ligand_importances}.
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'. Default: ncitations (variable contained in this package). See function \code{get_ncitations_genes} for a function that makes this data frame from current Pubmed information.
#'
#' @return A data.frame in which you can find several evaluation metrics of the ligand activity prediction performance and an indication of what percentage of most popular ligands has not been considered ($popularity_index).
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#'
#' ligand_activity_popularity_bias = lapply(0:3,ligand_activity_performance_top_i_removed, ligand_importances, ncitations) %>% bind_rows()
#' # example plot
#' # ligand_activity_popularity_bias %>% ggplot(aes(popularity_index,aupr)) + geom_smooth(method = "lm") + geom_point()
#'
#' }
#'
#' @export
#'
ligand_activity_performance_top_i_removed = function(i, importances, ncitations){
  if(i < 0)
    stop("i should be 0 or higher")
  if (!is.data.frame(importances))
    stop("importances must be a data frame")
  if(!is.character(importances$setting) | !is.character(importances$test_ligand) | !is.character(importances$ligand))
    stop("importances$setting, importances$test_ligand and importances$ligand should be character vectors")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")

  all_ligands = importances$test_ligand %>% unique()
  total_ligands = length(all_ligands)
  if (0.8 * total_ligands < i)
    warning("Be aware that most validation datasets are removed now and that ligand activty performances might not be robust anymore because only a few datasets left")
  if (i >= total_ligands)
    stop("The number of ligands removed cannot exceed the total number of different ligands")

  requireNamespace("dplyr")

  if (i > 0){
    ncitations_filtered = ncitations %>% distinct(symbol, ncitations) %>% filter(symbol %in% all_ligands)
    top_i_ligands = ncitations_filtered %>% top_n(i, ncitations) %>% arrange(-ncitations) %>% .$symbol
    top_ligands_active = are_ligands_oi_active(importances, top_i_ligands)
    importances = importances %>% mutate(active = top_ligands_active)
    importances_filtered = importances %>% filter(active == FALSE) %>% select(-active)
  } else {
    importances_filtered = importances
  }

  # evaluation = suppressWarnings(evaluate_importances_ligand_prediction(importances_filtered,"median","lda",cv_number = 3, cv_repeats = 20))
  # performances = evaluation$performances %>% select(-Resample) %>% mutate_all(mean) %>% distinct() %>% mutate(popularity_index = (i)/(total_ligands))
  # print(importances_filtered)

  performances_ligand_prediction_single = evaluate_single_importances_ligand_prediction(importances_filtered, "median")
  performances = performances_ligand_prediction_single %>% filter(importance_measure == "pearson" | importance_measure == "ipa_pval") %>% mutate(popularity_index = (i)/(total_ligands))
  # performances = performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% mutate(popularity_index = (i)/(total_ligands))
  return(performances)
}
#' @title Regression analysis between popularity of left-out ligands for ligand activity prediction performance
#'
#' @description \code{get_ligand_slope_ligand_prediction_popularity}: Performs regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of left-out ligands for ligand activity prediction performance.
#'
#' @usage
#' get_ligand_slope_ligand_prediction_popularity(metric,performances)
#'
#' @param metric The name of the performance metric of which the trend with the popularity of the target genes should be calculated.
#' @param performances A data.frame in which the performance measures for target gene predictions of ligands are denoted together with the "popularity index" indicating which percentage of most popular ligands were left out.)
#'
#' @return A data.frame in which the regression coefficient estimate, p-value and corresponding R-squared value are shown for the regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of the left out ligands.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' settings = lapply(expression_settings_validation[1:5],convert_expression_settings_evaluation)
#' settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands = unlist(extract_ligands_from_settings(settings,combination = FALSE)), validation = TRUE, single = TRUE)
#'
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' ligands = extract_ligands_from_settings(settings_ligand_pred,combination = FALSE)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_importances = dplyr::bind_rows(lapply(settings_ligand_pred,get_single_ligand_importances,ligand_target_matrix))
#'
#' ligand_activity_popularity_bias = lapply(0:3,ligand_activity_performance_top_i_removed, ligand_importances, ncitations) %>% bind_rows()
#' slopes_auroc = get_ligand_slope_ligand_prediction_popularity("auroc",ligand_activity_popularity_bias)
#' slopes_df = ligand_activity_popularity_bias %>% select(-importance_measure,-popularity_index) %>% colnames() %>% lapply(.,get_ligand_slope_ligand_prediction_popularity,ligand_activity_popularity_bias) %>% bind_rows()
#' }
#'
#' @export
#'
get_ligand_slope_ligand_prediction_popularity = function(metric,performances){
  if (!is.data.frame(performances))
    stop("performances must be a data frame")
  if (!is.character(metric))
    stop("metric must be a character vector")
  if((metric %in% colnames(performances)) == FALSE)
    stop("metric must be in colnames of performances")
  if(!is.numeric(performances$popularity_index))
    stop("performances$popularity_index should be a numeric vector")


  requireNamespace("dplyr")

  performances = performances %>% select(popularity_index, paste(metric))
  colnames(performances) = c("popularity_index", "metric")
#


  ligand_prediction_pop_regression = lm(metric ~ popularity_index,performances)
  ligand_prediction_slope =  summary(ligand_prediction_pop_regression) %>% .$coefficients %>% .[2,1]
  ligand_prediction_slope_pval = summary(ligand_prediction_pop_regression) %>% .$coefficients %>% .[2,4]
  ligand_prediction_slope_rsquared = summary(ligand_prediction_pop_regression) %>% .$r.squared
  output = tibble(metric = metric, ligand_prediction_slope = ligand_prediction_slope, ligand_prediction_slope_pval = ligand_prediction_slope_pval, ligand_prediction_slope_rsquared = ligand_prediction_slope_rsquared)
  return(output) # remark: the more negative the slope, the most extreme the decrease in performance if popular ligands are left out, so the higher the popularity bias
}
#' @title Evaluate ligand activity predictions for different bins/groups of targets genes
#'
#' @description \code{evaluate_ligand_prediction_per_bin}: Evaluate ligand activity predictions for different bins/groups of targets genes. Bins are constructed such that genes that are similarly frequently cited are grouped together and the different bins have similar size.
#'
#' @usage
#' evaluate_ligand_prediction_per_bin(nbins,settings,ligand_target_matrix,ncitations,ligands_position = "cols",...)
#'
#' @inheritParams  evaluate_target_prediction_per_bin
#' @param ... Additional arguments to \code{make_discrete_ligand_target_matrix}.
#'
#' @return A data.frame containing several classification evaluation metrics for ligand activity prediction. Predictions were evaluated for n different bins of target genes. The specific bin is indicated in the variable target_bin_id. target_bin_id = 1: target genes that are least mentioned in the Pubmed literature.
#'
#' @importFrom Hmisc cut2
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' # ncitations = get_ncitations_genes()
#' ligand_prediction_performances_target_bins_popularity = evaluate_ligand_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
#' }
#'
#' @export
#'
evaluate_ligand_prediction_per_bin = function(nbins,settings,ligand_target_matrix,ncitations,ligands_position = "cols",...){
  if (!is.numeric(nbins))
    stop("nbins must be a numeric vector")
  if (!is.list(settings))
    stop("settings must be a list")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (!is.data.frame(ncitations))
    stop("ncitations must be a data frame")
  if(!is.character(ncitations$symbol) | !is.numeric(ncitations$ncitations))
    stop("ncitations$symbol should be a character vector and ncitations$ncitations a numeric vector")

  if (ligands_position == "cols"){
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows") {
    targets = colnames(ligand_target_matrix)
  }
  ncitations = ncitations %>% filter(symbol %in% targets) %>% mutate(bin = Hmisc::cut2(log(ncitations), m = 1000, g = nbins) %>% as.numeric())
  performances = bind_rows(lapply(seq(nbins),evaluate_ligand_prediction_bins_direct,settings,ligand_target_matrix,ligands_position,ncitations,...))
}
#' @title Regression analysis between target gene popularity and ligand activity predictive performance
#'
#' @description \code{get_slope_target_gene_popularity_ligand_prediction}: Performs regression analysis to investigate the trend between a particular classficiation evaluation metric (ligand activity prediction) and the popularity of target genes.
#'
#' @usage
#' get_slope_target_gene_popularity_ligand_prediction(metric,performances)
#'
#' @param metric The name of the performance metric of which the trend with the popularity of the target genes should be calculated.
#' @param performances A data.frame in which the performance measures for ligand activity predictions of ligands are denoted together with the popularity bin of the target genes for which predictions were evaluated (should contain at least following variables: target_bin_id and the metric of interest)
#'
#' @return A data.frame in which the regression coefficient estimate, p-value and corresponding R-squared value are shown for the regression analysis to investigate the trend between a particular classficiation evaluation metric and the popularity of the target genes.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' # ncitations = get_ncitations_genes()
#' performances_target_bins_popularity = evaluate_ligand_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
#' slopes_auroc = get_slope_target_gene_popularity_ligand_prediction("auroc",performances_target_bins_popularity)
#' slopes_df = performances_target_bins_popularity %>% select(-importance_measure, -target_bin_id) %>% colnames() %>% lapply(.,get_slope_target_gene_popularity_ligand_prediction,performances_target_bins_popularity) %>% bind_rows()
#' }
#'
#' @export
#'
get_slope_target_gene_popularity_ligand_prediction = function(metric,performances){
  if (!is.data.frame(performances))
    stop("performances must be a data frame")
  if (!is.character(metric))
    stop("metric must be a character vector")
  if((metric %in% colnames(performances)) == FALSE)
    stop("metric must be in colnames of performances")
  if(!is.numeric(performances$target_bin_id))
    stop("performances$target_bin_id should be a numeric vector")


  requireNamespace("dplyr")
  # performance per bin!
  performances = performances %>% select(target_bin_id, paste(metric))
  colnames(performances) = c("bin_id", "metric")


  target_pop_regression = lm(metric ~ bin_id,performances)
  target_slope =  summary(target_pop_regression) %>% .$coefficients %>% .[2,1]
  target_slope_pval = summary(target_pop_regression) %>% .$coefficients %>% .[2,4]
  target_slope_rsquared = summary(target_pop_regression) %>% .$r.squared
  output = tibble(metric = metric, target_slope = target_slope, target_slope_pval = target_slope_pval, target_slope_rsquared = target_slope_rsquared)
  return(output)

}

