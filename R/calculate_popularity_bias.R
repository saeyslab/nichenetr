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
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'.
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
#' ncitations = get_ncitations_genes()
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
  performances = performances %>% tidyr::drop_na()
  performances_ligand_pop = performances %>% left_join(ncitations %>% rename(ligand = symbol) %>% select(ligand,ncitations), by = "ligand")
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
#' ncitations = get_ncitations_genes()
#' performances_ligand_popularity = add_ligand_popularity_measures_to_perfs(performances,ncitations)
#' slopes_df = performances_ligand_popularity %>% select(-setting,-ligand,-ncitations) %>% colnames() %>% lapply(.,get_slope_ligand_popularity,performances_ligand_popularity) %>% bind_rows()
#' }
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
#' @param ncitations A data frame denoting the number of times a gene is mentioned in the Pubmed literature. Should at least contain following variables: 'symbol' and 'ncitations'.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A data.frame containing several classification evaluation metrics for target gene prediction. Predictions were evaluated for n different bins of target genes. The specific bin is indicated in the variable target_bin_id. target_bin_id = 1: target genes that are least mentioned in the Pubmed literature.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' settings = lapply(expression_settings_validation[1:10],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(settings)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ncitations = get_ncitations_genes()
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
#' ncitations = get_ncitations_genes()
#' performances_target_bins_popularity = evaluate_target_prediction_per_bin(5,settings,ligand_target_matrix,ncitations)
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
