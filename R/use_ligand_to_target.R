#' @title Get the predicted top n percentage target genes of a ligand of interest
#'
#' @description \code{extract_top_fraction_targets} Get the predicted top n percentage target genes of a ligand of interest.
#'
#' @usage
#' extract_top_fraction_targets(ligand_oi,top_fraction,ligand_target_matrix,ligands_position = "cols")
#'
#' @param ligand_oi The ligand of interest of which top target genes should be returned
#' @param top_fraction A number between 0 and 1 indicating which top fraction of target genes should be returned.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A named numeric vector of ligand-target gene probability scores of the top target genes.
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' targets = extract_top_fraction_targets("BMP2",0.01,ligand_target_matrix)
#' }
#' @export
#'
extract_top_fraction_targets = function(ligand_oi,top_fraction,ligand_target_matrix,ligands_position = "cols"){
  # input check
  if (!is.character(ligand_oi) | length(ligand_oi) != 1)
    stop("ligand_oi must be a character vector of length 1")
  if (!is.numeric(top_fraction) | length(top_fraction) != 1 | top_fraction > 1 | top_fraction < 0)
    stop("top_fraction_oi must be a numeric vector of length 1 of which the value is between 0 and 1")
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (ligands_position == "cols"){
    ligands = colnames(ligand_target_matrix)
  } else if (ligands_position == "rows"){
    ligands = rownames(ligand_target_matrix)
  }
  if ((ligand_oi %in% ligands) == FALSE)
    stop("ligand_oi must be in ligand_target_matrix")

  requireNamespace("dplyr")
  # ligand should be in matrix
  if (ligands_position == "cols"){
    ligand_oi_vector = ligand_target_matrix[,ligand_oi]
  } else if (ligands_position == "rows"){
    ligand_oi_vector = ligand_target_matrix[ligand_oi,]
  }
  ligand_oi_vector[ligand_oi_vector >= quantile(ligand_oi_vector,1-top_fraction)] %>% sort(decreasing = T)
}
#' @title Get the predicted top n target genes of a ligand of interest
#'
#' @description \code{extract_top_n_targets} Get the predicted top n  target genes of a ligand of interest.
#'
#' @usage
#' extract_top_n_targets(ligand_oi,top_n,ligand_target_matrix,ligands_position = "cols")
#'
#' @param ligand_oi The ligand of interest of which top target genes should be returned
#' @param top_n A number between 0 and the total nr of target genes indicating which top n of target genes should be returned.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A named numeric vector of ligand-target gene probability scores of the top target genes.
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' targets = extract_top_n_targets("BMP2",50,ligand_target_matrix)
#' }
#' @export
#'
extract_top_n_targets = function(ligand_oi,top_n,ligand_target_matrix,ligands_position = "cols"){
  # input check
  if (!is.character(ligand_oi) | length(ligand_oi) != 1)
    stop("ligand_oi must be a character vector of length 1")
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (ligands_position == "cols"){
    ligands = colnames(ligand_target_matrix)
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows"){
    ligands = rownames(ligand_target_matrix)
    targets = colnames(ligand_target_matrix)
  }
  if ((ligand_oi %in% ligands) == FALSE)
    stop("ligand_oi must be in ligand_target_matrix")
  if (!is.numeric(top_n) | length(top_n) != 1 | top_n > length(targets))
    stop("top_n must be a numeric vector of length 1 of which the value is smaller than the total number of target genes")

  requireNamespace("dplyr")
  # ligand should be in matrix
  if (ligands_position == "cols"){
    ligand_oi_vector = ligand_target_matrix[,ligand_oi]
  } else if (ligands_position == "rows"){
    ligand_oi_vector = ligand_target_matrix[ligand_oi,]
  }
  ligand_oi_vector[ligand_oi_vector >= sort(ligand_oi_vector,decreasing = T)[top_n]] %>% sort(decreasing = T)
}
#' @title Get the predicted top n percentage ligands of a target of interest
#'
#' @description \code{extract_top_fraction_ligands} Get the predicted top n percentage ligands of a target of interest
#'
#' @usage
#' extract_top_fraction_ligands(target_oi,top_fraction,ligand_target_matrix,ligands_position = "cols")
#'
#' @param target_oi The target gene of interest of which top upstream ligands should be returned
#' @param top_fraction A number between 0 and 1 indicating which top fraction of ligands should be returned.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A named numeric vector of ligand-target gene probability scores of the top ligands.
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' targets = extract_top_fraction_ligands("ID3",0.01,ligand_target_matrix)
#' }
#' @export
#'
extract_top_fraction_ligands = function(target_oi,top_fraction,ligand_target_matrix,ligands_position = "cols"){
  # input check
  if (!is.character(target_oi) | length(target_oi) != 1)
    stop("target_oi must be a character vector of length 1")
  if (!is.numeric(top_fraction) | length(top_fraction) != 1 | top_fraction > 1 | top_fraction < 0)
    stop("top_fraction_oi must be a numeric vector of length 1 of which the value is between 0 and 1")
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (ligands_position == "cols"){
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows"){
    targets = colnames(ligand_target_matrix)
  }
  if ((target_oi %in% targets) == FALSE)
    stop("target_oi must be in ligand_target_matrix")

  requireNamespace("dplyr")

  if (ligands_position == "cols"){
    target_oi_vector = ligand_target_matrix[target_oi,]
  } else if (ligands_position == "rows"){
    target_oi_vector = ligand_target_matrix[,target_oi]
  }
  target_oi_vector[target_oi_vector >= quantile(target_oi_vector,1-top_fraction)] %>% sort(decreasing = T)
}
#' @title Get the predicted top n ligands of a target gene of interest
#'
#' @description \code{extract_top_n_ligands} Get the predicted top n ligands of a target gene of interest.
#'
#' @usage
#' extract_top_n_ligands(target_oi,top_n,ligand_target_matrix,ligands_position = "cols")
#'
#' @param target_oi The target gene of interest of which top upstream ligands should be returned
#' @param top_n A number between 0 and the total nr of ligands indicating which top n of ligands should be returned.
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A named numeric vector of ligand-target gene probability scores of the top ligands.
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' targets = extract_top_n_ligands("BMP2",2,ligand_target_matrix)
#' }
#' @export
extract_top_n_ligands = function(target_oi,top_n,ligand_target_matrix,ligands_position = "cols"){

  # input check
  if (!is.character(target_oi) | length(target_oi) != 1)
    stop("target_oi must be a character vector of length 1")
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (ligands_position == "cols"){
    ligands = colnames(ligand_target_matrix)
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows"){
    ligands = rownames(ligand_target_matrix)
    targets = colnames(ligand_target_matrix)
  }
  if ((target_oi %in% targets) == FALSE)
    stop("target_oi must be in ligand_target_matrix")
  if (!is.numeric(top_n) | length(top_n) != 1 | top_n > length(ligands))
    stop("top_n must be a numeric vector of length 1 of which the value is smaller than the total number of ligands")

  requireNamespace("dplyr")

  if (ligands_position == "cols"){
    target_oi_vector = ligand_target_matrix[target_oi,]
  } else if (ligands_position == "rows"){
    target_oi_vector = ligand_target_matrix[,target_oi]
  }
  target_oi_vector[target_oi_vector >= sort(target_oi_vector,decreasing = T)[top_n]] %>% sort(decreasing = T)
}
#' @title Get a set of predicted target genes of a ligand of interest
#'
#' @description \code{get_target_genes_ligand_oi} Get a set of predicted target genes of a ligand of interest.
#'
#' @usage
#' get_target_genes_ligand_oi(ligand_oi, ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", fdr_method = "global", output = "logical",ligands_position = "cols")
#'
#' @param ligand_oi The ligand of interest of which top target genes should be returned
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param error_rate FDR for cutoff_method "fdrtool" and "distribution"; number between 0 and 1 indicating which top fraction of target genes should be returned for cutoff_method "quantile". Default: 0.1
#' @param cutoff_method Method to determine which genes can be considered as a target of a ligand and which genes not, based on the ligand-target probability scores. Possible options: "distribution", "fdrtool" and "quantile". Default: "distribution".
#' @param fdr_method Only relevant when cutoff_method is "fdrtool". Possible options: "global" and "local". Default: "global".
#' @param output Determines whether a vector with target gene names should be returned ("gene_symbols") or a logical vector indicating for every target gene whether or not it is a target ("logical").
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A vector with target gene names should be returned ("gene_symbols") or a logical vector indicating for every target gene whether or not it is a target ("logical").
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' targets = get_target_genes_ligand_oi("BMP2", ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", fdr_method = "global", output = "logical",ligands_position = "cols")
#' }
#' @importFrom fdrtool fdrtool
#'
#' @export
#'
get_target_genes_ligand_oi = function(ligand_oi, ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", fdr_method = "global", output = "logical",ligands_position = "cols"){
  # input check
  if (!is.character(ligand_oi) | length(ligand_oi) != 1)
    stop("ligand_oi must be a character vector of length 1")
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if (output != "logical" & output != "gene_symbols")
    stop("output must be 'logical' or 'gene_symbols'")
  if (cutoff_method != "distribution" &  cutoff_method != "fdrtool" & cutoff_method != "quantile")
    stop("cutoff_method must be 'distribution' or 'fdrtool' or 'quantile'")
  if (!is.numeric(error_rate) | length(error_rate) != 1 | error_rate > 1 | error_rate < 0)
    stop("error_rate must be a numeric vector of length 1 of which the value is between 0 and 1")
  if (fdr_method != "global" &  fdr_method != "local")
    stop("fdr_method must be 'global' or 'local'")

  if (ligands_position == "cols"){
    ligands = colnames(ligand_target_matrix)
    targets = rownames(ligand_target_matrix)
  } else if (ligands_position == "rows"){
    ligands = rownames(ligand_target_matrix)
    targets = colnames(ligand_target_matrix)
  }
  if ((ligand_oi %in% ligands) == FALSE)
    stop("ligand_oi must be in ligand_target_matrix")

  requireNamespace("dplyr")

  if (ligands_position == "cols"){
    ligand_target_vector= ligand_target_matrix[,ligand_oi]
  } else if (ligands_position == "rows"){
    ligand_target_vector = ligand_target_matrix[ligand_oi,]
  }
  ## FDRtools method
  if (cutoff_method == "fdrtool"){
    if (!requireNamespace("fdrtool", quietly = TRUE)) {
      stop("Package 'fdrtool' needed for this function. Please install it.")
    }

    lognormal_ligand_target = ligand_target_vector %>% log() %>%  .[is.finite(.)]
    z = lognormal_ligand_target + abs(mean(lognormal_ligand_target)) # to convert log scores to "z-scores" with mean=0
    fdr = fdrtool::fdrtool(z,verbose = F,plot = F)
    if (fdr_method == "local"){
      target_gene_scores = z[fdr$lfdr <= 2*error_rate] %>% .[.>0]  # two-sided --> 2x fdr one-sided (greater than); l
    } else if (fdr_method == "global") {
      target_gene_scores = z[fdr$qval <= 2*error_rate] %>% .[.>0]  # two-sided --> 2x fdr one-sided (greater than); l
    }
    ## Own method based on estimated background distribution
  } else if (cutoff_method == "distribution"){
    lognormal_ligand_target = ligand_target_vector %>% log() %>%  .[is.finite(.)]
    background_dist = rnorm(n = length(lognormal_ligand_target),mean = mean(lognormal_ligand_target), sd = sd(lognormal_ligand_target))
    target_gene_scores = lognormal_ligand_target[lognormal_ligand_target >= quantile(background_dist,1-error_rate, na.rm = TRUE)] %>% sort(decreasing = T)
    ## Quantile-based method: targets = genes in top 1-error rate percentage
  } else if (cutoff_method == "quantile"){
    target_gene_scores = ligand_target_vector[ligand_target_vector >= quantile(ligand_target_vector,1-error_rate, na.rm = TRUE)] %>% sort(decreasing = T)}

  ## Output
  if (output == "gene_symbols"){
    return(names(target_gene_scores))
  } else if (output == "logical") {
    output = targets %in% names(target_gene_scores)
    names(output) = targets
    return(output)
  }

}
