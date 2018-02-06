## introduction to package, which functions to import,...

#' nichenetr: Linking Extracellular Protein Signals to Target Genes by data-integration.
#'
#' This package allows you the investigate intercellular communication from a computational perspective. Functionalities of this package (e.g. including predicting extracellular upstream regulators) build upon a probabilistic model of ligand-target links that was inferred by data-integration.
#'
#' @section Construction of the probabilistic model:
#' \code{\link{construct_weighted_networks}}, \code{\link{construct_ligand_target_matrix}}
#'
#' @section Evaluation functions:
#' \code{\link{evaluate_target_prediction}}, \code{\link{evaluate_ligand_prediction}}
#'
#' @section Application functions:
#' \code{\link{prioritize_ligands}}, \code{\link{visualize_signaling}}
#'
#' @docType package
#' @name nichenetr
#'
#' @import tidyverse
#' @import purrr
#' @import magrittr
NULL
