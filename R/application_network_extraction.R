#' @title Get active ligand-receptor network for cellular interaction between a sender and receiver cell.
#'
#' @description \code{get_active_ligand_receptor_network} Get active ligand-receptor network by looking for which ligands are expressed in a sender/signaling cell and which receptors are expressed in the receiver cell. Instead of looking at absolute expression, it is possible as well to extract a ligand-receptor network of differentially expressed ligands and receptors if a vector of log2 fold change values is used as input.
#'
#' @usage
#' get_active_ligand_receptor_network(expression_sender, expression_receiver, lr_network, expression_cutoff_sender = 0, expression_cutoff_receiver = 0)
#'
#' @param expression_sender A named numeric vector of gene expression levels (absolute or logfc) for the signaling cell that sends extracellular signals to the receiver cell
#' @param expression_receiver A named numeric vector of gene expression levels (absolute or logfc) for the receiver cell that receives extracellular signals from the sender cell
#' @param lr_network A data frame / tibble containing ligand-receptor interactions (required columns: from, to). Can be both unweighted and weighted.
#' @param expression_cutoff_sender The cutoff on expression value for the sender cell: ligands will be considered active if their expression is higher than the cutoff. Default: 0.
#' @param expression_cutoff_receiver The cutoff on expression value for the receiver cell: receptors will be considered active if their expression is higher than the cutoff. Default: 0.
#'
#' @return A data frame containing at least the variables from, to, sender_expression, receiver_expression. In this network, the active ligand-receptor interactions in the system of interest are shown.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' expression_vector_sender = rnorm(n = 10000, mean = 6, sd = 3)
#' expression_vector_receiver = rnorm(n = 10000, mean = 6, sd = 3)
#' names(expression_vector_sender) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)
#' names(expression_vector_receiver) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)
#' weighted_lr_network = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df, n_output_networks = 3) %>% .$lr
#' sender_cell_receiver_lr_network = get_active_ligand_receptor_network(expression_vector_sender,expression_vector_receiver,weighted_lr_network,expression_cutoff_sender = 4, expression_cutoff_receiver = 4)
#' }
#' @export
#'
get_active_ligand_receptor_network = function(expression_sender, expression_receiver, lr_network, expression_cutoff_sender = 0, expression_cutoff_receiver = 0){

  if(!is.numeric(expression_sender) | is.null(names(expression_sender)))
    stop("expression_sender should be a named numeric vector indicating the expression value for every gene")
  if(!is.numeric(expression_receiver) | is.null(names(expression_receiver)))
    stop("expression_receiver should be a named numeric vector indicating the expression value for every gene")
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if(!is.numeric(expression_cutoff_sender) | length(expression_cutoff_sender) != 1 )
    stop("expression_cutoff_sender should be a single number")
  if(!is.numeric(expression_cutoff_receiver) | length(expression_cutoff_receiver) != 1 )
    stop("expression_cutoff_receiver should be a single number")

  requireNamespace("dplyr")

  sender_expression_df = tibble(from = names(expression_sender), sender_expression = expression_sender)
  receiver_expression_df = tibble(to = names(expression_receiver), receiver_expression = expression_receiver)

  active_lr_network = lr_network %>% inner_join(sender_expression_df, by = "from") %>% inner_join(receiver_expression_df, by = "to") %>% filter(sender_expression > expression_cutoff_sender & receiver_expression > expression_cutoff_receiver)

  return(active_lr_network)
}
#' @title Get active signaling network in a receiver cell.
#'
#' @description \code{get_active_signaling_network} Get active signaling network by looking for which signaling-related genes are expressed in the receiver cell. Instead of looking at absolute expression, it is possible as well to extract a signaling network of differentially expressed signaling mediators if a vector of log2 fold change values is used as input.
#'
#' @usage
#' get_active_signaling_network(expression_receiver, sig_network, expression_cutoff_receiver = 0)
#'
#' @param sig_network A data frame / tibble containing signaling interactions (required columns: from, to). Can be both unweighted and weighted.
#' @inheritParams get_active_ligand_receptor_network
#'
#' @return A data frame containing at least the variables from, to, receiver_expression. In this network, the active signaling interactions in the system of interest are shown.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' expression_vector_receiver = rnorm(n = 10000, mean = 6, sd = 3)
#' names(expression_vector_receiver) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)
#' receiver_sig_network = get_active_signaling_network(expression_vector_receiver,sig_network,expression_cutoff_receiver = 4)
#' }
#' @export
#'
get_active_signaling_network = function(expression_receiver, sig_network, expression_cutoff_receiver = 0){

  if(!is.numeric(expression_receiver) | is.null(names(expression_receiver)))
    stop("expression_receiver should be a named numeric vector indicating the expression value for every gene")
  if (!is.data.frame(sig_network))
    stop("lr_network must be a data frame or tibble object")
  if(!is.numeric(expression_cutoff_receiver) | length(expression_cutoff_receiver) != 1 )
    stop("expression_cutoff_receiver should be a single number")

  requireNamespace("dplyr")

  receiver_expression_df = tibble(from = names(expression_receiver), receiver_expression = expression_receiver)
  receiver_expression_df2 = tibble(to = names(expression_receiver), receiver_expression2 = expression_receiver)

  active_sig_network = sig_network %>% inner_join(receiver_expression_df, by = "from") %>% inner_join(receiver_expression_df2, by = "to") %>% filter(receiver_expression > expression_cutoff_receiver & receiver_expression2 > expression_cutoff_receiver)

  return(active_sig_network)
}
#' @title Get active gene regulatory network in a receiver cell.
#'
#' @description \code{get_active_regulatory_network} Get active regulatory network by looking for which TFs/regulatory genes are expressed in the receiver cell. Instead of looking at absolute expression, it is possible as well to extract a gene regulatory network of differentially expressed TFs if a vector of log2 fold change values is used as input.
#'
#' @usage
#' get_active_regulatory_network(expression_receiver, gr_network, expression_cutoff_receiver = 0)
#'
#' @param gr_network A data frame / tibble containing gene regulatory interactions (required columns: from, to). Can be both unweighted and weighted.
#' @inheritParams get_active_ligand_receptor_network
#'
#' @return A data frame containing at least the variables from, to, receiver_expression. In this network, the active gene regulatory interactions in the system of interest are shown.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' expression_vector_receiver = rnorm(n = 10000, mean = 6, sd = 3)
#' names(expression_vector_receiver) = sample(x = geneinfo_human$symbol,size = 10000,replace = FALSE)
#' receiver_gr_network = get_active_regulatory_network(expression_vector_receiver,gr_network,expression_cutoff_receiver = 4)
#' }
#' @export
#'
get_active_regulatory_network = function(expression_receiver, gr_network, expression_cutoff_receiver = 0){

  if(!is.numeric(expression_receiver) | is.null(names(expression_receiver)))
    stop("expression_receiver should be a named numeric vector indicating the expression value for every gene")
  if (!is.data.frame(gr_network))
    stop("lr_network must be a data frame or tibble object")
  if(!is.numeric(expression_cutoff_receiver) | length(expression_cutoff_receiver) != 1 )
    stop("expression_cutoff_receiver should be a single number")

  requireNamespace("dplyr")

  receiver_expression_df = tibble(from = names(expression_receiver), receiver_expression = expression_receiver)

  active_gr_network = gr_network %>% inner_join(receiver_expression_df, by = "from") %>%  filter(receiver_expression > expression_cutoff_receiver)

  return(active_gr_network)
}
#' @title Get active ligand-target matrix.
#'
#' @description \code{get_active_ligand_target_matrix} Get active ligand-target matrix, meaning that only target genes part of the response of interest will be kept as target genes in the input ligand-target matrix.
#'
#' @usage
#' get_active_ligand_target_matrix(response,ligand_target_matrix, ligands_position = "cols")
#'
#' @param response A named logical vector indicating whether a gene is responding in a biological system or not (e.g. DE after cell-cell interaction)
#' @inheritParams evaluate_target_prediction
#'
#' @return A matrix with ligand-target probability scores (or discrete ligand-target assignments) for the active target genes in the system of interest.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = lapply(expression_settings_validation[1:2],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(setting)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' active_lt = get_active_ligand_target_matrix(setting[[1]] %>% .$response, ligand_target_matrix)
#' }
#' @export
#'
get_active_ligand_target_matrix = function(response,ligand_target_matrix, ligands_position = "cols"){

  # input check
  if(!is.logical(response) | is.null(names(response)))
    stop("response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")

  requireNamespace("dplyr")

  affected_targets = response %>% .[. == TRUE] %>% names()
  if (ligands_position == "cols"){
    affected_targets = affected_targets[affected_targets %in% rownames(ligand_target_matrix)]
    active_ligand_target = ligand_target_matrix[affected_targets,]
  } else if (ligands_position == "rows") {
    affected_targets = affected_targets[affected_targets %in% colnames(ligand_target_matrix)]
    active_ligand_target = ligand_target_matrix[,affected_targets] %>% t()
  }
  return(active_ligand_target)
}
#' @title Get active ligand-target network in data frame format.
#'
#' @description \code{get_active_ligand_target_df} Get active ligand-target network, meaning that only target genes that are part of the response of interest will be kept as target genes. In addition, target genes with probability scores beneath a predefined cutoff will be removed from the nework.
#'
#' @usage
#' get_active_ligand_target_df(response,ligand_target_matrix, ligands_position = "cols", cutoff = 0)
#'
#' @param cutoff A number indicating how high ligand-target probability scores should be in order to be kept. Default: 0. When 0 is used as cutoff and the input matrix is a discrete ligand-target assingment matrix, only TRUE interactions will be kept.
#' @inheritParams get_active_ligand_target_matrix
#'
#' @return A data frame representing the active ligand-target network; with variables $ligand, $target and $score.
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df)
#' setting = lapply(expression_settings_validation[1:2],convert_expression_settings_evaluation)
#' ligands = extract_ligands_from_settings(setting)
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_target_matrix_discrete = make_discrete_ligand_target_matrix(ligand_target_matrix)
#' active_lt_df = get_active_ligand_target_df(setting[[1]] %>% .$response, ligand_target_matrix_discrete)
#' }
#' @export
#'
get_active_ligand_target_df = function(response,ligand_target_matrix, ligands_position = "cols", cutoff = 0){

  # input check
  if(!is.logical(response) | is.null(names(response)))
    stop("response should be named logical vector containing class labels of the response that needs to be predicted ")
  if(!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix should be a matrix")
  if(!is.double(ligand_target_matrix) & !is.logical(ligand_target_matrix))
    stop("ligand_target matrix should be of type double if it contains numeric probabilities as predictions; or of type logical when it contains categorical target predictions (TRUE or FALSE)")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
  if(!is.numeric(cutoff) | length(cutoff) != 1 )
    stop("cutoff should be a single number")

  requireNamespace("dplyr")

  active_ligand_target = get_active_ligand_target_matrix(response, ligand_target_matrix, ligands_position)
  ligand_target_network = active_ligand_target %>% data.frame() %>% rownames_to_column("target") %>%  gather("ligand","score",-target) %>% as_tibble() %>% select(ligand,target,score) %>% filter(score > cutoff)
  return(ligand_target_network)
}

