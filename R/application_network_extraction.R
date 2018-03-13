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
#' library(Biobase)
#' library(dplyr)
#' mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")
#' expression_vector_sender = Exprs_lsec[,Exprs_lsec$celltype == "LSEC_12h"] %>% apply(1,mean)
#' expression_vector_receiver = Exprs_mono_kc[,Exprs_mono_kc$celltype == "BM_mono"] %>% apply(1,mean)
#' names(expression_vector_sender) = names(expression_vector_sender) %>% mousesymbol2humansymbol[.]
#' names(expression_vector_receiver) = names(expression_vector_receiver) %>% mousesymbol2humansymbol[.]
#' expression_vector_sender = expression_vector_sender %>% .[!is.na(names(.))]
#' expression_vector_receiver = expression_vector_receiver %>% .[!is.na(names(.))]
#' weighted_lr_network = construct_weighted_networks(lr_network, sig_network, gr_network, source_weights_df, n_output_networks = 3) %>% .$lr
#' lsec_mono_lr_network = get_active_ligand_receptor_network(expression_vector_sender,expression_vector_receiver,weighted_lr_network,expression_cutoff_sender = 0, expression_cutoff_receiver = 4)
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
#' library(Biobase)
#' library(dplyr)
#' mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")
#' expression_vector_receiver = Exprs_mono_kc[,Exprs_mono_kc$celltype == "BM_mono"] %>% apply(1,mean)
#' names(expression_vector_receiver) = names(expression_vector_receiver) %>% mousesymbol2humansymbol[.]
#' expression_vector_receiver = expression_vector_receiver %>% .[!is.na(names(.))]
#' mono_sig_network = get_active_signaling_network(expression_vector_receiver,sig_network,expression_cutoff_receiver = 4)
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
#' get_active_regulatory_network(expression_receiver, sig_network, expression_cutoff_receiver = 0)
#'
#' @param gr_network A data frame / tibble containing gene regulatory interactions (required columns: from, to). Can be both unweighted and weighted.
#' @inheritParams get_active_ligand_receptor_network
#'
#' @return A data frame containing at least the variables from, to, receiver_expression. In this network, the active gene regulatory interactions in the system of interest are shown.
#'
#' @examples
#' library(Biobase)
#' library(dplyr)
#' mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")
#' expression_vector_receiver = Exprs_mono_kc[,Exprs_mono_kc$celltype == "BM_mono"] %>% apply(1,mean)
#' names(expression_vector_receiver) = names(expression_vector_receiver) %>% mousesymbol2humansymbol[.]
#' expression_vector_receiver = expression_vector_receiver %>% .[!is.na(names(.))]
#' mono_gr_network = get_active_regulatory_network(expression_vector_receiver,gr_network,expression_cutoff_receiver = 4)
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


