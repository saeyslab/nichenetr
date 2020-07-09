#' @title Construct weighted layer-specific networks
#'
#' @description \code{construct_weighted_networks} construct layer-specific weighted integrated networks from input source networks via weighted aggregation.
#'
#' @usage
#' construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df, n_output_networks = 2)
#'
#' @param lr_network A data frame / tibble containing ligand-receptor interactions (required columns: from, to, source)
#' @param sig_network A data frame / tibble containing signaling interactions (required columns: from, to, source)
#' @param gr_network A data frame / tibble containing gene regulatory interactions (required columns: from, to, source)
#' @param source_weights_df A data frame / tibble containing the weights associated to each individual data source. Sources with higher weights will contribute more to the final model performance (required columns: source, weight). Note that only interactions described by sources included here, will be retained during model construction.
#' @param n_output_networks The number of output networks to return: 2 (ligand-signaling and gene regulatory; default) or 3 (ligand-receptor, signaling and gene regulatory).
#'
#' @return A list containing 2 elements (lr_sig and gr) or 3 elements (lr, sig, gr): the integrated weighted ligand-signaling and gene regulatory networks or ligand-receptor, signaling and gene regulatory networks in data frame / tibble format with columns: from, to, weight.
#'
#'
#' @examples
#' \dontrun{
#' ## Generate the weighted networks from input source networks
#' wn = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' }
#'
#' @export
#'
construct_weighted_networks = function(lr_network, sig_network, gr_network,source_weights_df,n_output_networks = 2) {
  # input check
  if (!is.data.frame(lr_network))
    stop("lr_network must be a data frame or tibble object")
  if (!is.data.frame(sig_network))
    stop("sig_network must be a data frame or tibble object")
  if (!is.data.frame(gr_network))
    stop("gr_network must be a data frame or tibble object")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")
  if (n_output_networks < 2 | n_output_networks > 3)
    stop("The number of required output networks must be 2 (ligand-signaling and regulatory) or 3 (ligand-receptor, signaling and regulatory)")

  requireNamespace("dplyr")

  # remove data sources for which weight equals 0
  source_weights_df = source_weights_df %>% filter(weight > 0)
  # perform weighted network aggregation
  gr_network_w = gr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
  if (n_output_networks == 2) {
    ligand_signaling_w = bind_rows(lr_network, sig_network) %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
    weighted_networks = list(lr_sig = ligand_signaling_w %>% ungroup(), gr = gr_network_w %>% ungroup())
  }
  else if (n_output_networks == 3) {
    lr_network_w = lr_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
    sig_network_w = sig_network %>% inner_join(source_weights_df,by = "source") %>% group_by(from, to) %>% summarize(weight = sum(weight)) %>% ungroup()
    weighted_networks = list(lr = lr_network_w %>% ungroup(), sig = sig_network_w %>% ungroup(), gr = gr_network_w %>% ungroup())
  }
}
#' @title Add a new data source to the model
#'
#' @description \code{add_new_datasource} adds a new data source to one of the ligand-receptor, signaling and gene regulatory data sources.
#'
#' @usage
#' add_new_datasource(new_source, network, new_weight,source_weights_df)
#'
#' @param new_source A data frame / tibble containing novel interactions of the ligand-receptor, signaling or gene regulatory layer (required columns: from, to, source)
#' @param network NULL or a data frame / tibble containing the base network to which you want to add the new data source (required columns: from, to, source)
#' @param new_weight a weight value between 0 and 1 to assign to your new data source
#' @param source_weights_df A data frame / tibble containing the weights associated to each already included individual data source (required columns: source, weight).
#'
#' @return A list containing 2 elements (network and source_weights_df): the updated network containing your data source and the updated source_weights_df containing the weight of the newly added data source.
#'
#'
#' @examples
#'  \dontrun{
#' ## Update the lr_network with a new ligand-receptor data source
#' library(dplyr)
#' lr_toy = tibble(from = "A", to = "B", source = "toy")
#' new_lr_network = add_new_datasource(lr_toy, lr_network,1,source_weights_df)
#' }
#' @export
#'
add_new_datasource = function(new_source, network, new_weight,source_weights_df) {
  # input check
  if (!is.data.frame(new_source))
    stop("new_source must be a data frame or tibble object")
  if (!is.data.frame(network) && !is.null(network))
    stop("network must be a data frame or tibble object")
  if (!is.numeric(new_weight) && length(new_weight) != 1 && (new_weight >= 0 | new_weight <= 1))
    stop("new_weight must be a one number between 0 and 1")
  if (!is.data.frame(source_weights_df) || sum((source_weights_df$weight > 1)) != 0 || new_weight > 1)
    stop("source_weights_df must be a data frame or tibble object and no data source weight may be higher than 1")


  requireNamespace("dplyr")


  if (is.null(network)){
    return(list(network = new_source, source_weights_df = tibble(source = new_source$source %>% unique(),weight = new_weight)))
  }
  network = network %>% bind_rows(new_source)
  source_weights_df = source_weights_df %>% bind_rows(tibble(source = new_source$source %>% unique(),weight = new_weight))

  list(network = network, source_weights_df = source_weights_df)
}
#' @title Apply hub corrections to the weighted integrated ligand-signaling and gene regulatory network
#'
#' @description \code{apply_hub_corrections} downweighs the importance of nodes with a lot of incoming links in the ligand-signaling and/or gene regulatory network. Hub correction method according to following equation: \eqn{Wcor =W * D^-h} with \eqn{D} the indegree matrix of the respective network and \eqn{h} the correction factor.
#'
#' @usage
#' apply_hub_corrections(weighted_networks, lr_sig_hub, gr_hub)
#'
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param lr_sig_hub a number between 0 and 1. 0: no correction for hubiness; 1: maximal correction for hubiness.
#' @param gr_hub a number between 0 and 1. 0: no correction for hubiness;  1: maximal correction for hubiness.
#'
#' @return A list containing 2 elements (lr_sig and gr): the hubiness-corrected integrated weighted ligand-signaling and gene regulatory networks in data frame / tibble format with columns: from, to, weight.
#'
#' @examples
#' \dontrun{
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' wn = apply_hub_corrections(weighted_networks, lr_sig_hub= 0.5, gr_hub= 0.5)
#' }
#' @export
#'
apply_hub_corrections = function(weighted_networks,lr_sig_hub, gr_hub) {
  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

  if (lr_sig_hub < 0 | lr_sig_hub > 1)
    stop("lr_sig_hub must be a number between 0 and 1 (0 and 1 included)")
  if (gr_hub < 0 | gr_hub > 1)
    stop("gr_hub must be a number between 0 and 1 (0 and 1 included)")

  requireNamespace("dplyr")


  # load in weighted networks
  ligand_signaling_network = weighted_networks$lr_sig
  regulatory_network = weighted_networks$gr

  # apply hub correction ligand-signaling network
  if (lr_sig_hub > 0){
    ligand_signaling_network = ligand_signaling_network %>% group_by(to) %>% count(to) %>% ungroup() %>% inner_join(ligand_signaling_network, ., by = "to") %>% group_by(from) %>% mutate(weight = weight/(n**lr_sig_hub)) %>% dplyr::select(-n)
  }
  # apply hub correction gene regulatory network
  if (gr_hub > 0){
    regulatory_network = regulatory_network %>% group_by(to) %>% count(to) %>% ungroup() %>% inner_join(regulatory_network, ., by = "to") %>% group_by(from) %>% mutate(weight = weight/(n**gr_hub)) %>% dplyr::select(-n)
  }
  return(list(lr_sig = ligand_signaling_network %>% ungroup(), gr = regulatory_network %>% ungroup()))
}
#' @title Construct a ligand-tf signaling probability matrix for ligands of interest.
#'
#' @description \code{construct_ligand_tf_matrix} Convert integrated weighted networks into a matrix containg ligand-tf probability scores. The higher this score, the more likely a particular ligand can signal to a downstream gene.
#'
#' @usage
#' construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, ligands_as_cols = FALSE)
#'
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param ligands A list of all ligands and ligand-combinations of which target gene probability scores should be calculated. Example format: list("TNF","BMP2",c("IL4","IL13")).
#' @param ltf_cutoff Ligand-tf scores beneath the "ltf_cutoff" quantile will be set to 0. Default: 0.99 such that only the 1 percent closest tfs will be considered as possible tfs downstream of the ligand of choice.
#' @param algorithm Selection of the algorithm to calculate ligand-tf signaling probability scores. Different options: "PPR" (personalized pagerank), "SPL" (shortest path length) and "direct"(just take weights of ligand-signaling network as ligand-tf weights + give the ligand itself the max score). Default and recommended: PPR.
#' @param damping_factor Only relevant when algorithm is PPR. In the PPR algorithm, the damping factor is the probability that the random walker will continue its walk on the graph; 1-damping factor is the probability that the walker will return to the seed node. Default: 0.5.
#' @param ligands_as_cols Indicate whether ligands should be in columns of the matrix and target genes in rows or vice versa. Default: FALSE
#'
#' @return A matrix containing ligand-target probability scores.
#' @importFrom igraph page_rank V graph_from_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#'  \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_tf = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5,ligands_as_cols = TRUE)
#' }
#' @export
#'
construct_ligand_tf_matrix = function(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, ligands_as_cols = FALSE) {

  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

  if (!is.list(ligands))
    stop("ligands must be a list object")
  if ( sum((unique(unlist(ligands)) %in% unique(c(lr_network$from,lr_network$to))) == FALSE) > 0)
    warning("One or more ligands of interest not present in the ligand-receptor network 'lr_network'. You can possibly ignore this warning if you provided your own ligand_receptor network to the weighted networks." )

  if(is.null(ltf_cutoff)){
    if( algorithm == "PPR" | algorithm == "SPL" )
      warning("Did you not forget to give a value to ltf_cutoff?")
  } else {
    if (ltf_cutoff < 0 | ltf_cutoff > 1)
      stop("ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }

  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(algorithm == "PPR"){
    if (damping_factor < 0 | damping_factor >= 1)
      stop("damping_factor must be a number between 0 and 1 (0 included, 1 not)")
    }
  if (!is.logical(ligands_as_cols) | length(ligands_as_cols) != 1)
    stop("ligands_as_cols must be a logical vector of length 1")

  requireNamespace("dplyr")

  # load in weighted networks
  ligand_signaling_network = weighted_networks$lr_sig
  regulatory_network = weighted_networks$gr

  # convert ids to numeric for making Matrix::sparseMatrix later on
  allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")

  ligand_signaling_network = ligand_signaling_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)

  if (algorithm == "PPR"){
    # Make Matrix::sparse signaling weighted matrix and graph to apply personalized pagerank
    ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
    signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
  # personalized pagerank
  # prepare preference vector E
  E = rep(0,times = length(igraph::V(signaling_igraph)))
  # ppr for every ligand individual
  complete_matrix = lapply(ligands,PPR_wrapper,E,signaling_igraph,damping_factor,id2allgenes,ltf_cutoff)
  } else if (algorithm == "SPL"){
    # Make Matrix::sparse signaling weighted matrix and graph to apply shortest path length
    ligand_signaling_network = ligand_signaling_network %>% mutate(weight = 1/weight) # reverse the weight because spl: find shortest path, with lowest weight, so original weights needed to be reversed
    ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
    signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")
    # SPL
    complete_matrix = lapply(ligands,SPL_wrapper,signaling_igraph,id2allgenes,ltf_cutoff)
  } else if (algorithm == "direct"){
    ligand_ligand_network = tibble::tibble(from_allgenes = id2allgenes[unlist(ligands)], to_allgenes = id2allgenes[unlist(ligands)]) %>% inner_join(ligand_signaling_network %>% filter(from_allgenes %in% id2allgenes[unlist(ligands)]) %>% group_by(from_allgenes) %>% top_n(1,weight) %>% ungroup() %>% distinct(from_allgenes,weight), by = "from_allgenes")
    ligand_signaling_network = ligand_signaling_network %>% bind_rows(ligand_ligand_network)
    ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
    signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")

    complete_matrix = lapply(ligands,direct_wrapper,ligand_signaling_network_matrix,id2allgenes,ltf_cutoff)
  }
  ltf_matrix = matrix(unlist(complete_matrix), ncol = length(igraph::V(signaling_igraph)), byrow = TRUE)

  rownames(ltf_matrix) = sapply(ligands,function(x){paste0(x,collapse = "-")})
  colnames(ltf_matrix) = allgenes

  if (ligands_as_cols == TRUE){
    tf_matrix = ltf_matrix %>% as.matrix()
    ltf_matrix = ltf_matrix %>% t()
  }

  return(ltf_matrix)
}
#' @title Construct a tf-target matrix.
#'
#' @description \code{construct_tf_target_matrix} Convert integrated gene regulatory weighted network into matrix format.
#'
#' @usage
#' construct_tf_target_matrix(weighted_networks, tfs_as_cols = FALSE, standalone_output = FALSE)
#'
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param tfs_as_cols Indicate whether ligands should be in columns of the matrix and target genes in rows or vice versa. Default: FALSE
#' @param standalone_output Indicate whether the ligand-tf matrix should be formatted in a way convenient to use alone (with gene symbols as row/colnames). Default: FALSE
#'
#' @return A matrix containing tf-target regulatory weights.
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#'  \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' tf_target = construct_tf_target_matrix(weighted_networks, tfs_as_cols = TRUE, standalone_output = TRUE)
#' }
#' @export
#'
construct_tf_target_matrix = function(weighted_networks, tfs_as_cols = FALSE, standalone_output = FALSE) {
  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

  if (!is.logical(tfs_as_cols) | length(tfs_as_cols) != 1)
    stop("ligands_as_cols must be a logical vector of length 1")
  if (!is.logical(standalone_output) | length(standalone_output) != 1)
    stop("standalone_output must be a logical vector of length 1")
  requireNamespace("dplyr")

  # load in weighted networks
  ligand_signaling_network = weighted_networks$lr_sig
  regulatory_network = weighted_networks$gr

  # convert ids to numeric for making Matrix::sparseMatrix later on
  allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")

  regulatory_network = regulatory_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)
  grn_matrix = Matrix::sparseMatrix(regulatory_network$from_allgenes %>% as.integer, regulatory_network$to_allgenes %>% as.integer, x=regulatory_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))

  rownames(grn_matrix) = allgenes
  colnames(grn_matrix) = allgenes
  if (standalone_output == TRUE){
    # keep only regulators with gene regulatory interactions
    regulators = weighted_networks$gr %>% .$from %>% unique()
    grn_matrix = grn_matrix[regulators,]

  }
  if (tfs_as_cols == TRUE){
    grn_matrix = grn_matrix %>% as.matrix()
    grn_matrix = grn_matrix %>% t()
  }

  return(grn_matrix)
}
#' @title Construct a ligand-target probability matrix for ligands of interest.
#'
#' @description \code{construct_ligand_target_matrix} Convert integrated weighted networks into a matrix containg ligand-target probability scores. The higher this score, the more likely a particular ligand can induce the expression of a particular target gene.
#'
#' @usage
#' construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE,ligands_as_cols = TRUE, remove_direct_links = "no")
#'
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param ligands A list of all ligands and ligand-combinations of which target gene probability scores should be calculated. Example format: list("TNF","BMP2",c("IL4","IL13")).
#' @param ltf_cutoff Ligand-tf scores beneath the "ltf_cutoff" quantile will be set to 0. Default: 0.99 such that only the 1 percent closest tfs will be considered as possible tfs downstream of the ligand of choice.
#' @param algorithm Selection of the algorithm to calculate ligand-tf signaling probability scores. Different options: "PPR" (personalized pagerank), "SPL" (shortest path length) and "direct"(just take weights of ligand-signaling network as ligand-tf weights). Default and recommended: PPR.
#' @param damping_factor Only relevant when algorithm is PPR. In the PPR algorithm, the damping factor is the probability that the random walker will continue its walk on the graph; 1-damping factor is the probability that the walker will return to the seed node. Default: 0.5.
#' @param secondary_targets Indicate whether a ligand-target matrix should be returned that explicitly includes putative secondary targets of a ligand (by means of an additional matrix multiplication step considering primary targets as possible regulators). Default: FALSE
#' @param ligands_as_cols Indicate whether ligands should be in columns of the matrix and target genes in rows or vice versa. Default: TRUE
#' @param remove_direct_links Indicate whether direct ligand-target and receptor-target links in the gene regulatory network should be kept or not. "no": keep links; "ligand": remove direct ligand-target links; "ligand-receptor": remove both direct ligand-target and receptor-target links. Default: "no"
#'
#' @return A matrix containing ligand-target probability scores.
#' @importFrom igraph page_rank V graph_from_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#'  \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE, remove_direct_links = "no")
#' }
#' @export
#'
construct_ligand_target_matrix = function(weighted_networks, ligands, ltf_cutoff = 0.99, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE,ligands_as_cols = TRUE, remove_direct_links = "no") {

  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

  if (!is.list(ligands))
    stop("ligands must be a list object")

  if(is.null(ltf_cutoff)){
    if( algorithm == "PPR" | algorithm == "SPL" )
      warning("Did you not forget to give a value to ltf_cutoff?")
  } else {
    if (ltf_cutoff < 0 | ltf_cutoff > 1)
      stop("ltf_cutoff must be a number between 0 and 1 (0 and 1 included)")
  }

  if (algorithm != "PPR" & algorithm != "SPL" & algorithm != "direct")
    stop("algorithm must be 'PPR' or 'SPL' or 'direct'")
  if(algorithm == "PPR"){
    if (damping_factor < 0 | damping_factor >= 1)
      stop("damping_factor must be a number between 0 and 1 (0 included, 1 not)")
  }
  if (!is.logical(secondary_targets) | length(secondary_targets) != 1)
    stop("secondary_targets must be a logical vector of length 1")
  if (!is.logical(ligands_as_cols) | length(ligands_as_cols) != 1)
    stop("ligands_as_cols must be a logical vector of length 1")
  if (remove_direct_links != "no" & remove_direct_links != "ligand" & remove_direct_links != "ligand-receptor")
    stop("remove_direct_links must be 'no' or 'ligand' or 'ligand-receptor'")

  requireNamespace("dplyr")

  ## workflow; first: give probability score to genes downstream in signaling path starting from ligand. second: multiply this matrix with gene regulatory matrix to get the probabilty scores of downstream target genes

  # construct ligand-tf matrix
  # remove direct links if required
  if (remove_direct_links != "no"){
    ligands_ = lr_network$from %>% unique()
    receptors_ = lr_network$to %>% unique()
    if (remove_direct_links == "ligand"){
      weighted_networks$gr = weighted_networks$gr %>% filter((from %in% ligands_) == FALSE)
    } else if (remove_direct_links == "ligand-receptor"){
      weighted_networks$gr = weighted_networks$gr %>% filter((from %in% c(ligands_,receptors_)) == FALSE)
    }
  }

  ltf_matrix = construct_ligand_tf_matrix(weighted_networks, ligands, ltf_cutoff, algorithm, damping_factor)

  # preparing the gene regulatory matrix
  grn_matrix = construct_tf_target_matrix(weighted_networks)

  # Multiply ligand-tf matrix with tf-target matrix
  ligand_to_target = (ltf_matrix %*% grn_matrix)

  # Secondary targets
  if (secondary_targets == TRUE) {
    ltf_matrix = ligand_to_target

    if (ltf_cutoff > 0){
      ltf_matrix_TRUE = apply(ltf_matrix,1,function(x){x <= quantile(x,ltf_cutoff)}) %>% t()
      ltf_matrix[ltf_matrix_TRUE] = 0
    }

    ligand_to_target_primary = ligand_to_target
    ligand_to_target_secondary = ltf_matrix %*% grn_matrix

    # set 0's to number higher than 0 to avoid Inf when inverting in sum
    ligand_to_target_primary[ligand_to_target_primary == 0] = Inf
    ligand_to_target_primary[is.infinite(ligand_to_target_primary)] = min(ligand_to_target_primary)

    ligand_to_target_secondary[ligand_to_target_secondary == 0] = Inf
    ligand_to_target_secondary[is.infinite(ligand_to_target_secondary)] = min(ligand_to_target_secondary)

    # inverting in sum to emphasize primary targets more (scores secondary targets matrix tend to be higher )
    ligand_to_target = (ligand_to_target_primary **-1 + ligand_to_target_secondary **-1) ** -1
  }

  ligand_to_target = ligand_to_target %>% as.matrix()

  if (ligands_as_cols == TRUE){
    ligand_to_target = ligand_to_target %>% t()
  }

  return(ligand_to_target)
}
#' @title Adapt a ligand-target probability matrix construced via PPR by correcting for network topolgoy.
#'
#' @description \code{correct_topology_ppr} The ligand-target probability scores of a matrix constructed via personalized pagerank will be subtracted by target probability scores calculated via global pagerank; these latter scores can be considered as scores solely attributed to network topology and not by proximity to the ligand of interest. Recommended to use this function in combination with a ligand-target matrix constructed without applying a cutoff on the ligand-tf matrix.
#'
#' @usage
#' correct_topology_ppr(ligand_target_matrix,weighted_networks,ligands_position = "cols")
#'
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param weighted_networks A list of two elements: lr_sig: a data frame/ tibble containg weighted ligand-receptor and signaling interactions (from, to, weight); and gr: a data frame/tibble containng weighted gene regulatory interactions (from, to, weight)
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A matrix containing ligand-target probability scores, after subtracting target scores solely due to network topology.
#' @importFrom igraph page_rank V graph_from_adjacency_matrix
#' @importFrom Matrix sparseMatrix
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands, ltf_cutoff = 0, algorithm = "PPR", damping_factor = 0.5, secondary_targets = FALSE)
#' ligand_target_matrix = correct_topology_ppr(ligand_target_matrix,weighted_networks,ligands_position = "cols")
#' }
#' @export
#'
correct_topology_ppr = function(ligand_target_matrix,weighted_networks,ligands_position = "cols"){
  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  if (!is.numeric(weighted_networks$lr_sig$weight))
    stop("lr_sig must contain a column named data source weights")
  if (!is.numeric(weighted_networks$gr$weight))
    stop("gr must contain a column named data source weights")

  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")

  requireNamespace("dplyr")

  background_pr = get_pagerank_target(weighted_networks)

  if (ligands_position == "rows"){
    ligand_target_matrix_new = apply(ligand_target_matrix,1,function(x){x - background_pr}) %>% t()
  } else {
    ligand_target_matrix_new = apply(ligand_target_matrix,2,function(x){x - background_pr})
  }

  ligand_target_matrix_new[ligand_target_matrix_new < 0] = 0

  return(ligand_target_matrix_new)
}
#' @title Convert probabilistic ligand-target matrix to a discrete one.
#'
#' @description \code{make_discrete_ligand_target_matrix} Convert probabilistic ligand-target matrix to a discrete one. This means that for every ligand, genes will be labeled as target (TRUE) or non-target (FALSE).
#'
#' @usage
#' make_discrete_ligand_target_matrix(ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", fdr_method = "global",ligands_position = "cols")
#'
#' @param ligand_target_matrix A matrix of ligand-target probabilty scores.
#' @param error_rate FDR for cutoff_method "fdrtool" and "distribution"; number between 0 and 1 indicating which top fraction of target genes should be returned for cutoff_method "quantile". Default: 0.1
#' @param cutoff_method Method to determine which genes can be considered as a target of a ligand and which genes not, based on the ligand-target probability scores. Possible options: "distribution", "fdrtool" and "quantile". Default: "distribution".
#' @param fdr_method Only relevant when cutoff_method is "fdrtool". Possible options: "global" and "local". Default: "global".
#' @param ligands_position Indicate whether the ligands in the ligand-target matrix are in the rows ("rows") or columns ("cols"). Default: "cols"
#'
#' @return A matrix of ligand-target assignments. TRUE: gene is a target of the ligand of interest; FALSE: gene is not a target of the ligand of interest.
#'
#' @examples
#' \dontrun{
#' ## Generate the ligand-target matrix from loaded weighted_networks
#' weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
#' ligands = list("TNF","BMP2",c("IL4","IL13"))
#' ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)
#' ligand_target_matrix = make_discrete_ligand_target_matrix(ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", ligands_position = "cols")
#'}
#' @export
#'
make_discrete_ligand_target_matrix = function(ligand_target_matrix, error_rate = 0.1, cutoff_method = "distribution", fdr_method = "global",ligands_position = "cols"){
  # input check
  if (!is.matrix(ligand_target_matrix))
    stop("ligand_target_matrix must be a matrix object")
  if (ligands_position != "cols" & ligands_position != "rows")
    stop("ligands_position must be 'cols' or 'rows'")
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
    ligand_target_matrix = ligand_target_matrix %>% t()
  }

  requireNamespace("dplyr")

  list_targets = lapply(ligands,get_target_genes_ligand_oi,ligand_target_matrix,cutoff_method = cutoff_method, fdr_method = fdr_method,error_rate = error_rate)
  names(list_targets) = ligands

  ligand_target_matrix_discrete = list_targets %>% bind_rows() %>% as.matrix()

  if(nrow(ligand_target_matrix_discrete) == length(targets) & ncol(ligand_target_matrix_discrete) == length(ligands)){
    ligand_target_matrix_discrete = ligand_target_matrix_discrete
  } else if (nrow(ligand_target_matrix_discrete) == length(ligands) & ncol(ligand_target_matrix_discrete) == length(targets)) {
    ligand_target_matrix_discrete = ligand_target_matrix_discrete %>% t()
  }


  if (ligands_position == "cols"){
    rownames(ligand_target_matrix_discrete) = targets
    colnames(ligand_target_matrix_discrete) = ligands
  } else if (ligands_position == "rows"){
    ligand_target_matrix_discrete = ligand_target_matrix_discrete %>% t()
    colnames(ligand_target_matrix_discrete) = targets
    rownames(ligand_target_matrix_discrete) = ligands
  }

  return(ligand_target_matrix_discrete)
}




