# geneinfo_human = load("data/geneinfo_human.rda")
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
# entrez2symbol = mapper(geneinfo_human, "symbol_mouse", "entrez_mouse")
# symbol2entrez = mapper(geneinfo_human, "entrez_mouse", "symbol_mouse")
# humanentrez2humansymbol = mapper(geneinfo_human,"symbol","entrez")
# humansymbol2humanentrez = mapper(geneinfo_human, "entrez", "symbol")
# mouseentrez2humansymbol = mapper(geneinfo_human, "symbol","entrez_mouse")
# humanentrez2mousesymbol = mapper(geneinfo_human, "symbol_mouse", "entrez")
# humansymbol2mouseentrez = mapper(geneinfo_human, "entrez_mouse", "symbol")
# mousesymbol2humanentrez = mapper(geneinfo_human,"entrez","symbol_mouse")
# humansymbol2mousesymbol = mapper(geneinfo_human,"symbol_mouse","symbol")
# mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")

#' @title Convert human gene symbols to their mouse one-to-one orthologs.
#'
#' @description \code{convert_human_to_mouse_symbols} Convert human gene symbols to their mouse one-to-one orthologs.
#'
#' @usage
#' convert_human_to_mouse_symbols(symbols)
#'
#' @param symbols A character vector of official human gene symbols
#'
#' @return A character vector of official mouse gene symbols (one-to-one orthologs of the input human gene symbols).
#'
#' @examples
#' library(dplyr)
#' human_symbols = c("TNF","IFNG")
#' mouse_symbols = human_symbols %>% convert_human_to_mouse_symbols()
#' @export
#'
convert_human_to_mouse_symbols = function(symbols){

  if(!is.character(symbols))
    stop("symbols should be a character vector of human gene symbols")

  requireNamespace("dplyr")

  unambiguous_mouse_genes = geneinfo_human %>% drop_na() %>% group_by(symbol_mouse) %>% count() %>% filter(n<2) %>% .$symbol_mouse
  ambiguous_mouse_genes = geneinfo_human %>% drop_na() %>% group_by(symbol_mouse) %>% count() %>% filter(n>=2) %>% .$symbol_mouse

  geneinfo_ambiguous_solved = geneinfo_human %>% filter(symbol_mouse %in% ambiguous_mouse_genes) %>% filter(symbol==toupper(symbol_mouse))
  geneinfo_human = geneinfo_human %>% filter(symbol_mouse %in% unambiguous_mouse_genes) %>% bind_rows(geneinfo_ambiguous_solved) %>% drop_na()

  humansymbol2mousesymbol = mapper(geneinfo_human,"symbol_mouse","symbol")

  converted_symbols = symbols %>% humansymbol2mousesymbol[.]

  return(converted_symbols)
}
#' @title Convert mouse gene symbols to their human one-to-one orthologs.
#'
#' @description \code{convert_mouse_to_human_symbols} Convert mouse gene symbols to their human one-to-one orthologs.
#'
#' @usage
#' convert_mouse_to_human_symbols(symbols)
#'
#' @param symbols A character vector of official mouse gene symbols
#'
#' @return A character vector of official human gene symbols (one-to-one orthologs of the input mouse gene symbols).
#'
#' @examples
#' library(dplyr)
#' mouse_symbols = c("Tnf","Ifng")
#' human_symbols = mouse_symbols %>% convert_mouse_to_human_symbols()
#' @export
#'
convert_mouse_to_human_symbols = function(symbols){

  if(!is.character(symbols))
    stop("symbols should be a character vector of mouse gene symbols")

  requireNamespace("dplyr")

  unambiguous_mouse_genes = geneinfo_human %>% drop_na() %>% group_by(symbol_mouse) %>% count() %>% filter(n<2) %>% .$symbol_mouse
  ambiguous_mouse_genes = geneinfo_human %>% drop_na() %>% group_by(symbol_mouse) %>% count() %>% filter(n>=2) %>% .$symbol_mouse

  geneinfo_ambiguous_solved = geneinfo_human %>% filter(symbol_mouse %in% ambiguous_mouse_genes) %>% filter(symbol==toupper(symbol_mouse))
  geneinfo_human = geneinfo_human %>% filter(symbol_mouse %in% unambiguous_mouse_genes) %>% bind_rows(geneinfo_ambiguous_solved) %>% drop_na()

  mousesymbol2humansymbol = mapper(geneinfo_human,"symbol","symbol_mouse")

  converted_symbols = symbols %>% mousesymbol2humansymbol[.]

  return(converted_symbols)
}

get_design = function(E) { # make design matrix for differential expression between celltypes
  TS <- phenoData(E)$celltype
  TS <- factor(TS, levels=unique(TS))
  design <- model.matrix(~0+TS)
  colnames(design) <- levels(TS)
  design
}
cal_diffexp = function(E, c1, c2, design=NULL, v=NULL) {
  celltypemapper = setNames(phenoData(E)$celltype %>% unique %>% make.names, phenoData(E)$celltype %>% unique)
  c1 = celltypemapper[c1]
  c2 = celltypemapper[c2]

  phenoData(E)$celltype = celltypemapper[phenoData(E)$celltype]

  if (is.null(design)) design = get_design(E)
  if (is.null(v)) v = E

  fit = lmFit(v,design)#fit linear model for each gene - cf Limma

  contrast = paste0(c1, "-", c2)
  cont.matrix =eval(parse(text=paste("makeContrasts(", contrast, ",levels=design)",sep="")))
  fit2 = contrasts.fit(fit, cont.matrix)
  #fit2 = contrasts.fit(fit, cont.matrix)
  #?why twice?-> error

  fit.eb = eBayes(fit2)#based on linear model fit and contrast, compute test statistics, DE values, etc...
  options(digits=3)
  toptable=topTable(fit.eb, adjust="BH", sort.by="none",number=Inf)# get a table of DEvalues
  toptable %>% tibble::rownames_to_column("gene") %>% rename(lfc=logFC, qval=adj.P.Val) %>% dplyr::select(lfc, qval, gene)
}
make_diffexp_timeseries = function(toptable) {
  # get a table of DEvalues
  toptable %>% rename(lfc=logFC, qval=adj.P.Val) %>% dplyr::select(lfc, qval, gene)
}





SPL_wrapper = function(ligand,G,id2allgenes,spl_cutoff) {
  # calculate spl distance between ligand and every other node in graph
  distances = igraph::distances(graph = G, v = id2allgenes[ligand], to = igraph::V(G),mode = "out")
  # the shorter the better
  # change na and infinite predictions to -1 in order to later replace them by the maximum (= no probability)
  distances[is.nan(distances)] = -1
  distances[is.infinite(distances)] = -1

  distances[distances == -1] = max(distances)

  # reverse distances to weights: short distance: high weight
  spl_matrix = max(distances) - distances

  if (spl_cutoff > 0){
    spl_matrix_TRUE = apply(spl_matrix,1,function(x){x <= quantile(x,spl_cutoff)}) %>% t()
    spl_matrix[spl_matrix_TRUE] = 0
  }

  spl_vector = apply(spl_matrix, 2, mean)
  return(spl_vector)
}

PPR_wrapper = function(ligand,E,G,delta,id2allgenes,ppr_cutoff) {
  partial_matrix = lapply(ligand,single_ligand_ppr_wrapper,E,G,delta,id2allgenes)
  ppr_matrix = matrix(unlist(partial_matrix), ncol = length(E), byrow = TRUE)
  if (delta == 0 & is.null(ppr_cutoff)){
    ppr_cutoff = 0
  }
  # put cutoff
  if (ppr_cutoff > 0){
    ppr_matrix_TRUE = apply(ppr_matrix,1,function(x){x <= quantile(x,ppr_cutoff)}) %>% t()
    ppr_matrix[ppr_matrix_TRUE] = 0
  }
  # if multiple ligands: total ligand-tf score = maximum score a particular ligand-tf interaction; if only one ligand: just the normal score
  return(apply(ppr_matrix, 2, mean))
}
single_ligand_ppr_wrapper = function(l,E,G,delta,id2allgenes){
  # set ligand as seed in preference vector E
  E[id2allgenes[l]] = 1
  return(igraph::page_rank(G, algo = c("prpack"), vids = igraph::V(G),directed = TRUE, damping = delta, personalized = E) %>% .$vector)
}
direct_wrapper = function(ligand,G,id2allgenes,cutoff) {
  ltf_matrix = G[id2allgenes[ligand],]
  if (length(ligand) == 1){
    ltf_matrix = matrix(ltf_matrix, nrow = 1)
  }
  if(is.null(cutoff)){
    cutoff = 0
  }
  if (cutoff > 0){
    ltf_matrix_TRUE = apply(ltf_matrix,1,function(x){x <= quantile(x,cutoff)}) %>% t()
    ltf_matrix[ltf_matrix_TRUE] = 0
  }
  # if multiple ligands: total ligand-tf score = maximum score a particular ligand-tf interaction; if only one ligand: just the normal score
  return(apply(ltf_matrix, 2, mean))
}
get_pagerank_target = function(weighted_networks, secondary_targets = FALSE) {

  # input check
  if (!is.list(weighted_networks))
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$lr_sig))
    stop("lr_sig must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$gr))
    stop("gr must be a data frame or tibble object")

  requireNamespace("dplyr")

  # load in weighted networks
  ligand_signaling_network = weighted_networks$lr_sig
  regulatory_network = weighted_networks$gr

  # convert ids to numeric for making Matrix::sparseMatrix later on
  allgenes = c(ligand_signaling_network$from, ligand_signaling_network$to, regulatory_network$from, regulatory_network$to) %>% unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% as_tibble()
  id2allgenes = mapper(allgenes_id_tbl,"allgenes_integer","allgenes")

  ligand_signaling_network = ligand_signaling_network %>% mutate(from_allgenes = id2allgenes[from], to_allgenes = id2allgenes[to]) %>% arrange(from_allgenes) %>% dplyr::select(from_allgenes,to_allgenes,weight)

  # Make Matrix::sparse signaling weighted matrix and graph to apply personalized pagerank
  ligand_signaling_network_matrix = Matrix::sparseMatrix(ligand_signaling_network$from_allgenes %>% as.integer, ligand_signaling_network$to_allgenes %>% as.integer, x=ligand_signaling_network$weight %>% as.numeric, dims = c(length(allgenes), length(allgenes)))
  signaling_igraph = igraph::graph_from_adjacency_matrix(ligand_signaling_network_matrix, weighted=TRUE, mode="directed")


  # background pr
  background_pr = igraph::page_rank(signaling_igraph, algo = c("prpack"), vids = igraph::V(signaling_igraph),directed = TRUE)
  background_pr = background_pr$vector

  ltf_matrix = matrix(background_pr, ncol = length(igraph::V(signaling_igraph)), byrow = TRUE)
  colnames(ltf_matrix) = allgenes

  # preparing the gene regulatory matrix
  grn_matrix = construct_tf_target_matrix(weighted_networks)

  # Multiply ligand-tf matrix with tf-target matrix
  ligand_to_target = (ltf_matrix %*% grn_matrix)

  # Secondary targets
  if (secondary_targets == TRUE) {
    ltf_matrix = ligand_to_target

    ligand_to_target_primary = ligand_to_target
    ligand_to_target_secondary = ltf_matrix_%*% grn_matrix

    # set 0's to number higher than 0 to avoid Inf when inverting in sum
    ligand_to_target_primary[ligand_to_target_primary == 0] = Inf
    ligand_to_target_primary[is.infinite(ligand_to_target_primary)] = min(ligand_to_target_primary)

    ligand_to_target_secondary[ligand_to_target_secondary == 0] = Inf
    ligand_to_target_secondary[is.infinite(ligand_to_target_secondary)] = min(ligand_to_target_secondary)

    # inverting in sum to emphasize primary targets more (scores secondary targets matrix tend to be higher )
    ligand_to_target = (ligand_to_target_primary **-1 + ligand_to_target_secondary **-1) ** -1
  }

  ligand_to_target = ligand_to_target %>% as.numeric()
  names(ligand_to_target) = allgenes


  return(ligand_to_target)
}
get_split_auroc = function(observed, known) {
  prediction = ROCR::prediction(observed, known)
  performance = ROCR::performance(prediction, measure="spec", x.measure="rec")
  metrics = tibble(tpr = performance@x.values[[1]], spec = performance@y.values[[1]])
  meetingpoint = which(-(1-metrics$spec) + 1 <= metrics$tpr)[[1]] # < or <= ?
  specs = c((1-metrics$spec)[seq_len(meetingpoint)-1],1-metrics$tpr[meetingpoint], 1)
  recs = c(metrics$tpr[seq_len(meetingpoint)], 0)
  auroc_specificity = caTools::trapz(specs, recs)
  auroc = caTools::trapz(1-metrics$spec, metrics$tpr)
  auroc_sensitivity = auroc - auroc_specificity
  tibble(auroc=auroc, auroc_sensitivity=auroc_sensitivity, auroc_specificity=auroc_specificity)
}
evaluate_target_prediction_strict = function(response,prediction,continuous = TRUE, prediction_response_df = FALSE){
  response_df = tibble(gene = names(response), response = response)
  prediction_df = tibble(gene = names(prediction), prediction = prediction)
  combined = inner_join(response_df,prediction_df, by = "gene")
  if (nrow(combined) == 0)
    stop("Gene names in response don't accord to gene names in ligand-target matrix (did you consider differences human-mouse namings?)")
  prediction_vector = combined$prediction
  names(prediction_vector) = combined$gene
  response_vector = combined$response
  names(response_vector) = combined$gene
  if (continuous == TRUE){
    performance = classification_evaluation_continuous_pred(prediction_vector,response_vector)

  } else{
    performance = classification_evaluation_categorical_pred(prediction_vector,response_vector)
  }
  if (prediction_response_df == TRUE){
    output = list(performance = performance, prediction_response_df = combined)
    return(output)
  } else {
    return(performance)
  }

}
classification_evaluation_continuous_pred = function(prediction,response, iregulon = TRUE){

  if ((sd(response) == 0 & sd(prediction) == 0) | is.null(prediction) | is.null(response)){ # problems can occur otherwise in these leave-one-in models
    return(dplyr::tibble(auroc = NA,
                         aupr = NA,
                         aupr_corrected = NA,
                         sensitivity_roc = NA,
                         specificity_roc = NA,
                         mean_rank_GST_log_pval = NA,
                         auc_iregulon = NA,
                         auc_iregulon_corrected = NA,
                         pearson = NA,
                         spearman = NA))
  }
  prediction_ROCR = ROCR::prediction(prediction, response)
  performance1 = ROCR::performance(prediction_ROCR, measure="tpr", x.measure="fpr")

  performance2 = ROCR::performance(prediction_ROCR, measure="prec", x.measure="rec")
  performance = tibble(tpr = performance1@y.values[[1]], fpr=performance1@x.values[[1]], precision=performance2@y.values[[1]], recall=performance2@x.values[[1]])

  performance = performance %>% replace_na(list(recall=0, precision=1))

  aupr = caTools::trapz(performance$recall, performance$precision)
  pos_class = sum(response)
  total = length(response)
  aupr_random = pos_class/total


  metrics = get_split_auroc(prediction, response)

  sensitivity = metrics$auroc_sensitivity
  specificity = metrics$auroc_specificity
  auroc = metrics$auroc

  cor_p = cor(prediction, response)
  cor_s = cor(prediction, response, method = "s")

  cor_p_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response))) %>% .$p.value
  cor_s_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response), method =  "s")) %>% .$p.value

  mean_rank_GST = limma::wilcoxGST(response, prediction)
  #### now start calculating the AUC-iRegulon
  tbl_perf = tibble(auroc = auroc,
                    aupr = aupr,
                    aupr_corrected = aupr - aupr_random,
                    sensitivity_roc = sensitivity,
                    specificity_roc = specificity,
                    mean_rank_GST_log_pval = -log(mean_rank_GST),
                    pearson_log_pval = -log10(cor_p_pval),
                    spearman_log_pval = -log10(cor_s_pval),
                    pearson = cor_p,
                    spearman = cor_s)
  if (iregulon == TRUE){
    output_iregulon = calculate_auc_iregulon(prediction,response)
    tbl_perf = tibble(auroc = auroc,
                      aupr = aupr,
                      aupr_corrected = aupr - aupr_random,
                      sensitivity_roc = sensitivity,
                      specificity_roc = specificity,
                      mean_rank_GST_log_pval = -log(mean_rank_GST),
                      auc_iregulon = output_iregulon$auc_iregulon,
                      auc_iregulon_corrected = output_iregulon$auc_iregulon_corrected,
                      pearson_log_pval = -log10(cor_p_pval),
                      spearman_log_pval = -log10(cor_s_pval),
                      pearson = cor_p,
                      spearman = cor_s)
  }
  return(tbl_perf)
}
classification_evaluation_categorical_pred = function(predictions, response) {
  # print(head(predictions))
  # print(length(predictions))

  if (sd(response) == 0){ # if all response is the same
    return(dplyr::tibble(accuracy = NA,
                         recall = NA,
                         specificity = NA,
                         precision = NA,
                         F1 =  NA,
                         F05 = NA,
                         F2 = NA,
                         mcc = NA,
                         informedness = NA,
                         markedness = NA,
                         fisher_pval_log = NA,
                         fisher_odds = NA))
  }

  num_positives = sum(response)
  num_total = length(response)

  # calculate base statistics
  pos_preds = sum(predictions)

  tp = sum(response[predictions])
  fp = pos_preds - tp

  num_negatives = num_total - num_positives

  tp = tp
  fp = fp
  fn = num_positives - tp
  tn = num_negatives - fp
  npv = tn / (tn + fn)
  if (sd(predictions) == 0){
    fisher = list(p.value = NA, estimate = NA)
  } else {
    fisher = fisher.test(as.factor(response), predictions)
  }

  mcc_S = (tp + fn)/num_total
  mcc_P = (tp + fp)/num_total

  metrics = dplyr::tibble(
    accuracy = (tp + tn) / (num_positives + num_negatives),
    recall = tp / num_positives,
    specificity = tn / num_negatives,
    precision = tp / (tp + fp),
    F1 =  (2 * precision * recall) / (precision + recall),
    F05 = (1.25 * precision * recall) / (0.25 * precision + recall),
    F2 = (5 * precision * recall) / (4 * precision + recall),
    mcc = (tp/num_total - mcc_S * mcc_P)/sqrt(mcc_P * mcc_S * (1-mcc_S) * (1-mcc_P)) ,
    informedness = recall + specificity - 1,
    markedness = precision + npv - 1,
    fisher_pval_log = -log(fisher$p.value),
    fisher_odds = fisher$estimate
  )
  if (sd(predictions) == 0){ # all predictions are the same!
    return(dplyr::tibble(accuracy = metrics$accuracy,
                         recall = metrics$recall,
                         specificity = metrics$specificity,
                         precision = 0,
                         F1 =  0,
                         F05 = 0,
                         F2 = 0,
                         mcc = 0,
                         informedness = 0,
                         markedness = 0,
                         fisher_pval_log = NA,
                         fisher_odds = NA))
  }
  return(metrics)
}
rank_desc = function(x){rank(desc(x), ties.method = "max")}
# rank_desc = function(x){rank(desc(x), ties.method = "random")}

.calcAUC = function(oneRanking, aucThreshold, maxAUC)
{
  x = unlist(oneRanking)
  x = sort(x[x<aucThreshold])
  y = 1:length(x)
  a = diff(c(x, aucThreshold)) * y
  return(sum(a)/maxAUC)
}

calculate_auc_iregulon = function(prior,response){

  genes_prior = names(prior)
  dim(prior) = c(1,length(prior))
  colnames(prior) = genes_prior
  rownames(prior) = "ligand"

  prior_rank = apply(prior,1,rank_desc)
  rankings = tibble(prior=prior_rank[,1], rn = rownames(prior_rank))

  fake_rankings = rankings %>% mutate(rn = sample(rn))
  rankings = data.table::data.table(rankings)
  fake_rankings = data.table::data.table(fake_rankings)

  aucMaxRank = 0.03*nrow(rankings)

  # calculate enrichment over the expression settings

  geneSet = response[response == TRUE] %>% names()
  geneSet = unique(geneSet)
  nGenes = length(geneSet)
  geneSet = geneSet[which(geneSet %in% rankings$rn)]

  missing = nGenes-length(geneSet)

  gSetRanks = subset(rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed

  aucThreshold = round(aucMaxRank)
  maxAUC = aucThreshold * nrow(gSetRanks)

  auc_iregulon = sapply(gSetRanks, .calcAUC, aucThreshold, maxAUC)

  gSetRanks_fake = subset(fake_rankings, rn %in% geneSet)[,-"rn", with=FALSE] # gene names are no longer needed
  auc_iregulon_fake = sapply(gSetRanks_fake, .calcAUC, aucThreshold, maxAUC)

  auc_iregulon_corrected = auc_iregulon - auc_iregulon_fake
  return(list(auc_iregulon = auc_iregulon, auc_iregulon_corrected = auc_iregulon_corrected))
}

evaluate_target_prediction_bins_direct = function(bin_id,settings,ligand_target_matrix,ligands_position,ncitations){
  requireNamespace("dplyr")
  targets_bin_id = ncitations %>% filter(bin == bin_id) %>% .$symbol %>% unique()
  if (ligands_position == "cols"){
    ligand_target_matrix = ligand_target_matrix[targets_bin_id,]
  } else if (ligands_position == "rows") {
    ligand_target_matrix = ligand_target_matrix[,targets_bin_id]
  }
  performances = bind_rows(lapply(settings, evaluate_target_prediction, ligand_target_matrix,ligands_position)) %>% mutate(target_bin_id = bin_id)
}
caret_classification_evaluation_continuous = function(data, lev = NULL, model = NULL){
  # print(data)
  # print(dim(data))
  if (length(unique(data$obs)) != 2)
    stop(paste("Your outcome has no 2 different classes as required here"))
  prediction = data[,4]
  response = data$obs %>% as.character() %>% gsub("\\.","",.) %>% as.logical()
  out_tibble = classification_evaluation_continuous_pred(prediction, response, iregulon = FALSE)
  # print(out_tibble)
  out = out_tibble[1,] %>% as.numeric()
  names(out) = colnames(out_tibble)
  out
}
caret_classification_evaluation_categorical = function(data, lev = NULL, model = NULL){
  # print(data)
  # print(dim(data))
  if (length(unique(data$obs)) != 2)
    stop(paste("Your outcome has no 2 different classes as required here"))
  prediction = data$pred %>% as.character() %>% gsub("\\.","",.) %>% as.logical()
  response = data$obs %>% as.character() %>% gsub("\\.","",.) %>% as.logical()
  # print(sum(prediction))
  # print(length(response))
  out_tibble = classification_evaluation_categorical_pred(prediction, response)
  # print(out_tibble)
  out = out_tibble[1,] %>% as.numeric()
  names(out) = colnames(out_tibble)
  out
}
wrapper_caret_classification = function(train_data, algorithm, continuous = TRUE, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4,prediction_response_df = NULL,ignore_errors = FALSE, return_model = FALSE){


  if (sum( table(train_data$obs) >= cv_number) != 2 )
    stop("Make sure that there are more instances of each class than the folds in the cross-validation scheme")
  if (parallel == TRUE){
    requireNamespace("doMC", quietly = TRUE)
    doMC::registerDoMC(cores = n_cores)
  }

  if (continuous == TRUE){
    if (cv == TRUE){
      control =  caret::trainControl(method="repeatedcv",
                                     number=cv_number,
                                     repeats=cv_repeats,
                                     preProcOptions = NULL,
                                     classProbs = TRUE,
                                     summaryFunction = caret_classification_evaluation_continuous)
    } else if (cv == FALSE) {
      control =  caret::trainControl(method="none",
                                     preProcOptions = NULL,
                                     classProbs = TRUE,
                                     summaryFunction = caret_classification_evaluation_continuous)
    }
    if (ignore_errors == TRUE){
      # avoid errors due to bad splits during cross-validation that make that not both classes are present
      caret_train = purrr::safely(caret::train)
      result = NULL
      while(is.null(result)){
        model = caret_train(y = train_data$obs,
                            x = train_data[,-(which(colnames(train_data) == "obs"))],
                            method=algorithm,
                            trControl=control,
                            metric = 'auroc')
        result = model$result
        # if(!is.null(model$error)){print(model$error)}

      }
      model = model$result
    } else {

      model = caret::train(y = train_data$obs,
                           x = train_data[,-(which(colnames(train_data) == "obs"))],
                           method=algorithm,
                           trControl=control,
                           metric = 'auroc')
    }

  } else if (continuous == FALSE) {
    # print("here")

    if( cv == TRUE){
      control =  caret::trainControl(method="repeatedcv",
                                     number=cv_number,
                                     repeats=cv_repeats,
                                     preProcOptions = NULL,
                                     classProbs = FALSE,
                                     summaryFunction = caret_classification_evaluation_categorical)
    } else if (cv == FALSE) {
      control =  caret::trainControl(method="none",
                                     preProcOptions = NULL,
                                     classProbs = FALSE,
                                     summaryFunction = caret_classification_evaluation_categorical)
    }
    if (ignore_errors == TRUE){
      # avoid errors due to bad splits during cross-validation that make that not both classes are present
      caret_train = purrr::safely(caret::train)
      result = NULL
      while(is.null(result)){
        model = caret_train(y = train_data$obs,
                            x = train_data[,-(which(colnames(train_data) == "obs"))],
                            method=algorithm,
                            trControl=control,
                            metric = 'mcc') #mcc
        result = model$result
      }
      model = model$result
    } else {
      model = caret::train(y = train_data$obs,
                           x = train_data[,-(which(colnames(train_data) == "obs"))],
                           method=algorithm,
                           trControl=control,
                           metric = 'mcc') #mcc
    }
    # if(!is.null(model$error)){print(model$error)}


  }

  performances =  model$resample
  if(continuous == TRUE){
    final_model_predictions = predict(model,newdata = train_data, type = "prob") %>% .[,2]
    response_class = train_data$obs %>% as.character() %>% gsub("\\.","",.) %>% as.logical()
    performances_training_continuous = classification_evaluation_continuous_pred(prediction = final_model_predictions, response = response_class,iregulon = FALSE)
    performances_training = classification_evaluation_categorical_pred(predictions = final_model_predictions >= 0.5, response = response_class)
  } else {
    final_model_predictions = predict(model,newdata = train_data)  %>% gsub("\\.","",.) %>% as.logical()
    response_class = train_data$obs %>% as.character() %>% gsub("\\.","",.) %>% as.logical()
    performances_training = classification_evaluation_categorical_pred(predictions = final_model_predictions, response = response_class)
    performances_training_continuous = NULL
  }

  if (var_imps == TRUE) {
    imps =  caret::varImp(model, scale = FALSE)
    var_imps_df = imps$importance  %>% tibble::rownames_to_column("feature") %>% as_tibble() %>% .[,1:2]
    colnames(var_imps_df) = c("feature","importance")
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, performances_training_continuous = performances_training_continuous ,var_imps = var_imps_df, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    } else {
    output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, performances_training_continuous = performances_training_continuous, var_imps = var_imps_df)}

  } else {
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, performances_training_continuous = performances_training_continuous, var_imps = NULL, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    } else {
    output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training,performances_training_continuous = performances_training_continuous, var_imps = NULL)}
  }

  if (return_model == TRUE){
    output_list$model = model
  }
  return(output_list)

}
# validation
make_new_setting_ligand_prediction_multi_validation = function(setting,all_ligands){

  if (length(setting$from) > 1) {
    setting$from = paste0(setting$from,collapse = "-")
  }

  new_setting = list()
  new_setting$name = setting$name
  new_setting$ligand = setting$from
  new_setting$from = all_ligands
  new_setting$response = setting$response
  return(new_setting)
}
make_new_setting_ligand_prediction_single_validation = function(setting,test_ligand){

  if (length(setting$from) > 1) {
    setting$from = paste0(setting$from,collapse = "-")
  }

  new_setting = list()
  new_setting$name = setting$name
  new_setting$ligand = setting$from
  new_setting$from = test_ligand
  new_setting$response = setting$response
  return(new_setting)
}


# application
make_new_setting_ligand_prediction_multi_application = function(setting,all_ligands){
  new_setting = list()
  new_setting$name = setting$name
  new_setting$from = all_ligands
  new_setting$response = setting$response
  return(new_setting)
}
make_new_setting_ligand_prediction_single_application = function(setting,test_ligand){
  new_setting = list()
  new_setting$name = setting$name
  new_setting$from = test_ligand
  new_setting$response = setting$response
  return(new_setting)
}

filter_genes_ligand_target_matrix = function(ligand_target_matrix, ligands_position = cols){
  if (ligands_position == "cols"){
    target_genes = rownames(ligand_target_matrix)
    sd_genes = apply(ligand_target_matrix,1,sd)
    ligand_target_matrix_ = ligand_target_matrix[sd_genes > quantile(sd_genes,0.5),]
  } else if (ligands_position == "rows") {
    target_genes = colnames(ligand_target_matrix)
    sd_genes = apply(ligand_target_matrix,2,sd)
    ligand_target_matrix_ = ligand_target_matrix[,sd_genes > quantile(sd_genes,0.5)]
  }
  return(ligand_target_matrix_)
}

is_ligand_active = function(importances){
  if(nrow(importances) == 0){
    return(NULL)
  }

  test_ligand = importances$test_ligand
  real_ligand = strsplit(importances$ligand,"[-]")
  true_ligand = rep(NULL, times = length(test_ligand))
  for (i in seq(length(test_ligand))){ #length ligand instead test_ligand??
    test = test_ligand[i]
    real = real_ligand[[i]]
    true_ligand[i] = test %in% real
  }
  return(true_ligand)
}
are_ligands_oi_active = function(importances, ligands_oi){
  test_ligand = importances$test_ligand
  real_ligand = strsplit(importances$ligand,"[-]")
  true_ligand = rep(NULL, times = length(test_ligand))
  for (i in seq(length(test_ligand))){
    test = test_ligand[i]
    real = real_ligand[[i]]
    true_ligand[i] = sum(ligands_oi %in% real) > 0
  }
  return(true_ligand)
}
evaluate_target_prediction_regression_strict = function(response,prediction, prediction_response_df = FALSE){
  response_df = tibble(gene = names(response), response = response)
  prediction_df = tibble(gene = names(prediction), prediction = prediction)
  combined = inner_join(response_df,prediction_df, by = "gene")
  if (nrow(combined) == 0)
    stop("Gene names in response don't accord to gene names in ligand-target matrix (did you consider differences human-mouse namings?)")
  prediction_vector = combined$prediction
  names(prediction_vector) = combined$gene
  response_vector = combined$response
  names(response_vector) = combined$gene

  performance = regression_evaluation(prediction_vector,response_vector)

  if (prediction_response_df == TRUE){
    output = list(performance = performance, prediction_response_df = combined)
    return(output)
  } else {
    return(performance)
  }
}
regression_evaluation = function(prediction,response){
  model = lm(response~prediction)
  model_summary = summary(model)

  cor_p = cor(prediction, response)
  cor_s = cor(prediction, response, method = "s")

  cor_p_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response))) %>% .$p.value
  cor_s_pval = suppressWarnings(cor.test(as.numeric(prediction), as.numeric(response), method =  "s")) %>% .$p.value

  tbl_perf = tibble(r_squared = model_summary$r.squared,
                    adj_r_squared = model_summary$adj.r.squared,
                    f_statistic = model_summary$fstatistic[1],
                    lm_coefficient_abs_t = model_summary$coefficients %>% .[2,3] %>% abs(),
                    inverse_rmse = 1/(model_summary$sigma),
                    reverse_aic = AIC(model) * -1,
                    reverse_bic = BIC(model) * -1,
                    inverse_mae = 1/(mae(model$residuals)),
                    pearson_log_pval = -log10(cor_p_pval),
                    spearman_log_pval = -log10(cor_s_pval),
                    pearson_regression = cor_p,
                    spearman_regression = cor_s)

  return(tbl_perf)
}

# rmse = function(error){
#   sqrt(mean(error^2))
# }
mae = function(error){
  mean(abs(error))
}
wrapper_caret_regression = function(train_data, algorithm, var_imps = TRUE, cv = TRUE, cv_number = 4, cv_repeats = 2, parallel = FALSE, n_cores = 4,prediction_response_df = NULL,ignore_errors = FALSE, return_model = FALSE){

  if (parallel == TRUE){
    requireNamespace("doMC", quietly = TRUE)
    doMC::registerDoMC(cores = n_cores)
  }

  if (cv == TRUE){
    control =  caret::trainControl(method="repeatedcv",
                                   number=cv_number,
                                   repeats=cv_repeats,
                                   preProcOptions = NULL,
                                   summaryFunction = caret_regression_evaluation_continuous)
  } else if (cv == FALSE) {
    control =  caret::trainControl(method="none",
                                   preProcOptions = NULL,
                                   summaryFunction = caret_regression_evaluation_continuous)
  }
  if (ignore_errors == TRUE){
    # avoid errors due to bad splits during cross-validation that make that not both classes are present
    caret_train = purrr::safely(caret::train)
    result = NULL
    while(is.null(result)){
      model = caret_train(y = train_data$obs,
                          x = train_data[,-(which(colnames(train_data) == "obs"))],
                          method=algorithm,
                          trControl=control,
                          metric = 'inverse_rmse',
                          maximize = TRUE,
                          importance = TRUE) # RMSE should be minimized - so maximize its inverse
      result = model$result
      # if(!is.null(model$error)){print(model$error)}

    }
    model = model$result
  } else {
    model = caret::train(y = train_data$obs,
                         x = train_data[,-(which(colnames(train_data) == "obs"))],
                         method=algorithm,
                         trControl=control,
                         metric = 'inverse_rmse',
                         maximize = TRUE,
                         importance = TRUE) # RMSE should be minimized - so maximize its inverse
  }


  performances =  model$resample

  final_model_predictions = predict(model,newdata = train_data)
  response = train_data$obs
  performances_training = regression_evaluation(prediction = final_model_predictions, response = response)

  if (var_imps == TRUE) {
    imps =  caret::varImp(model, scale = FALSE)
    var_imps_df = imps$importance  %>% tibble::rownames_to_column("feature") %>% as_tibble() %>% .[,1:2]
    colnames(var_imps_df) = c("feature","importance")
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, var_imps = var_imps_df, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    }
    output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, var_imps = var_imps_df)

  } else {
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, var_imps = NULL, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    }
    output_list = list(performances = performances %>% as_tibble(), performances_training = performances_training, var_imps = NULL)
  }

  if (return_model == TRUE){
    output_list$model = model
  }
  return(output_list)

}
caret_regression_evaluation_continuous = function(data, lev = NULL, model = NULL){
  # print(data)
  # print(dim(data))
  prediction = data$pred
  response = data$obs
  out_tibble = regression_evaluation(prediction, response)
  # print(out_tibble)
  out = out_tibble[1,] %>% as.numeric()
  names(out) = colnames(out_tibble)
  out
}
get_shortest_path_signaling = function(ligand_oi, signaling_df, signaling_igraph){
  ligand_signaling = signaling_df %>% filter(ligand == ligand_oi)
  tfs = ligand_signaling$TF %>% unique()
  sp = igraph::shortest_paths(signaling_igraph, from = ligand_oi, to = tfs, mode ="out")
  tf_nodes = unlist(sp$vpath) %>% names() %>% unique() %>% .[. != ligand_oi] # ligand should not belong to tf nodes
}

construct_ligand_signaling_df = function(ligands_all,targets_all,k,weighted_networks,ligand_tf_matrix){
  final_combined_df = bind_rows(expand.grid(ligands_all,targets_all) %>% apply(.,1,wrappper_visualization,k,weighted_networks,ligand_tf_matrix))
}
wrappper_visualization = function(grid,k,weighted_networks,ligand_tf_matrix){
  ligand = grid[1]
  target = grid[2]
  network = get_network_df(ligand,target,k,weighted_networks,ligand_tf_matrix)
}

get_network_df = function(ligand_to_vis,target_to_vis,k,weighted_networks,ligand_tf_matrix){

  ## prepare TFs downstream of ligands

  ## ligands should be in columns

  ligand_tf_matrix_visualization = ligand_tf_matrix[,ligand_to_vis] %>% as.matrix()

  colnames(ligand_tf_matrix_visualization) = "V1"

  ligand_tf_matrix_visualization_df = as_tibble(ligand_tf_matrix_visualization) %>% rename(weight = V1) %>% mutate(TF = rownames(ligand_tf_matrix)) %>% select(TF,weight)
  ligand_tf_matrix_visualization_df_filtered = ligand_tf_matrix_visualization_df %>% filter(weight > 0) %>% mutate(ligand = ligand_to_vis)

  ## prepare TFs upstream of target
  regulatory_network_filtered = weighted_networks$gr %>% filter(to == target_to_vis) %>% rename(TF = from, weight_grn = weight)

  ## combine both
  combined_df = inner_join(ligand_tf_matrix_visualization_df_filtered,regulatory_network_filtered, by = "TF")
  combined_df = combined_df %>% mutate(total_weight = weight*weight_grn)
  final_combined_df = combined_df %>% arrange(-total_weight) %>% .[1:min(k,nrow(combined_df)),]

  return(final_combined_df)
}
cal_celltype_average_wrapper = function(E){
  avexprs = lapply(levels(E$celltype), cal_celltype_average, E)
  names(avexprs) = levels(E$celltype)
  avexprs = avexprs %>% bind_cols() %>% as.matrix()
  rownames(avexprs) = rownames(exprs(E))
  return(avexprs)
}
cal_celltype_average = function(cell,E){
  E = E[, E$celltype == cell]
  expression = exprs(E)
  average_expression = apply(expression, 1, mean)

}
evaluate_ligand_prediction_bins_direct = function(bin_id,settings,ligand_target_matrix,ligands_position,ncitations,...){
  requireNamespace("dplyr")
  targets_bin_id = ncitations %>% filter(bin == bin_id) %>% .$symbol %>% unique()

  if (ligands_position == "cols"){
    ligand_target_matrix = ligand_target_matrix[targets_bin_id,]
  } else if (ligands_position == "rows") {
    ligand_target_matrix = ligand_target_matrix[,targets_bin_id]
    ligand_target_matrix = ligand_target_matrix %>% t()
  }
  ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix(...)
  ligand_target_matrix_discrete = ligand_target_matrix %>% make_discrete_ligand_target_matrix()


  # performances = bind_rows(lapply(settings, evaluate_target_prediction, ligand_target_matrix,ligands_position)) ##### checking!!
  all_ligands = unlist(extract_ligands_from_settings(settings, combination = FALSE))
  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = TRUE)

  # print(all_ligands)
  # print(colnames(ligand_target_matrix_discrete))
  # print(colnames(ligand_target_matrix))

  ligand_importances = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix[, all_ligands]))
  ligand_importances_discrete = bind_rows(lapply(settings_ligand_pred, get_single_ligand_importances, ligand_target_matrix_discrete[, all_ligands]))
  # ligand_importances_discrete = ligand_importances_discrete %>% select_if(.predicate = function(x){sum(is.na(x)) == 0})
  # if(sum(is.na(ligand_importances_discrete$fisher_odds)) > 0){
  #   ligand_importances_discrete = ligand_importances_discrete %>% select(-fisher_odds) %>% select(-fisher_pval_log) # because contains too much NA sometimes in leave one in models
  # }

  settings_ligand_pred = convert_settings_ligand_prediction(settings, all_ligands, validation = TRUE, single = FALSE)
  ligand_importances_glm = bind_rows(lapply(settings_ligand_pred, get_multi_ligand_importances, ligand_target_matrix[,all_ligands], algorithm = "glm", cv = FALSE)) %>% rename(glm_imp = importance)

  all_importances = full_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand")) %>% full_join(ligand_importances_discrete, by = c("setting","test_ligand", "ligand"))
  # all_importances = inner_join(ligand_importances, ligand_importances_glm, by = c("setting","test_ligand","ligand"))

  all_importances = all_importances %>% select_if(.predicate = function(x){sum(is.na(x)) == 0}) # importance measures full of NAs possible dependent on input networks - remove these columns for further analysis

  # evaluation = suppressWarnings(evaluate_importances_ligand_prediction(all_importances, "median","lda",cv_number = 3, cv_repeats = 20))
  # performances_ligand_prediction = evaluation$performances %>% select(-Resample) %>% mutate_all(mean) %>% distinct()
  # # warning lda here: variables are collinear --> not problematic but logical here

  performances_ligand_prediction_single = evaluate_single_importances_ligand_prediction(all_importances, "median")

  performances = performances_ligand_prediction_single %>% filter(auroc == max(auroc)) %>% mutate(target_bin_id = bin_id)
  return(performances)
}

get_affected_targets_output = function(diff_matrix, rankings){
  return(lapply(colnames(diff_matrix),function(ligand_oi,diff_matrix, rankings){
    vector_oi = diff_matrix[,ligand_oi]
    if(rankings == TRUE){
      sorted_vector_oi_high2low = sort(vector_oi)
      sorted_vector_oi_low2high = sort(vector_oi, decreasing = TRUE)
    } else {
      sorted_vector_oi_high2low = sort(vector_oi, decreasing = TRUE)
      sorted_vector_oi_low2high = sort(vector_oi)
    }

    return(list(ligand = ligand_oi, targets_higher = sorted_vector_oi_high2low %>% head(100), targets_lower = sorted_vector_oi_low2high %>% head(100)))
  },diff_matrix,rankings ))
}
wrapper_evaluate_multi_ligand_target_prediction = function(...){
  average_performances = evaluate_multi_ligand_target_prediction(...) %>% .$performances %>% select(-Resample) %>% summarise_all(mean)
}
subtract_performances = function(performances, metric){
  performances_base = performances %>% filter(model_name == "complete_model")
  performances_oi = performances %>% filter(model_name != "complete_model")
  performances_diff = (performances_oi %>% select(metric) %>% unlist()) - (performances_base %>% select(metric) %>% unlist())
  performances_diff_df = tibble(metric = performances_diff, model_name = performances_oi$model_name)
  colnames(performances_diff_df)  = c(paste0(metric), "model_name")
  return(performances_diff_df )
}
combine_loi_loo_performances = function(loi_performances,loo_performances){
  loi_diff_auroc = subtract_performances(loi_performances, "auroc") %>% rename(auroc_loi = auroc)
  loi_diff_aupr = subtract_performances(loi_performances, "aupr") %>% rename(aupr_loi = aupr)
  loi_diff_pearson = subtract_performances(loi_performances, "pearson") %>% rename(pearson_loi = pearson)
  loo_diff_auroc = subtract_performances(loo_performances, "auroc") %>% rename(auroc_loo = auroc)
  loo_diff_aupr = subtract_performances(loo_performances, "aupr") %>% rename(aupr_loo = aupr)
  loo_diff_pearson = subtract_performances(loo_performances, "pearson") %>% rename(pearson_loo = pearson)
  input_df = purrr::reduce(list(loi_diff_aupr, loi_diff_auroc, loi_diff_pearson, loo_diff_aupr, loo_diff_auroc, loo_diff_pearson), inner_join, by = "model_name") %>% select(model_name, aupr_loi, auroc_loi:pearson_loo)
}
regression_characterization_optimization = function(loi_performances, loo_performances,source_weights_df, random_forest = FALSE){
  input_df = combine_loi_loo_performances(loi_performances,loo_performances)
  combined_df = inner_join(input_df %>% rename(source = model_name), source_weights_df, by = "source") %>% select(-source)
  if (random_forest == FALSE){
    model = lm(weight ~ ., data = combined_df)
  } else {
    model = randomForest::randomForest(data = combined_df %>% drop_na(), weight ~ ., ntree = 1000)
  }
  return(model)
}
assign_new_weight = function(loi_performances, loo_performances, output_regression_model, source_weights_df){
  performances_oi = combine_loi_loo_performances(loi_performances, loo_performances)
  source_oi = performances_oi$model_name
  performances_oi = performances_oi %>% select(-model_name)
  weight_oi =  predict(output_regression_model,performances_oi)
  weight_oi = min(weight_oi,1) # weight should be lower than 1
  weight_oi = max(weight_oi,0) # weight should be higher than 0
  source_weights_df = source_weights_df %>% bind_rows(tibble(source = source_oi, weight = weight_oi ))
}
train_rf = function(setting,ligand_target_matrix, ligands_position = "cols", ntrees = 1000, mtry = 2,  continuous = TRUE, known = TRUE, filter_genes = FALSE){
  setting_name = setting$name
  ligands_oi = setting$from

  prediction_matrix = ligand_target_matrix[,ligands_oi]
  target_genes = rownames(ligand_target_matrix)

  response_vector = setting$response
  response_df = tibble(gene = names(response_vector), response = response_vector %>% make.names() %>% as.factor())

  prediction_df = prediction_matrix %>% data.frame() %>% as_tibble()

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  train_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()
  #
  # w_pos = table(train_data$obs)["TRUE."]/length(train_data$obs)
  # w_neg = table(train_data$obs)["FALSE."]/length(train_data$obs)
  #
  # classwt = c(w_neg,w_pos) %>% set_names(c("FALSE.","TRUE."))
  # set.seed(1)
  # rf_model = randomForest::randomForest(y = train_data$obs,
  #                                       x = train_data[,-(which(colnames(train_data) == "obs"))],
  #                                       ntree = ntrees,
  #                                       mtry = ncol(train_data[,-(which(colnames(train_data) == "obs"))])**(1/mtry) %>% ceiling(),
  #                                       importance = FALSE,
  #                                       classwt = classwt
  #
  # )
  set.seed(1)
  rf_model = randomForest::randomForest(y = train_data$obs,
                                        x = train_data[,-(which(colnames(train_data) == "obs"))],
                                        ntree = ntrees,
                                        mtry = ncol(train_data[,-(which(colnames(train_data) == "obs"))])**(1/mtry) %>% ceiling(),
                                        importance = FALSE

  )
}
#' @import e1071
test_rf = function(setting,rf_model, ligand_target_matrix, ligands_position = "cols"){
  requireNamespace("e1071")
  setting_name = setting$name
  ligands_oi = setting$from

  prediction_matrix = ligand_target_matrix[,ligands_oi]
  target_genes = rownames(ligand_target_matrix)

  response_vector = setting$response
  response_df = tibble(gene = names(response_vector), response = response_vector %>% make.names() %>% as.factor())

  prediction_df = prediction_matrix %>% data.frame() %>% as_tibble()

  prediction_df = tibble(gene = target_genes) %>% bind_cols(prediction_df)
  combined = inner_join(response_df,prediction_df, by = "gene")
  test_data = combined %>% select(-gene) %>% rename(obs = response) %>% data.frame()

  rf_prediction = predict(rf_model,test_data[,-(which(colnames(test_data) == "obs"))], type = "prob")

  return(tibble(gene = combined$gene, response = combined$response, prediction = rf_prediction[,"TRUE."]))

  # ModelMetrics::confusionMatrix(predicted = rf_prediction[,"TRUE."],
  #                               actual = test_data$obs,
  #                 cutoff = quantile(rf_prediction[,"TRUE."],0.5))

}
rf_target_prediction = function(i,affected_genes_grouped,strict_background_expressed_genes_grouped,best_upstream_ligands,ligand_target_matrix){
  training_indices = seq(length(affected_genes_grouped)) %>% .[. != i]
  training_affected_genes = affected_genes_grouped[training_indices] %>% unlist() %>% unique()
  training_background_expressed_genes = c(training_affected_genes,strict_background_expressed_genes_grouped[training_indices] %>% unlist() %>% unique()) %>% unique()

  test_affected_genes = affected_genes_grouped[[i]]
  test_background_expressed_genes = c(affected_genes_grouped[[i]],strict_background_expressed_genes_grouped[[i]]) %>% unique()

  setting = lapply(list(training_affected_genes), convert_gene_list_settings_evaluation, name = "target_pred", ligands_oi = best_upstream_ligands, background = training_background_expressed_genes)
  rf_model = setting %>% lapply(train_rf,ligand_target_matrix[,best_upstream_ligands]) %>% .[[1]]

  # Predicting response variable
  setting = lapply(list(test_affected_genes), convert_gene_list_settings_evaluation, name = "target_pred", ligands_oi = best_upstream_ligands, background = test_background_expressed_genes)
  rf_prediction_confusion_matrix = setting %>% lapply(test_rf,rf_model,ligand_target_matrix[,best_upstream_ligands]) %>% .[[1]]

}
apply_quantile_scale = function (x, addend, multiplier)
  # same function as apply_quantile_scale from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the great (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)!
{
  if (is.null(dim(x))) {
    sc <- apply_quantile_scale(matrix(x, ncol = 1), addend = addend,
                               multiplier = multiplier)
    out <- sc[, 1]
    names(out) <- names(x)
    attr(out, "addend") <- attr(sc, "addend")
    attr(out, "multiplier") <- attr(sc, "multiplier")
    out
  }
  else {
    y <- apply_uniform_scale(x, addend, multiplier)
    y[y > 1] <- 1
    y[y < 0] <- 0
    y
  }
}
apply_uniform_scale = function (x, addend, multiplier)
  # same function as apply_uniform_scale from dynutils (copied here for use in vignette to avoid having dynutils as dependency)
  # credits to the great (w/z)outer and r(obrecht)cannood(t) from dynverse (https://github.com/dynverse)!
{
  if (is.null(dim(x))) {
    sc <- apply_uniform_scale(matrix(x, ncol = 1), addend = addend,multiplier = multiplier)
    out <- sc[, 1]
    names(out) <- names(x)
    attr(out, "addend") <- attr(sc, "addend")
    attr(out, "multiplier") <- attr(sc, "multiplier")
    out
  }
  else {
    y <- x %>% sweep(2, addend, "+") %>% sweep(2, multiplier,
                                               "*")
    attr(y, "addend") <- addend
    attr(y, "multiplier") <- multiplier
    y
  }
}

