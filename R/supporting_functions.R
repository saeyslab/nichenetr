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
  allgenes_id_tbl = data.frame(allgenes,allgenes_integer) %>% tbl_df()
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
mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
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
evaluate_target_prediction_strict = function(response,prediction,continuous = TRUE){
  response_df = tibble(gene = names(response), response = response)
  prediction_df = tibble(gene = names(prediction), prediction = prediction)
  combined = inner_join(response_df,prediction_df, by = "gene")

  prediction_vector = combined$prediction
  names(prediction_vector) = combined$gene
  response_vector = combined$response
  names(response_vector) = combined$gene

  if (continuous == TRUE){
    performance = classification_evaluation_continuous_pred(prediction_vector,response_vector)

  } else{
    performance = classification_evaluation_categorical_pred(prediction_vector,response_vector)
  }
  return(performance)
}
classification_evaluation_continuous_pred = function(prediction,response, iregulon = FALSE){

  prediction_ROCR = ROCR::prediction(prediction, response)
  performance1 = ROCR::performance(prediction_ROCR, measure="tpr", x.measure="fpr")

  performance2 = ROCR::performance(prediction_ROCR, measure="prec", x.measure="rec")
  performance = tibble(tpr = performance1@y.values[[1]], fpr=performance1@x.values[[1]], precision=performance2@y.values[[1]], recall=performance2@x.values[[1]])

  performance = performance %>% replace_na(list(recall=0, precision=1))

  aupr = caTools::trapz(performance$recall, performance$precision) # i hope this is correct, but still gotta check
  pos_class = sum(response)
  total = length(response)
  aupr_random = pos_class/total


  metrics = get_split_auroc(prediction, response)
  sensitivity = metrics$auroc_sensitivity
  specificity = metrics$auroc_specificity
  auroc = metrics$auroc

  cor_p = cor(prediction, response)
  cor_s = cor(prediction, response, method = "s")

  mean_rank_GST = limma::wilcoxGST(response, prediction)
  #### now start calculating the AUC-iRegulon
  tbl_perf = tibble(auroc = auroc,
                    aupr = aupr,
                    aupr_corrected = aupr - aupr_random,
                    sensitivity_roc = sensitivity,
                    specificity_roc = specificity,
                    mean_rank_GST_log_pval = -log(mean_rank_GST),
                    pearson = cor_p,
                    spearman = cor_s)
  if (iregulon == TRUE){
    output_iregulon = calculate_auc_iregulon(prediction,response)
    tbl_perf = tbl_perf %>% bind_cols(tibble(
      auc_iregulon = output_iregulon$auc_iregulon,
      auc_iregulon_corrected = output_iregulon$auc_iregulon_corrected
    ))
  }
  return(tbl_perf)
}
classification_evaluation_categorical_pred = function(predictions, response) {
  if (sum(predictions) == 0){
    return(accuracy = 0,
           recall = 0,
           specificity = 0,
           precision = 0,
           F1 =  0,
           F05 = 0,
           F2 = 0,
           mcc = 0,
           informedness = 0,
           markedness = 0,
           fisher_pval_log = 0,
           fisher_odds = 1)
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

  fisher = fisher.test(as.factor(response), predictions)

  mcc_S = (tp + fn)/num_total
  mcc_P = (tp + fp)/num_total

  metrics = dplyr::tibble(
    accuracy = (tp + tn) / (num_positives + num_negatives),
    recall = tp / num_positives,
    specificity = tn / num_negatives,
    precision = tp / (tp + fp),
    F1 =  (2 * precision * recall) / (precision + recall),
    F05 = (1.5 * precision * recall) / (0.5 * precision + recall),
    F2 = (3 * precision * recall) / (2 * precision + recall),
    mcc = (tp/num_total - mcc_S * mcc_P)/sqrt(mcc_P * mcc_S * (1-mcc_S) * (1-mcc_P)) ,
    informedness = recall + specificity - 1,
    markedness = precision + npv - 1,
    fisher_pval_log = -log(fisher$p.value),
    fisher_odds = fisher$estimate
  )

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















# # popularity table
# gene2pubmed = read_tsv("ftp://ftp.ncbi.nih.gov/gene/DATA/gene2pubmed.gz", col_names = c("taxid", "entrez", "pubmedid"), skip = 1, col_types = cols(
#   taxid = col_integer(),
#   entrez = col_character(),
#   pubmedid = col_integer()
# ))
# ncitations = gene2pubmed %>% filter(taxid==9606) %>% group_by(entrez) %>% summarize(ncitations=n())
# ncitations_bins = ncitations %>% left_join(geneinfo_human, "entrez") %>% arrange(-ncitations) %>% filter(ncitations > 3) %>% mutate(bin = Hmisc::cut2(log(ncitations), m = 1000, g = 5) %>% as.numeric()) # 0,1 and 2 dont fall within

cal_path_performance_cluster_multiple_binned_targets = function(setting, ligand_to_target, ncitations, bin_id){
  if (length(setting$from) == 1){
    performance = cal_path_performance(setting$diffexp %>% filter(gene %in% (ncitations %>% filter(bin == bin_id) %>% .$entrez)), ligand_to_target[humansymbol2humanentrez[[setting$from]],])
    performance$setting = setting$name
    performance$ligand = setting$from
    performance$bin_id = bin_id
    performance$expression_type = setting$type

    performance = performance %>% tbl_df()
    performance
  } else {
    multi_ligands = humansymbol2humanentrez[setting$from]

    multi_ligands_name = paste0(multi_ligands,collapse = "-")

    multi_ligands_symbol = paste0(setting$from,collapse = "-")

    performance = cal_path_performance(setting$diffexp %>% filter(gene %in% (ncitations %>% filter(bin == bin_id) %>% .$entrez)), ligand_to_target[multi_ligands_name,])

    performance$setting = setting$name
    performance$ligand = multi_ligands_symbol
    performance$bin_id = bin_id
    performance$expression_type = setting$type

    performance = performance %>% tbl_df()

    performance
  }
}
add_popularity_measures_to_perfs = function(performances, settings,ligand_target_matrix, bin_function = cal_path_performance_cluster_multiple_binned_targets,ncitations_bins){
  performances = performances %>% drop_na(-ligand)
  if (nrow(performances) == 0) {
    performances = performances %>% mutate(ligand_pop_auc = NA,
                                           ligand_pop_aupr = NA,
                                           target_auc = NA,
                                           target_aupr = NA,
                                           target_specificity = NA,
                                           target_sensitivity = NA)
    return(performances)
  }



  # print(performances)

  performances_ligand_pop = performances %>% left_join(ncitations_bins %>% rename(ligand = symbol), by = "ligand")
  ligand_pop_regression_auc_log = lm(auc ~ log(ncitations),performances_ligand_pop %>% select(ncitations,auc))
  ligand_slope_auc_log =  summary(ligand_pop_regression_auc_log) %>% .$coefficients %>% .[2,1]
  ligand_pop_regression_aupr_log = lm(aupr ~ log(ncitations),performances_ligand_pop %>% select(ncitations,aupr))
  ligand_slope_aupr_log =  summary(ligand_pop_regression_aupr_log) %>% .$coefficients %>% .[2,1]

  # print(performances_ligand_pop)


  calculate_all_performances_lapply = function(i,settings,fun,ligand_to_target,ncitations_bins) {bind_rows(lapply(X = settings, FUN = fun, ligand_to_target, ncitations_bins,i))}
  performances_bins = bind_rows(lapply(1:5,calculate_all_performances_lapply,settings,bin_function,ligand_target_matrix,ncitations_bins))
  performances_bins = performances_bins %>% drop_na()

  if (nrow(performances_bins %>% group_by(bin_id) %>% count()) != 5){
    performances = performances %>% mutate(ligand_pop_auc = NA,
                                           ligand_pop_aupr = NA,
                                           target_auc = NA,
                                           target_aupr = NA,
                                           target_specificity = NA,
                                           target_sensitivity = NA)
    return(performances)
  }

  # print(performances_bins)
  target_pop_regression_auc = lm(auc ~ bin_id, performances_bins %>% select(auc,bin_id))
  target_slope_auc = summary(target_pop_regression_auc) %>% .$coefficients %>% .[2,1]
  target_pop_regression_aupr = lm(aupr ~ bin_id, performances_bins %>% select(aupr,bin_id))
  target_slope_aupr = summary(target_pop_regression_aupr) %>% .$coefficients %>% .[2,1]
  target_pop_regression_sensitivity = lm(sensitivity ~ bin_id, performances_bins %>% select(sensitivity,bin_id))
  target_slope_sensitivity = summary(target_pop_regression_sensitivity) %>% .$coefficients %>% .[2,1]
  target_pop_regression_specificity = lm(specificity ~ bin_id, performances_bins %>% select(specificity,bin_id))
  target_slope_specificity = summary(target_pop_regression_specificity) %>% .$coefficients %>% .[2,1]

  performances = performances %>% mutate(ligand_pop_auc = ligand_slope_auc_log,
                                         ligand_pop_aupr = ligand_slope_aupr_log,
                                         target_auc = target_slope_auc,
                                         target_aupr = target_slope_aupr,
                                         target_specificity = target_slope_specificity,
                                         target_sensitivity = target_slope_sensitivity)

}
add_ligand_popularity_measures_to_perfs = function(performances,ligand_target_matrix,ncitations_bins){
  performances = performances %>% drop_na()
  if (nrow(performances) == 0) {
    performances = performances %>% mutate(ligand_pop_auc = NA,
                                           ligand_pop_aupr = NA)
    return(performances)
  }
  performances_ligand_pop = performances %>% left_join(ncitations_bins %>% rename(ligand = symbol), by = "ligand")
  ligand_pop_regression_auc_log = lm(auc ~ log(ncitations),performances_ligand_pop %>% select(ncitations,auc))
  ligand_slope_auc_log =  summary(ligand_pop_regression_auc_log) %>% .$coefficients %>% .[2,1]
  ligand_pop_regression_aupr_log = lm(aupr ~ log(ncitations),performances_ligand_pop %>% select(ncitations,aupr))
  ligand_slope_aupr_log =  summary(ligand_pop_regression_aupr_log) %>% .$coefficients %>% .[2,1]
  performances = performances %>% mutate(ligand_pop_auc = ligand_slope_auc_log,
                                         ligand_pop_aupr = ligand_slope_aupr_log)
}

