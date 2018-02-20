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
classification_evaluation_continuous_pred = function(prediction,response, iregulon = TRUE){

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

evaluate_target_prediction_bins_direct = function(bin_id,settings,ligand_target_matrix,ligands_position,ncitations){
  requireNamespace("dplyr")
  targets_bin_id = ncitations %>% filter(bin == bin_id) %>% .$symbol
  if (ligands_position == "cols"){
    ligand_target_matrix = ligand_target_matrix[targets_bin_id,]
  } else if (ligands_position == "rows") {
    ligand_target_matrix = ligand_target_matrix[,targets_bin_id]
  }
  performances = bind_rows(lapply(settings, evaluate_target_prediction, ligand_target_matrix,ligands_position)) %>% mutate(target_bin_id = bin_id)
}


# cal_path_performance_cellcell = function(diffexp, prior){
#   # print(head(diffexp))
#   # print(prior[1:15,1:8])
#   LT_df = prior %>% as.matrix() %>% t() %>% data.frame() %>% rownames_to_column("gene") %>% tbl_df()
#
#   E_df_classification = diffexp %>% tbl_df() %>% mutate(diffexp = abs(lfc) > 1 & qval <= 0.1) %>% rename(lfc_cellcell = diffexp) %>% select(gene,lfc_cellcell)
#   # E_df_classification = E_df_classification %>% mutate(gene = [gene]) %>% drop_na()
#   all_data_classification = inner_join(E_df_classification, LT_df)
#
#   Y_classification = all_data_classification[,2:ncol(E_df_classification)] %>% as.matrix()
#   X_classification = all_data_classification[,(ncol(E_df_classification)+1):ncol(all_data_classification)] %>% as.matrix()
#
#   # print(all_data_classification)
#   # print(all_data_classification %>% filter(gene == "1"))
#   #
#   # print(X_classification[1:5,1:5])
#   # print(Y_classification[1:5,])
#   # print(Y_classification)
#   # print(sum(Y_classification))
#
#   # doMC::registerDoMC(cores = 2)
#   # algorithm = "glmnet"
#   # control =  caret::trainControl(method="repeatedcv", number=5, repeats=1, classProbs = TRUE, summaryFunction = twoClassSummary)
#   # model =  caret::train(y = Y_classification %>% make.names() %>% as.factor(),x = X_classification, method=algorithm, trControl=control, metric = 'ROC')
#   # # imps =  caret::varImp(model, scale = FALSE)
#   # auroc_glmnet =  model$results$ROC %>% max()
#
#   # algorithm = "rf"
#   Grid <-  expand.grid(mtry = round(length(colnames(X_classification)) ** 0.33))
#   control =  caret::trainControl(method="repeatedcv", number=5, repeats=1, classProbs = TRUE, summaryFunction = twoClassSummary)
#   model =  caret::train(y = Y_classification %>% make.names() %>% as.factor(),x = X_classification, method=algorithm, trControl=control, metric = 'ROC',tuneGrid = Grid)
#   auroc_rf =  model$results$ROC %>% max()
#   # return(list(auroc_glmnet = auroc_glmnet, auroc_rf = auroc_rf))
#   # return(list(auroc_glmnet = auroc_glmnet))
#   return(list(auroc_glmnet = auroc_rf, auroc_rf = auroc_rf))
#
# }



# trainIndex = caret::createDataPartition(y, p = 0.8) # y should be factor; p=0.8: 0.8-0.2 split
# caret::train: estimate model performance on training set: bootstrap instead repeatedcv?
# want specific splits --> index argument of trainControl
# index and indexOut: list with elements for each resampling iteration: sample rowss used for traning
# preProcess in caret function train
# traincontrol: none? no cv; but better doing so ... (altought takes longer --> so single-cell: no CV)
# allowParallel in trainControl
# user-specific performance metrics in trainControl: elements of function: "data" with columns: obs and pred; "lev": outcome factor levels; "model"; output = named vector of numeric summary metrics
##
# caret::twoClassSummary
# function (data, lev = NULL, model = NULL)
# {
#   lvls <- levels(data$obs)
#   if (length(lvls) > 2)
#     stop(paste("Your outcome has", length(lvls), "levels. The twoClassSummary() function isn't appropriate."))
#   requireNamespaceQuietStop("ModelMetrics")
#   if (!all(levels(data[, "pred"]) == lvls))
#     stop("levels of observed and predicted data do not match")
#   rocAUC <- ModelMetrics::auc(ifelse(data$obs == lev[2], 0,
#                                      1), data[, lvls[1]])
#   out <- c(rocAUC, sensitivity(data[, "pred"], data[, "obs"],
#                                lev[1]), specificity(data[, "pred"], data[, "obs"], lev[2]))
#   names(out) <- c("ROC", "Sens", "Spec")
#   out
# }

## extracting predictions from train object in finalModel slot!
# predict(object, type = "prob") # or type = "class"
# if windows = FALSE: library(doMC) registerDoMC

#########
# variable importances: varImp.train: set to FALSE
# linear models - RF - PLS -  gbm - pam?
# classification without varimp: nbayes - lda - nnet - PenalizedLDA
# glmnet?
# varImp(fit, scale = FALSE)
