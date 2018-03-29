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
    tbl_perf = tibble(auroc = auroc,
                      aupr = aupr,
                      aupr_corrected = aupr - aupr_random,
                      sensitivity_roc = sensitivity,
                      specificity_roc = specificity,
                      mean_rank_GST_log_pval = -log(mean_rank_GST),
                      auc_iregulon = output_iregulon$auc_iregulon,
                      auc_iregulon_corrected = output_iregulon$auc_iregulon_corrected,
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
    F05 = (1.5 * precision * recall) / (0.5 * precision + recall),
    F2 = (3 * precision * recall) / (2 * precision + recall),
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
  targets_bin_id = ncitations %>% filter(bin == bin_id) %>% .$symbol
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
wrapper_caret_classification = function(train_data, algorithm, continuous = TRUE, var_imps = TRUE, cv = TRUE, cv_number = 5, cv_repeats = 2, parallel = FALSE, n_cores = 4,prediction_response_df = NULL,ignore_errors = FALSE, return_model = FALSE){

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
    var_imps_df = imps$importance  %>% tibble::rownames_to_column("feature") %>% tbl_df() %>% .[,1:2]
    colnames(var_imps_df) = c("feature","importance")
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, performances_training_continuous = performances_training_continuous ,var_imps = var_imps_df, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    } else {
    output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, performances_training_continuous = performances_training_continuous, var_imps = var_imps_df)}

  } else {
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, performances_training_continuous = performances_training_continuous, var_imps = NULL, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    } else {
    output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training,performances_training_continuous = performances_training_continuous, var_imps = NULL)}
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
  test_ligand = importances$test_ligand
  real_ligand = strsplit(importances$ligand,"[-]")
  true_ligand = rep(NULL, times = length(test_ligand))
  for (i in seq(length(test_ligand))){
    test = test_ligand[i]
    real = real_ligand[[i]]
    true_ligand[i] = test %in% real
  }
  return(true_ligand)
}

scaling_zscore = function(x){
  if (typeof(x) == "double"){
    if(sd(x) > 0){
      return((x - mean(x))/sd(x))
    } else{
      return((x - mean(x)))
    }
  } else {return(x)}
}
scaling_modified_zscore = function(x){
  if (typeof(x) == "double"){
    if(mad(x) > 0){
      return(0.6745*(x - median(x))/mad(x))
    } else{
      return(0.6745*(x - median(x)))
    }
  } else {return(x)}
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

  tbl_perf = tibble(r_squared = model_summary$r.squared,
                    adj_r_squared = model_summary$adj.r.squared,
                    f_statistic = model_summary$fstatistic[1],
                    lm_coefficient_abs_t = model_summary$coefficients %>% .[2,3] %>% abs(),
                    inverse_rmse = 1/(model_summary$sigma),
                    reverse_aic = AIC(model) * -1,
                    reverse_bic = BIC(model) * -1,
                    inverse_mae = 1/(mae(model$residuals)),
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
wrapper_caret_regression = function(train_data, algorithm, var_imps = TRUE, cv = TRUE, cv_number = 5, cv_repeats = 2, parallel = FALSE, n_cores = 4,prediction_response_df = NULL,ignore_errors = FALSE, return_model = FALSE){

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
                            maximize = TRUE) # RMSE should be minimized - so maximize its inverse
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
                           maximize = TRUE) # RMSE should be minimized - so maximize its inverse
      }


  performances =  model$resample

  final_model_predictions = predict(model,newdata = train_data)
  response = train_data$obs
  performances_training = regression_evaluation(prediction = final_model_predictions, response = response)

  if (var_imps == TRUE) {
    imps =  caret::varImp(model, scale = FALSE)
    var_imps_df = imps$importance  %>% tibble::rownames_to_column("feature") %>% tbl_df() %>% .[,1:2]
    colnames(var_imps_df) = c("feature","importance")
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, var_imps = var_imps_df, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    }
    output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, var_imps = var_imps_df)

  } else {
    if(!is.null(prediction_response_df)){
      output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, var_imps = NULL, prediction_response_df = prediction_response_df %>% mutate(model = final_model_predictions))
    }
    output_list = list(performances = performances %>% tbl_df(), performances_training = performances_training, var_imps = NULL)
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



  ligand_tf_matrix_visualization_df = tbl_df(ligand_tf_matrix_visualization) %>% rename(weight = V1) %>% mutate(TF = rownames(ligand_tf_matrix)) %>% select(TF,weight)
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
