scale_quantile_adapted = function(x){
  y = nichenetr::scale_quantile(x,outlier_cutoff = 0)
  y = y + 0.001
  return(y)
}
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

# scale_quantile_adapted = function(x){
#   y = nichenetr::scale_quantile(x,outlier_cutoff = 0)
#   y = y + 0.001
#   return(y)
# }
calculate_niche_de = function(seurat_obj, niches, type, assay_oi = "SCT"){

  if (type == "sender"){
    #  Determine DE between sender cell types across the different niches
    DE_sender = niches %>% lapply(function(niche, seurat_obj) {

      senders_niche = niche$sender %>% unlist() %>% unique()
      senders_other = niches %>% sapply(function(niche){niche$sender})  %>% unlist() %>% unique() %>% setdiff(senders_niche)

      DE_sender = senders_niche %>% lapply(function(sender_oi, seurat_obj, senders_other) {

        print(paste0("Calculate Sender DE between: ",sender_oi, " and ", senders_other))

        DE_subtable = senders_other %>% lapply(function(sender_other_niche, seurat_obj, sender_oi) {
          DE_sender_oi = FindMarkers(object = seurat_obj, ident.1 = sender_oi, ident.2 = sender_other_niche, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
          DE_sender_oi = DE_sender_oi %>% mutate(sender = sender_oi, sender_other_niche = sender_other_niche) %>% arrange(-avg_log2FC)
        }, seurat_obj, sender_oi) %>% bind_rows()

      }, seurat_obj, senders_other) %>% bind_rows()

    }, seurat_obj) %>% bind_rows()

    return(DE_sender)

  }
  if (type == "receiver"){
    #  Determine DE between receiver cell types across the different niches
    # # DE receiver: niche one vs other niches
    all_receivers = niches %>% purrr::map("receiver") %>% unlist() %>% unique()
    #
    DE_receiver = all_receivers %>% lapply(function(receiver_oi, seurat_obj) {

      receivers_other = setdiff(all_receivers,receiver_oi)

      print(paste0("Calculate Receiver DE between: ",receiver_oi, " and ", receivers_other))

      DE_subtable = receivers_other %>% lapply(function(receiver_other_niche, seurat_obj, receiver_oi) {

        DE_receiver_oi = FindMarkers(object = seurat_obj, ident.1 = receiver_oi, ident.2 = receiver_other_niche, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
        DE_receiver_oi = DE_receiver_oi %>% mutate(receiver = receiver_oi, receiver_other_niche = receiver_other_niche) %>% arrange(-avg_log2FC)

      }, seurat_obj, receiver_oi) %>% bind_rows()
    }, seurat_obj) %>% bind_rows()

    return(DE_receiver)

  }


}

calculate_receiver_de_targets  = function(seurat_obj, niches, assay_oi = "SCT", lfc_cutoff, expression_pct){
  all_receivers = niches %>% purrr::map("receiver") %>% unlist() %>% unique()
  #
  DE_receiver = all_receivers %>% lapply(function(receiver_oi, seurat_obj) {

    receivers_other = setdiff(all_receivers,receiver_oi)

    print(paste0("Calculate Receiver DE between: ",receiver_oi, " and ", receivers_other))

    DE_subtable = receivers_other %>% lapply(function(receiver_other_niche, seurat_obj, receiver_oi) {

      DE_receiver_oi = FindMarkers(object = seurat_obj, ident.1 = receiver_oi, ident.2 = receiver_other_niche, min.pct = expression_pct, logfc.threshold = lfc_cutoff, only.pos = TRUE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
      DE_receiver_oi = DE_receiver_oi %>% mutate(receiver = receiver_oi, receiver_other_niche = receiver_other_niche) %>% arrange(-avg_log2FC)

    }, seurat_obj, receiver_oi) %>% bind_rows()
  }, seurat_obj) %>% bind_rows()

}

process_niche_de = function(DE_table, niches, type, expression_pct){
  ### process the DE_tables of senders and receiver

  if(type == "sender"){
    DE_sender = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC) %>% arrange(-score)
    DE_sender_processed = DE_sender %>% group_by(gene, sender) %>% summarise(mean_avg_log2FC = mean(avg_log2FC), min_avg_log2FC = min(avg_log2FC), mean_significant = mean(significant), mean_present = mean(present), mean_score = mean(score), min_score = min(score)) %>% arrange(-min_avg_log2FC)
    DE_sender_processed = names(niches) %>% lapply(function(niche_name, niches){
      tibble(niche = niche_name, sender = niches[[niche_name]]$sender)
    }, niches) %>% bind_rows() %>% inner_join(DE_sender_processed, by = c("sender"))
    return(DE_sender_processed)

  }
  if(type == "receiver"){
    DE_receiver = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC) %>% arrange(-score)
    DE_receiver_processed = DE_receiver %>% group_by(gene, receiver) %>% summarise(mean_avg_log2FC = mean(avg_log2FC), min_avg_log2FC = min(avg_log2FC), mean_significant = mean(significant), mean_present = mean(present), mean_score = mean(score), min_score = min(score)) %>% arrange(-min_avg_log2FC)
    DE_receiver_processed = names(niches) %>% lapply(function(niche_name, niches){
      tibble(niche = niche_name, receiver = niches[[niche_name]]$receiver)
    }, niches) %>% bind_rows() %>% inner_join(DE_receiver_processed, by = c("receiver"))
    return(DE_receiver_processed)

  }

}

combine_sender_receiver_de = function(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = "min_lfc"){

  if(specificity_score == "min_lfc"){
    DE_sender_processed_ligands = DE_sender_processed %>% rename(ligand = gene, ligand_score = min_avg_log2FC, ligand_significant = mean_significant, ligand_present = mean_present)
    DE_receiver_processed_receptors = DE_receiver_processed %>% rename(receptor = gene, receptor_score = min_avg_log2FC, receptor_significant = mean_significant, receptor_present = mean_present)
  }
  if(specificity_score == "mean_lfc"){
    DE_sender_processed_ligands = DE_sender_processed %>% rename(ligand = gene, ligand_score = mean_avg_log2FC, ligand_significant = mean_significant, ligand_present = mean_present)
    DE_receiver_processed_receptors = DE_receiver_processed %>% rename(receptor = gene, receptor_score = mean_avg_log2FC, receptor_significant = mean_significant, receptor_present = mean_present)
  }
  if(specificity_score == "min_score"){
    DE_sender_processed_ligands = DE_sender_processed %>% rename(ligand = gene, ligand_score = min_score, ligand_significant = mean_significant, ligand_present = mean_present)
    DE_receiver_processed_receptors = DE_receiver_processed %>% rename(receptor = gene, receptor_score = min_score, receptor_significant = mean_significant, receptor_present = mean_present)
  }
  if(specificity_score == "mean_score"){
    DE_sender_processed_ligands = DE_sender_processed %>% rename(ligand = gene, ligand_score = mean_score, ligand_significant = mean_significant, ligand_present = mean_present)
    DE_receiver_processed_receptors = DE_receiver_processed %>% rename(receptor = gene, receptor_score = mean_score, receptor_significant = mean_significant, receptor_present = mean_present)
  }

  DE_sender_processed_ligands = DE_sender_processed_ligands %>% ungroup() %>% filter(ligand %in% lr_network$ligand) %>% dplyr::mutate(scaled_ligand_score = scale_quantile_adapted(ligand_score))
  DE_receiver_processed_receptors = DE_receiver_processed_receptors %>% ungroup() %>% filter(receptor %in% lr_network$receptor) %>% dplyr::mutate(scaled_receptor_score = scale_quantile_adapted(receptor_score))

  DE_sender_receiver = lr_network %>%
    inner_join(DE_sender_processed_ligands, by = c("ligand")) %>%
    inner_join(DE_receiver_processed_receptors, by = c("receptor","niche"))

  DE_sender_receiver = DE_sender_receiver %>% select(niche, sender, receiver, ligand, receptor, ligand_score, ligand_significant, ligand_present, receptor_score, receptor_significant, receptor_present, scaled_ligand_score, scaled_receptor_score) %>% mutate(avg_score_ligand_receptor = 0.5*(ligand_score + receptor_score)) %>% arrange(-avg_score_ligand_receptor) %>% ungroup()
  DE_sender_receiver = DE_sender_receiver %>% dplyr::mutate(scaled_avg_score_ligand_receptor = scale_quantile_adapted(avg_score_ligand_receptor))
  return(DE_sender_receiver)
}
process_receiver_target_de = function(DE_receiver_targets, niches, expression_pct, specificity_score = "min_lfc"){

  DE_receiver = DE_receiver_targets %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC) %>% arrange(-score)

  DE_receiver_processed = DE_receiver %>% group_by(gene, receiver) %>% summarise(mean_avg_log2FC = mean(avg_log2FC), min_avg_log2FC = min(avg_log2FC), mean_significant = mean(significant), mean_present = mean(present), mean_score = mean(score), min_score = min(score)) %>% arrange(-min_avg_log2FC)
  DE_receiver_processed = names(niches) %>% lapply(function(niche_name, niches){
    tibble(niche = niche_name, receiver = niches[[niche_name]]$receiver)
  }, niches) %>% bind_rows() %>% inner_join(DE_receiver_processed, by = c("receiver"))


  if(specificity_score == "min_lfc"){
    DE_receiver_processed_targets = DE_receiver_processed %>% rename(target = gene, target_score = min_avg_log2FC, target_significant = mean_significant, target_present = mean_present)
  }
  if(specificity_score == "mean_lfc"){
    DE_receiver_processed_targets = DE_receiver_processed %>% rename(target = gene, target_score = mean_avg_log2FC, target_significant = mean_significant, target_present = mean_present)
  }
  if(specificity_score == "min_score"){
    DE_receiver_processed_targets = DE_receiver_processed %>% rename(target = gene, target_score = min_score, target_significant = mean_significant, target_present = mean_present)
  }
  if(specificity_score == "mean_score"){
    DE_receiver_processed_targets = DE_receiver_processed %>% rename(target = gene, target_score = mean_score, target_significant = mean_significant, target_present = mean_present)
  }

  DE_receiver_processed_targets = DE_receiver_processed_targets %>% select(niche, receiver, target, target_score, target_significant, target_present)  %>% arrange(-target_score)


  return(DE_receiver_processed_targets)
}
get_ligand_activities_targets = function(niche_geneset_list, ligand_target_matrix, top_n_target){

  ligand_activities_targets = niche_geneset_list %>% lapply(function(niche_oi, ligand_target_matrix, top_n_target){
    receiver_oi = niche_oi$receiver

    background_expressed_genes= niche_oi$background
    geneset_oi = niche_oi$geneset

    background_expressed_genes = background_expressed_genes %>% dplyr::intersect(rownames(ligand_target_matrix))
    geneset_oi = geneset_oi %>% dplyr::intersect(rownames(ligand_target_matrix))

    ligand_target_matrix = ligand_target_matrix[rownames(ligand_target_matrix) %in% background_expressed_genes, ]
    ligands = colnames(ligand_target_matrix)

    print(paste0("Calculate Ligand activities for: ",receiver_oi))

    ligand_activities = nichenetr::predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands)
    ligand_activities = ligand_activities  %>% dplyr::rename(ligand = test_ligand, activity = pearson) %>% dplyr::select(-aupr, -auroc) %>% filter(!is.na(activity))
    ligand_target_df = ligand_activities$ligand %>% unique() %>% lapply(nichenetr::get_weighted_ligand_target_links, geneset_oi, ligand_target_matrix, top_n_target) %>% dplyr::bind_rows()  %>% dplyr::rename(ligand_target_weight = weight)
    ligand_activities = ligand_activities %>% dplyr::inner_join(ligand_target_df, by = c("ligand")) %>% dplyr::mutate(receiver = receiver_oi) %>% dplyr::group_by(receiver) %>% dplyr::mutate(activity_normalized = nichenetr::scaling_zscore(activity))

  }, ligand_target_matrix, top_n_target) %>% bind_rows()

  scaled_ligand_activities_targets = ligand_activities_targets %>% ungroup() %>% distinct(ligand, activity, activity_normalized, receiver) %>% mutate(scaled_activity_normalized = scale_quantile_adapted(activity_normalized), scaled_activity = scale_quantile_adapted(activity))
  ligand_activities_targets = ligand_activities_targets %>% inner_join(scaled_ligand_activities_targets, by = c("ligand","activity","activity_normalized","receiver")) %>% ungroup()

  return(ligand_activities_targets)
}
calculate_spatial_DE = function(seurat_obj, spatial_info){
  spatial_info$celltype_region_oi %>% lapply(function(celltype_oi, seurat_obj, spatial_info){
    other_celltype = spatial_info %>% filter(celltype_region_oi == celltype_oi) %>% pull(celltype_other_region) %>% unique()
    niche_oi = spatial_info %>% filter(celltype_region_oi == celltype_oi) %>% pull(niche) %>% unique()

    print(paste0("Calculate Spatial DE between: ",celltype_oi, " and ", other_celltype))

    DE_table = FindMarkers(object = seurat_obj, ident.1 = celltype_oi, ident.2 = other_celltype, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
    DE_table = DE_table %>% mutate(celltype = celltype_oi, niche = niche_oi) %>% arrange(-avg_log2FC)
  }, seurat_obj, spatial_info) %>% bind_rows()
}
process_spatial_de = function(DE_table, type, lr_network, expression_pct, specificity_score = "lfc"){
  ### process the DE_tables of senders and receiver

  if(type == "sender"){
    DE_sender = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC)
    if(specificity_score == "lfc"){
      DE_sender = DE_sender %>% rename(ligand = gene, ligand_score_zonation = avg_log2FC, sender = celltype)
    }
    if(specificity_score == "score"){
      DE_sender = DE_sender %>% rename(ligand = gene, ligand_score_zonation = score, sender = celltype)
    }
    DE_sender = DE_sender %>% select(niche, sender, ligand, ligand_score_zonation)  %>% arrange(-ligand_score_zonation) %>% filter(ligand %in% lr_network$ligand)
    return(DE_sender)
  }
  if(type == "receiver"){
    DE_receiver = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC)
    if(specificity_score == "lfc"){
      DE_receiver = DE_receiver %>% rename(receptor = gene, receptor_score_zonation = avg_log2FC, receiver = celltype)
    }
    if(specificity_score == "score"){
      DE_receiver = DE_receiver %>% rename(receptor = gene, receptor_score_zonation = score, receiver = celltype)
    }
    DE_receiver = DE_receiver %>% select(niche, receiver, receptor, receptor_score_zonation)  %>% arrange(-receptor_score_zonation) %>% filter(receptor %in% lr_network$receptor)
    return(DE_receiver)
  }
}
get_non_spatial_de = function(niches, spatial_info, type, lr_network){

  if(type == "sender"){

    niche_df = niches %>% names() %>% lapply(function(niche_oi, niches){
      tibble(niche = niche_oi, sender = niches[[niche_oi]]$sender)
    }, niches) %>% bind_rows()
    spatial_df = niches %>% purrr::map("sender") %>% unlist() %>% unique() %>% setdiff(spatial_info %>% filter(celltype_type == "sender") %>% pull(celltype_region_oi)) %>% lapply(function(sender_oi, niches, lr_network){
      spatial_df =  tibble(ligand = lr_network$ligand %>% unique(), ligand_score_zonation = 0)
      spatial_df = spatial_df %>% mutate(sender = sender_oi)
    },  niches, lr_network) %>% bind_rows() %>% select(sender, ligand, ligand_score_zonation)  %>% inner_join(niche_df, by = "sender")

  }
  if(type == "receiver"){

    niche_df = niches %>% names() %>% lapply(function(niche_oi, niches){
      tibble(niche = niche_oi, receiver = niches[[niche_oi]]$receiver)
    }, niches) %>% bind_rows()
    spatial_df = niches %>% purrr::map("receiver") %>% unlist() %>% unique() %>% setdiff(spatial_info %>% filter(celltype_type == "receiver") %>% pull(celltype_region_oi)) %>% lapply(function(receiver_oi, niches, lr_network){
      spatial_df =  tibble(receptor = lr_network$receptor %>% unique(), receptor_score_zonation = 0)
      spatial_df = spatial_df %>% mutate(receiver = receiver_oi)
    },  niches, lr_network) %>% bind_rows() %>% select(receiver, receptor, receptor_score_zonation)  %>% inner_join(niche_df, by = "receiver")


  }
  return(spatial_df)
}
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}
get_prioritization_tables = function(output_nichenet_analysis, prioritizing_weights){
  combined_information = output_nichenet_analysis$DE_sender_receiver %>%
    inner_join(output_nichenet_analysis$ligand_scaled_receptor_expression_fraction_df, by = c("receiver", "ligand", "receptor")) %>%
    inner_join(output_nichenet_analysis$sender_spatial_DE_processed, by = c("niche", "sender", "ligand")) %>%
    inner_join(output_nichenet_analysis$receiver_spatial_DE_processed, by = c("niche", "receiver", "receptor")) %>%
    inner_join(output_nichenet_analysis$ligand_activities_targets, by = c("receiver", "ligand")) %>%
    inner_join(output_nichenet_analysis$DE_receiver_processed_targets, by = c("niche", "receiver", "target")) %>%
    inner_join(output_nichenet_analysis$exprs_tbl_ligand, by = c("sender", "ligand")) %>%
    inner_join(output_nichenet_analysis$exprs_tbl_receptor, by = c("receiver", "receptor")) %>%
    inner_join(output_nichenet_analysis$exprs_tbl_target, by = c("receiver", "target"))

  # reorder the columns



  combined_information = combined_information %>% mutate(ligand_receptor = paste(ligand, receptor, sep = "--"))  %>%  mutate(bonafide_score = 1) %>%  mutate_cond(bonafide == FALSE, bonafide_score = 0.5)


  combined_information = combined_information %>% select(
    niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
    ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_zonation,
    receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_zonation,
    ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor, bonafide_score,
    target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
    activity, activity_normalized,
    scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
    scaled_ligand_score_zonation, scaled_receptor_score_zonation,
    scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
    scaled_activity, scaled_activity_normalized
  ) %>% distinct()


  # combine in a prioritization score



  combined_information_prioritized = combined_information %>%
    dplyr::mutate(prioritization_score =
                    ((prioritizing_weights["scaled_ligand_score"] * scaled_ligand_score) +
                       (prioritizing_weights["scaled_ligand_expression_scaled"] * scaled_ligand_expression_scaled) +
                       (prioritizing_weights["scaled_receptor_score"] * scaled_receptor_score) +
                       (prioritizing_weights["scaled_receptor_expression_scaled"] * scaled_receptor_expression_scaled) +
                       (prioritizing_weights["scaled_ligand_score_zonation"] * scaled_ligand_score_zonation) +
                       (prioritizing_weights["scaled_receptor_score_zonation"] * scaled_receptor_score_zonation) +
                       (prioritizing_weights["ligand_scaled_receptor_expression_fraction"] * ligand_scaled_receptor_expression_fraction) +
                       (prioritizing_weights["scaled_activity"] * scaled_activity) +
                       (prioritizing_weights["scaled_activity_normalized"] * scaled_activity_normalized) +
                       (prioritizing_weights["ligand_fraction"] * scaled_ligand_fraction_adapted ) +
                       (prioritizing_weights["receptor_fraction"] * scaled_receptor_fraction_adapted  ) +
                       (prioritizing_weights["bona_fide"] * bonafide_score)
                    )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)

  prioritization_tbl_ligand_receptor = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide,
                                                                                   ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_zonation,
                                                                                   receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_zonation,
                                                                                   ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor,
                                                                                   activity, activity_normalized,
                                                                                   scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
                                                                                   scaled_ligand_score_zonation, scaled_receptor_score_zonation,
                                                                                   scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
                                                                                   scaled_activity, scaled_activity_normalized,
                                                                                   prioritization_score) %>% distinct()
  prioritization_tbl_ligand_target = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
                                                                                 target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
                                                                                 activity, activity_normalized, scaled_activity, scaled_activity_normalized, prioritization_score) %>% distinct()

  return(list(prioritization_tbl_ligand_receptor = prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target = prioritization_tbl_ligand_target))

}
