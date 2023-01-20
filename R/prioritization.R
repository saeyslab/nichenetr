scale_quantile_adapted = function(x, outlier_cutoff = 0){
  y = nichenetr::scale_quantile(x,outlier_cutoff  = outlier_cutoff)
  y = y + 0.001
  return(y)
}

#' @title Calculate differential expression of one cell type versus all other cell types
#'
#' @description \code{calculate_niche_de} Calculate differential expression of one cell type versus all other cell types. This is possible for sender cell types and receiver cell types.
#'
#' @usage
#' calculate_de(seurat_obj, celltype_id, senders_oi = NA, receivers_oi = NA, assay_oi = "RNA")
#'
#' @param seurat_obj Seurat object
#' @param celltype_id Name of the meta data column that indicates the cell type of a cell
#' @param group_id group_id Name of the meta data column that indicates from which group/condition a cell comes from
#' @param assay_oi Which assay need to be used for DE calculation via `FindMarkers`. Default RNA, alternatives: SCT.
#'
#' @return A tibble containing the DE results
#'
#' @examples
#' \dontrun{
#' TODO
#' seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
#' niches = list(
#' "KC_niche" = list(
#'   "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
#'   "receiver" = c("KCs")),
#' "MoMac2_niche" = list(
#'   "sender" = c("Cholangiocytes","Fibroblast 2"),
#'   "receiver" = c("MoMac2")),
#' "MoMac1_niche" = list(
#'   "sender" = c("Capsule fibroblasts","Mesothelial cells"),
#'  "receiver" = c("MoMac1"))
#' )
#' calculate_niche_de(seurat_obj, niches, "sender")
#' }
#'
#' @export
#'
calculate_de = function(seurat_obj, celltype_id,
                        condition_oi = NA, group_id = NA,
                        min.pct = 0.1, assay_oi = "RNA"){

  # if (!is.na(group_id) & (is.na(condition_reference) || is.na(condition_oi))){
  #   stop("Please input both condition_reference and condition_oi")
  # }

  if (any(!is.na(group_id), !is.na(condition_oi)) & !all(!is.na(group_id), !is.na(condition_oi))){
    stop("Please input both group_id and condition_oi")
  }

  celltypes <- unique(seurat_obj[[celltype_id, drop=TRUE]])

  if (!is.na(condition_oi)) {
    seurat_obj = seurat_obj[,seurat_obj[[group_id]] == condition_oi]
  }

  DE_table = FindAllMarkers(seurat_obj, min.pct = min.pct) %>%
                rename(cluster_id = cluster)


  # # If there is more than one condition, calculate condition-celltype specificity
  # if (!is.na(group_id)){
  #   DE_table = lapply(celltypes, function(ct) {
  #     seurat_obj_subset = subset(seurat_obj, idents = ct) %>% SetIdent(value = .[[group_id]])
  #     FindMarkers(object = seurat_obj_subset, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = min.pct) %>%
  #       rownames_to_column("gene") %>% mutate(cluster_id = ct)
  #   }) %>% do.call(rbind, .)
  # # If there is only one condition, calculate cell type specificity
  # } else {
  #   seurat_obj <- SetIdent(seurat_obj, value = .[[celltype_id]])
  #   DE_table = FindAllMarkers(seurat_obj, min.pct = min.pct) %>%
  #               rename(cluster_id = cluster)
  # }
  #
  # SeuratV4 = c("avg_log2FC") %in% colnames(DE_table)
  # if(!SeuratV4){
  #   DE_table = DE_table %>% dplyr::rename(avg_log2FC = avg_logFC)
  # }


  return(DE_table)


}
#' @title get_exprs_avg
#'
#' @description \code{get_exprs_avg}  Calculate (group-)average of gene expression per cell type.
#' @usage get_muscat_exprs_avg(sce, celltype_id, group_id)
#'
#' @return Data frame with average gene expression per sample and per group.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom Seurat NormalizeData
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' muscat_exprs_avg = get_muscat_exprs_avg(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' }
#'
#' @export
#'
get_exprs_avg = function(seurat_obj, celltype_id, lr_network,
                         condition_oi = NA, group_id = NA){

  requireNamespace("dplyr")

  seurat_obj[[celltype_id]] <- make.names(seurat_obj[[celltype_id, drop=TRUE]])
  if (!is.na(group_id)) seurat_obj[[group_id]] <- make.names(seurat_obj[[group_id, drop=TRUE]])

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  avg_list <- list()
  avg_celltype <- AverageExpression(seurat_obj, assay = "RNA", slot = "data",
                    group.by = celltype_id) %>%
                    .$RNA %>% data.frame()
  avg_list$celltype <- avg_celltype %>% rownames_to_column("gene") %>%
                        pivot_longer(!gene, names_to = "celltype", values_to = "avg_expr")

  # If there is more than one condition, calculate also condition-specific expressions
  if (!is.na(group_id)) {
    avg_group <- AverageExpression(seurat_obj, assay = "RNA", slot = "data",
                             group.by = c(celltype_id, group_id)) %>%
      .$RNA %>% data.frame()

    pattern <- paste0("(", paste0(unique(seurat_obj[[celltype_id, drop=TRUE]]), collapse="|"), ")_(",
                      paste0(unique(seurat_obj[[group_id, drop=TRUE]]), collapse="|"), ")")

    avg_list$group <- avg_group %>% rownames_to_column("gene") %>%
      pivot_longer(!gene, names_to = c("celltype", "group"), values_to = "avg_expr", names_pattern = pattern)
  }

  ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
  receptors = lr_network %>% dplyr::pull(receptor) %>% unique()

  sender_list <- lapply(avg_list, function(u) u %>% filter(gene %in% ligands) %>%
                         dplyr::rename(sender = celltype, ligand = gene, avg_ligand = avg_expr))

  receiver_list <- lapply(avg_list, function(u) u %>% filter(gene %in% receptors) %>%
                           dplyr::rename(receiver = celltype, receptor = gene, avg_receptor = avg_expr))

  avg_sender_receiver <- lapply(names(avg_list), function(u) {
    columns_select <- c("sender", "receiver", "ligand", "receptor", "avg_ligand", "avg_receptor", "ligand_receptor_prod")
    avg_df_sender_receiver = sender_list[[u]] %>% dplyr::inner_join(lr_network, by = "ligand") %>% dplyr::inner_join(receiver_list[[u]], by = c("receptor"))

    if (u == "group"){
      avg_df_sender_receiver = avg_df_sender_receiver %>% filter(group.x == group.y) %>% select(-group.y) %>% rename(group = group.x)
      columns_select <- c("group", columns_select)
    }

    avg_df_sender_receiver = avg_df_sender_receiver %>% dplyr::mutate(ligand_receptor_prod = avg_ligand * avg_receptor) %>% dplyr::arrange(-ligand_receptor_prod) %>%
      dplyr::select(all_of(columns_select)) %>% dplyr::distinct()

  }) %>% setNames(names(avg_list))


  return (avg_sender_receiver)


}
#' @title process_info_to_ic
#'
#' @description \code{process_info_to_ic}  Process cell type expression information into intercellular communication focused information. Only keep information of ligands for the sender cell type setting, and information of receptors for the receiver cell type.
#' @usage process_info_to_ic(info_object, ic_type = "sender", lr_network)
#'
#' @param info_object Output of `get_exprs_avg`
#' @param ic_type "sender" or "receiver": indicates whether we should keep ligands or receptors respectively.
#'
#' @return List with expression information of ligands (sender case) or receptors (receiver case).
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' }
#'
#' @export
#'
process_info_to_ic = function(info_object, ic_type, lr_network){

  requireNamespace("dplyr")

  ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
  receptors = lr_network %>% dplyr::pull(receptor) %>% unique()

  if(ic_type == "sender"){
    avg_df = info_object$avg_df %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, avg_ligand = average_sample)

    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% ligands) %>% dplyr::rename(sender = celltype, ligand = gene, avg_ligand_group = average_group)

  }
  if(ic_type == "receiver"){
    frq_df = info_object$frq_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor = fraction_sample)
    pb_df = info_object$pb_df %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, pb_receptor = pb_sample)

    avg_df_group = info_object$avg_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, avg_receptor_group = average_group)
    frq_df_group = info_object$frq_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, fraction_receptor_group = fraction_group)
    pb_df_group = info_object$pb_df_group %>% dplyr::filter(gene %in% receptors) %>% dplyr::rename(receiver = celltype, receptor = gene, pb_receptor_group = pb_group)

    rel_abundance_df = info_object$rel_abundance_df %>% dplyr::rename(receiver = celltype, rel_abundance_scaled_receiver = rel_abundance_scaled)
  }

  return(list(avg_df = avg_df, frq_df = frq_df, pb_df = pb_df,  avg_df_group = avg_df_group, frq_df_group = frq_df_group, pb_df_group = pb_df_group, rel_abundance_df = rel_abundance_df))
}

#' @title combine_sender_receiver_de
#'
#' @description \code{combine_sender_receiver_de}  Combine FindMarkers differential expression output for senders and receivers by linking ligands to receptors based on the prior knowledge ligand-receptor network.
#' @usage combine_sender_receiver_de(DE_table_sender, DE_table_receiver, lr_network)
#'
#' @param DE_table_sender Differential expression analysis output for the sender cell types from the output of `Seurat::FindMarkers`.
#' @param DE_table_receiver Differential expression analysis output for the receiver cell types from the output of `Seurat::FindMarkers`.
#'
#' @return Data frame combining Seurat:: FindMarkers DE output for sender and receiver linked to each other through joining by the ligand-receptor network.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to)
#' condition_col = "aggregate"
#' condition_oi = "LCMV"
#' condition_reference = "SS"
#' senders_oi = c("CD4 T","Treg", "Mono", "NK", "B", "DC")
#' receivers_oi = "CD8 T"
#' DE_table_receiver = FindMarkers(object = subset(seuratObj, idents = receivers_oi) %>% SetIdent(.[[condition_col]]),
#'                                 ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene") %>%
#'                                 mutate(cluster_id = receiver)
#' DE_table_sender = lapply(senders_oi, function(sender) {
#'                          seurat_obj_sender = subset(seuratObj, idents = sender) %>% SetIdent(value = .[[condition_col]])
#'                          FindMarkers(object = seurat_obj_sender, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene") %>%
#'                            mutate(cluster_id = sender)
#' }) %>% do.call(rbind, .)
#' sender_receiver_de = combine_sender_receiver_de(
#'  DE_table_sender = DE_table_sender,
#'  DE_table_receiver = DE_table_receiver,
#'  lr_network = lr_network)
#' }
#'
#' @export
#'
combine_sender_receiver_de = function(DE_table_sender, DE_table_receiver, senders_oi, receivers_oi, lr_network){

  requireNamespace("dplyr")

  de_output_tidy_sender = DE_table_sender %>% dplyr::filter(cluster_id %in% senders_oi)
  de_output_tidy_receiver = DE_table_receiver %>% dplyr::filter(cluster_id %in% receivers_oi)

  de_output_tidy_sender = de_output_tidy_sender %>% dplyr::rename(ligand = gene, lfc_ligand = avg_log2FC, p_val_ligand = p_val,  p_adj_ligand = p_val_adj, sender = cluster_id)
  de_output_tidy_receiver =  de_output_tidy_receiver %>% dplyr::rename(receptor = gene, lfc_receptor = avg_log2FC, p_val_receptor = p_val,  p_adj_receptor = p_val_adj, receiver = cluster_id)

  de_tbl_sender_receiver = de_output_tidy_sender %>% dplyr::inner_join(lr_network, by = "ligand") %>%
    dplyr::inner_join(de_output_tidy_receiver, by = c("receptor"))

  de_tbl_sender_receiver = de_tbl_sender_receiver %>% dplyr::mutate(ligand_receptor_lfc_avg = (lfc_receptor + lfc_ligand)/2) %>% dplyr::arrange(-ligand_receptor_lfc_avg) %>%
    dplyr::select(sender, receiver, ligand, receptor, lfc_ligand, lfc_receptor, ligand_receptor_lfc_avg, p_val_ligand, p_adj_ligand, p_val_receptor, p_adj_receptor) %>% dplyr::distinct()

  return(de_tbl_sender_receiver)
}


#' @title generate_prioritization_tables
#'
#' @description \code{generate_prioritization_tables}  Perform a prioritization of cell-cell interactions (similar to MultiNicheNet).
#' User can choose the importance attached to each of the following prioritization criteria: differential expression of ligand and receptor, cell-type-(condition-)specificity of expression of ligand and receptor, NicheNet ligand activity
#' @usage generate_prioritization_tables(sender_receiver_info, sender_receiver_de, ligand_activities_targets_DEgenes, contrast_tbl, sender_receiver_tbl, grouping_tbl, prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2), fraction_cutoff, abundance_data_receiver, abundance_data_sender)
#'
#' @inheritParams combine_sender_receiver_info_ic
#' @param sender_receiver_info Output of `combine_sender_receiver_info_ic`
#' @param sender_receiver_de Output of `combine_sender_receiver_de`
#' @param ligand_activities_targets_DEgenes Output of `get_ligand_activities_targets_DEgenes`
#' @param sender_receiver_tbl Data frame with all sender-receiver cell type combinations (columns: sender and receiver)
#' @param grouping_tbl Data frame showing the groups of each sample (and batches per sample if applicable) (columns: sample and group; and if applicable all batches of interest)
#' @param abundance_data_receiver Data frame with number of cells per cell type - sample combination;  output of `process_info_to_ic`
#' @param abundance_data_sender Data frame with number of cells per cell type - sample combination; output of `process_info_to_ic`
#'
#' @return List containing multiple data frames prioritized senderLigand-receiverReceptor interactions (with sample- and group-based expression information), ligand activities and ligand-target links.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' sample_id = "tumor"
#' group_id = "pEMT"
#' celltype_id = "celltype"
#' batches = NA
#' contrasts_oi = c("'High-Low','Low-High'")
#' contrast_tbl = tibble(contrast = c("High-Low","Low-High"), group = c("High","Low"))
#'
#' metadata_abundance = SummarizedExperiment::colData(sce)[,c(sample_id, group_id, celltype_id)]
#' colnames(metadata_abundance) =c("sample_id", "group_id", "celltype_id")
#' abundance_data = metadata_abundance %>% tibble::as_tibble() %>% dplyr::group_by(sample_id , celltype_id) %>% dplyr::count() %>% dplyr::inner_join(metadata_abundance %>% tibble::as_tibble() %>% dplyr::distinct(sample_id , group_id ))
#' abundance_data = abundance_data %>% dplyr::mutate(keep = n >= min_cells) %>% dplyr::mutate(keep = factor(keep, levels = c(TRUE,FALSE)))
#' abundance_data_receiver = process_info_to_ic(abund_data = abundance_data, ic_type = "receiver")
#' abundance_data_sender = process_info_to_ic(abund_data = abundance_data, ic_type = "sender")
#'
#' celltype_info = get_avg_frac_exprs_abund(sce = sce, sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id)
#'
#' receiver_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "receiver", lr_network = lr_network)
#' sender_info_ic = process_info_to_ic(info_object = celltype_info, ic_type = "sender", lr_network = lr_network)
#' senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
#' sender_receiver_info = combine_sender_receiver_info_ic(sender_info = sender_info_ic,receiver_info = receiver_info_ic,senders_oi = senders_oi,receivers_oi = receivers_oi,lr_network = lr_network)
#'
#' celltype_de = perform_muscat_de_analysis(
#'    sce = sce,
#'    sample_id = sample_id,
#'    celltype_id = celltype_id,
#'    group_id = group_id,
#'    batches = batches,
#'    contrasts = contrasts_oi)
#'
#' sender_receiver_de = combine_sender_receiver_de(
#'  sender_de = celltype_de,
#'  receiver_de = celltype_de,
#'  senders_oi = senders_oi,
#'  receivers_oi = receivers_oi,
#'  lr_network = lr_network)
#'
#' ligand_activities_targets_DEgenes = get_ligand_activities_targets_DEgenes(
#'    receiver_de = celltype_de,
#'    receivers_oi = receivers_oi,
#'    receiver_frq_df_group = celltype_info$frq_df_group,
#'    ligand_target_matrix = ligand_target_matrix)
#'
#'
#' sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)
#' metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()
#' grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% tibble::as_tibble() %>% dplyr::distinct()
#' colnames(grouping_tbl) = c("sample","group")
#'
#' prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 1,"exprs_ligand" = 1,"exprs_receptor" = 1, "frac_exprs_ligand_receptor" = 1,"abund_sender" = 0,"abund_receiver" = 0)
#' frac_cutoff = 0.05
#' prioritization_tables = generate_prioritization_tables(
#'     sender_receiver_info = sender_receiver_info,
#'     sender_receiver_de = sender_receiver_de,
#'     ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
#'     contrast_tbl = contrast_tbl,
#'     sender_receiver_tbl = sender_receiver_tbl,
#'     grouping_tbl = grouping_tbl,
#'     prioritizing_weights = prioritizing_weights,
#'     fraction_cutoff = frac_cutoff, abundance_data_receiver, abundance_data_sender)
#' }
#'
#' @export
#'
#'
#'
#'
generate_prioritization_tables = function(sender_receiver_info, sender_receiver_de, ligand_activities,
                                          prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2)){

  requireNamespace("dplyr")
  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
                                  dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                                                p_val_ligand_adapted = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand_adapted = rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_ligand_adapted, ties.method = "average", na.last = FALSE))) %>%
                                                                  dplyr::arrange(-lfc_pval_ligand)

  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
                                      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                                                    p_val_receptor_adapted = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor_adapted = rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE)/max(rank(p_val_receptor_adapted, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)

  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>% select(test_ligand, pearson, rank) %>% rename(activity = pearson, ligand=test_ligand) %>%
                                      dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                                                    scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)

  # Check if the group exists
  sender_receiver_info_subset <- sender_receiver_info[[last(names(sender_receiver_info))]]
  # Cell-type and condition specificity of expression of ligand:  per ligand: score each sender-condition combination based on expression and fraction
  ligand_celltype_specificity_prioritization = sender_receiver_info_subset %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::select(group, sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
                                               dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)

  # Cell-type and condition specificity of expression of receptor:  per receptor: score each receiver-condition combination based on expression and fraction
  receptor_celltype_specificity_prioritization = sender_receiver_info_subset %>% dplyr::inner_join(sender_receiver_tbl) %>% dplyr::select(group, receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)

  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    #dplyr::inner_join(ligand_activities %>% rename(ligand=test_ligand, activity=pearson) %>% select(ligand, activity) %>% dplyr::distinct()) %>%
    #dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
    dplyr::inner_join(sender_receiver_info_subset) %>%
    dplyr::inner_join(sender_ligand_prioritization) %>%
    dplyr::inner_join(ligand_activity_prioritization) %>%
    dplyr::inner_join(receiver_receptor_prioritization) %>%
    dplyr::inner_join(ligand_celltype_specificity_prioritization) %>%
    dplyr::inner_join(receptor_celltype_specificity_prioritization)

  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = 2*prioritizing_weights["de_ligand"] + 2*prioritizing_weights["de_receptor"] + prioritizing_weights["activity_scaled"] + prioritizing_weights["exprs_ligand"] + prioritizing_weights["exprs_receptor"]
  group_prioritization_tbl = group_prioritization_tbl %>%
    dplyr::mutate(prioritization_score =
                    (
                      (prioritizing_weights["de_ligand"] * scaled_lfc_ligand) +
                        (prioritizing_weights["de_receptor"] * scaled_lfc_receptor) +
                        (prioritizing_weights["de_ligand"] * scaled_p_val_ligand_adapted) +
                        (prioritizing_weights["de_receptor"] * scaled_p_val_receptor_adapted) +
                        (prioritizing_weights["activity_scaled"] * scaled_activity) +
                        (prioritizing_weights["exprs_ligand"] * scaled_avg_exprs_ligand) +
                        (prioritizing_weights["exprs_receptor"] * scaled_avg_exprs_receptor)
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score)

  return (group_prioritization_tbl)

}
