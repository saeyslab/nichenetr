#' @title Normalize values in a vector by quantile scaling and add a pseudovalue of 0.001
#'
#' @description \code{scale_quantile_adapted} Normalize values in a vector by quantile scaling. Add a pseudovalue of 0.001 to avoid having a score of 0 for the lowest value.
#'
#' @usage
#' scale_quantile_adapted(x)
#'
#' @param x A numeric vector.
#'
#' @return A quantile-scaled numeric vector.
#'
#' @examples
#' \dontrun{
#' scale_quantile_adapted(rnorm(5))
#' }
#'
#' @export
#'
scale_quantile_adapted = function(x){
  y = scale_quantile(x,outlier_cutoff = 0)
  y = y + 0.001
  return(y)
}
#' @title Change values in a tibble if some condition is fulfilled.
#'
#' @description \code{mutate_cond} Change values in a tibble if some condition is fulfilled. Credits: https://stackoverflow.com/questions/34096162/dplyr-mutate-replace-several-columns-on-a-subset-of-rows.
#'
#' @usage
#' mutate_cond(.data, condition, ..., envir = parent.frame())
#'
#' @param .data Data frame / tibble
#' @param condition A condition that need to be fulfilled.
#' @param ... The change that need to happen if condition fulfilled -- through the use of `dplyr::mutate()`
#' @param envir parent.frame() by default
#'
#' @return A tibble
#'
#' @examples
#' \dontrun{
#' mutate_cond(df, a == 3, b = 4)
#' }
#'
#' @export
#'
mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% dplyr::mutate(...)
  .data
}
#' @title Calculate differential expression of cell types in one niche versus all other niches of interest.
#'
#' @description \code{calculate_niche_de} Calculate differential expression of cell types in one niche versus all other niches of interest. This is possible for sender cell types and receiver cell types.
#'
#' @usage
#' calculate_niche_de(seurat_obj, niches, type, assay_oi = "SCT")
#'
#' @param seurat_obj Seurat object
#' @param niches a list of lists/niches giving the name, senders and receiver celltypes for each nice. Sender and receiver cell types should be part of Idents(seurat_obj).
#' @param type For what type of cellype is the DE analysis: "sender" or "receiver"?
#' @param assay_oi Which assay need to be used for DE calculation via `FindMarkers`. Default SCT, alternatives: RNA.
#'
#' @return A tibble containing the DE results of the niches versus each other.
#'
#' @examples
#' \dontrun{
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
calculate_niche_de = function(seurat_obj, niches, type, assay_oi = "SCT"){

  if (type == "sender"){
      sender_vs_sender_tbl = NULL
    #  Determine DE between sender cell types across the different niches
      for(niche_n in seq(length(niches))){
        niche = niches[[niche_n]]
        senders_niche = niche$sender %>% unlist() %>% unique()
        senders_other = niches %>% sapply(function(niche){niche$sender})  %>% unlist() %>% unique() %>% setdiff(senders_niche)
        for(sender_niche_n in seq(length(senders_niche))){
          sender_oi = senders_niche[sender_niche_n]
          for(sender_other_n in seq(length(senders_other))){
            sender_other_niche_oi = senders_other[sender_other_n]
            if(is.null(sender_vs_sender_tbl)){
              sender_vs_sender_tbl = tibble(sender = sender_oi, sender_other_niche = sender_other_niche_oi)
            } else {
                if (nrow(sender_vs_sender_tbl %>% filter(sender == sender_other_niche_oi & sender_other_niche == sender_oi)) != 1){
                  sender_vs_sender_tbl = bind_rows(sender_vs_sender_tbl, tibble(sender = sender_oi, sender_other_niche = sender_other_niche_oi))
              }
            }
          }
        }
      }
    DE_sender = niches %>% lapply(function(niche, seurat_obj) {

      senders_niche = niche$sender %>% unlist() %>% unique()
      senders_other = niches %>% sapply(function(niche){niche$sender})  %>% unlist() %>% unique() %>% setdiff(senders_niche)

      DE_sender = senders_niche %>% lapply(function(sender_oi, seurat_obj, senders_other) {

        print(paste0("Calculate Sender DE between: ",sender_oi, " and ", senders_other))

        DE_subtable = senders_other %>% lapply(function(sender_other_niche_oi, seurat_obj, sender_oi) {

          if (nrow(sender_vs_sender_tbl %>% filter(sender == sender_other_niche_oi & sender_other_niche == sender_oi)) == 1){
              DE_sender_oi = NULL
            } else {
              DE_sender_oi = FindMarkers(object = seurat_obj, ident.1 = sender_oi, ident.2 = sender_other_niche_oi, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
              DE_sender_oi = DE_sender_oi %>% mutate(sender = sender_oi, sender_other_niche = sender_other_niche_oi) %>% arrange(-avg_log2FC)
            }
        }, seurat_obj, sender_oi) %>% bind_rows()

      }, seurat_obj, senders_other) %>% bind_rows()

    }, seurat_obj) %>% bind_rows()

    DE_sender_reverse = DE_sender %>% mutate(avg_log2FC = avg_log2FC * -1) %>%
      rename(pct.1_old = pct.1, pct.2_old = pct.2, sender_old = sender, sender_other_niche_old = sender_other_niche) %>%
      rename(pct.1 = pct.2_old, pct.2 = pct.1_old, sender = sender_other_niche_old, sender_other_niche = sender_old) %>%
      select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, sender, sender_other_niche)
    DE_sender = bind_rows(DE_sender, DE_sender_reverse) %>% distinct()
    return(DE_sender)

  }

  if (type == "receiver"){
    receiver_vs_receiver_tbl = NULL
    #  Determine DE between receiver cell types across the different niches
    for(niche_n in seq(length(niches))){
      niche = niches[[niche_n]]
      receivers_niche = niche$receiver %>% unlist() %>% unique()
      receivers_other = niches %>% sapply(function(niche){niche$receiver})  %>% unlist() %>% unique() %>% setdiff(receivers_niche)
      for(receiver_niche_n in seq(length(receivers_niche))){
        receiver_oi = receivers_niche[receiver_niche_n]
        for(receiver_other_n in seq(length(receivers_other))){
          receiver_other_niche_oi = receivers_other[receiver_other_n]
          if(is.null(receiver_vs_receiver_tbl)){
            receiver_vs_receiver_tbl = tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi)
          } else {
            if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == receiver_other_niche_oi & receiver_other_niche == receiver_oi)) != 1){
              receiver_vs_receiver_tbl = bind_rows(receiver_vs_receiver_tbl, tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi))
            }
          }
        }
      }
    }
    print(receiver_vs_receiver_tbl)
    DE_receiver = niches %>% lapply(function(niche, seurat_obj) {

      receivers_niche = niche$receiver %>% unlist() %>% unique()
      receivers_other = niches %>% sapply(function(niche){niche$receiver})  %>% unlist() %>% unique() %>% setdiff(receivers_niche)

      DE_receiver = receivers_niche %>% lapply(function(receiver_oi, seurat_obj, receivers_other) {

        print(paste0("Calculate receiver DE between: ",receiver_oi, " and ", receivers_other))

        DE_subtable = receivers_other %>% lapply(function(receiver_other_niche_oi, seurat_obj, receiver_oi) {

          if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == receiver_other_niche_oi & receiver_other_niche == receiver_oi)) == 1){
            DE_receiver_oi = NULL
          } else {
            DE_receiver_oi = FindMarkers(object = seurat_obj, ident.1 = receiver_oi, ident.2 = receiver_other_niche_oi, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
            DE_receiver_oi = DE_receiver_oi %>% mutate(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi) %>% arrange(-avg_log2FC)
          }
        }, seurat_obj, receiver_oi) %>% bind_rows()

      }, seurat_obj, receivers_other) %>% bind_rows()

    }, seurat_obj) %>% bind_rows()

    DE_receiver_reverse = DE_receiver %>% mutate(avg_log2FC = avg_log2FC * -1) %>%
      rename(pct.1_old = pct.1, pct.2_old = pct.2, receiver_old = receiver, receiver_other_niche_old = receiver_other_niche) %>%
      rename(pct.1 = pct.2_old, pct.2 = pct.1_old, receiver = receiver_other_niche_old, receiver_other_niche = receiver_old) %>%
      select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, receiver, receiver_other_niche)
    DE_receiver = bind_rows(DE_receiver, DE_receiver_reverse) %>% distinct()
    return(DE_receiver)

  }


}
#' @title Calculate differential expression of receiver cell type in one niche versus all other niches of interest: focus on finding DE genes
#'
#' @description \code{calculate_niche_de_targets} Calculate differential expression of receiver cell type in one niche versus all other niches of interest: focus on finding DE genes
#'
#' @usage
#' calculate_niche_de_targets(seurat_obj, niches, expression_pct, ltf_cutoff, assay_oi = "SCT")
#'
#' @inheritParams calculate_niche_de
#' @param expression_pct input of `min.pct` of `Seurat::FindMarkers`
#' @param lfc_cutoff input of `logfc.threshold` of `Seurat::FindMarkers`
#'
#' @return A tibble containing the DE results of the niches versus each other.
#'
#' @examples
#' \dontrun{
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
#' calculate_niche_de_targets(seurat_obj, niches, expression_pct = 0.10, lfc_cutoff = 0.15)
#' }
#'
#' @export
#'
calculate_niche_de_targets = function(seurat_obj, niches, expression_pct, lfc_cutoff, assay_oi = "SCT"){

    receiver_vs_receiver_tbl = NULL
    all_gene_celltype_tbl = NULL
    #  Determine DE between receiver cell types across the different niches
    for(niche_n in seq(length(niches))){
      niche = niches[[niche_n]]
      receivers_niche = niche$receiver %>% unlist() %>% unique()
      receivers_other = niches %>% sapply(function(niche){niche$receiver})  %>% unlist() %>% unique() %>% setdiff(receivers_niche)
      for(receiver_niche_n in seq(length(receivers_niche))){
        receiver_oi = receivers_niche[receiver_niche_n]
        for(receiver_other_n in seq(length(receivers_other))){
          receiver_other_niche_oi = receivers_other[receiver_other_n]
          if(is.null(receiver_vs_receiver_tbl)){
            receiver_vs_receiver_tbl = tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi)
            all_gene_celltype_tbl = tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi, gene = rownames(seurat_obj))
          } else {
            if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == receiver_other_niche_oi & receiver_other_niche == receiver_oi)) != 1){
              receiver_vs_receiver_tbl = bind_rows(receiver_vs_receiver_tbl, tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi))
              new_gene_tbl = tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi, gene = rownames(seurat_obj))
              all_gene_celltype_tbl = bind_rows(all_gene_celltype_tbl, new_gene_tbl)

            }
          }
        }
      }
    }

    DE_receiver = niches %>% lapply(function(niche, seurat_obj) {

      receivers_niche = niche$receiver %>% unlist() %>% unique()
      receivers_other = niches %>% sapply(function(niche){niche$receiver})  %>% unlist() %>% unique() %>% setdiff(receivers_niche)

      DE_receiver = receivers_niche %>% lapply(function(receiver_oi, seurat_obj, receivers_other) {

        print(paste0("Calculate receiver DE between: ",receiver_oi, " and ", receivers_other))

        DE_subtable = receivers_other %>% lapply(function(receiver_other_niche_oi, seurat_obj, receiver_oi) {

          if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == receiver_other_niche_oi & receiver_other_niche == receiver_oi)) == 1){
            DE_receiver_oi = NULL
          } else {
            DE_receiver_oi = FindMarkers(object = seurat_obj, ident.1 = receiver_oi, ident.2 = receiver_other_niche_oi, min.pct = expression_pct, logfc.threshold = lfc_cutoff, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
            DE_receiver_oi = DE_receiver_oi %>% mutate(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi) %>% arrange(-avg_log2FC)
          }
        }, seurat_obj, receiver_oi) %>% bind_rows()

      }, seurat_obj, receivers_other) %>% bind_rows()

    }, seurat_obj) %>% bind_rows()

    DE_receiver_reverse = DE_receiver %>% mutate(avg_log2FC = avg_log2FC * -1) %>%
      rename(pct.1_old = pct.1, pct.2_old = pct.2, receiver_old = receiver, receiver_other_niche_old = receiver_other_niche) %>%
      rename(pct.1 = pct.2_old, pct.2 = pct.1_old, receiver = receiver_other_niche_old, receiver_other_niche = receiver_old) %>%
      select(gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, receiver, receiver_other_niche)
    DE_receiver = bind_rows(DE_receiver, DE_receiver_reverse) %>% distinct()


    all_gene_celltype_tbl_reverse = all_gene_celltype_tbl %>%
      rename(receiver_old = receiver, receiver_other_niche_old = receiver_other_niche) %>%
      rename(receiver = receiver_other_niche_old, receiver_other_niche = receiver_old) %>%
      select(receiver, receiver_other_niche, gene)
    all_gene_celltype_tbl = bind_rows(all_gene_celltype_tbl, all_gene_celltype_tbl_reverse) %>% distinct()
    DE_receiver = all_gene_celltype_tbl %>% left_join(DE_receiver, by = c("receiver", "receiver_other_niche", "gene"))
    DE_receiver = DE_receiver %>% mutate_cond(is.na(p_val_adj), p_val_adj = 1, p_val = 1, avg_log2FC = 0) %>% distinct()
    return(DE_receiver)

}
#' @title Process the DE output of `calculate_niche_de`
#'
#' @description \code{process_niche_de} Process the DE output of `calculate_niche_de`: define what the average and minimum logFC value is after comparing a celltype vs all the celltypes of the other niches.
#'
#' @usage
#' process_niche_de(DE_table, niches, type, expression_pct)
#'
#' @param DE_table Output of `calculate_niche_de`
#' @param expression_pct Percentage of cells of a cell type having a non-zero expression value for a gene such that a gene can be considered expressed by that cell type.
#' @inheritParams calculate_niche_de
#'
#' @return A tibble containing processed DE information
#'
#' @examples
#' \dontrun{
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
#' DE_table = calculate_niche_de(seurat_obj, niches, "sender")
#' process_niche_de(DE_table, niches, "sender",expression_pct = 0.10)
#' }
#'
#' @export
#'
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
#' @title Combine the differential expression information of ligands in the sender celltypes with the differential expression information of their cognate receptors in the receiver cell types
#'
#' @description \code{combine_sender_receiver_de} Combine the differential expression information of ligands in the sender celltypes with the differential expression information of their cognate receptors in the receiver cell types.
#'
#' @usage
#' combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = "min_lfc")
#'
#' @param DE_sender_processed Output of `process_niche_de` with `type = receiver`
#' @param DE_receiver_processed Output of `process_niche_de` with `type = receiver`
#' @param lr_network Ligand-Receptor Network in tibble format: ligand, receptor, bonafide as columns
#' @param specificity_score Defines which score will be used to prioritze ligand-receptor pairs and consider their differential expression. Default and recommended: "min_lfc".
#' "min_lfc" looks at the minimal logFC of the ligand/receptor between the celltype of interest and all the other celltypes.
#' Alternatives: "mean_lfc", "min_score", and "mean_score". Mean uses the average/mean instead of minimum.
#' score is the product of the logFC and the ratio of fraction of expressing cells.
#'
#' @return A tibble giving the differential expression information of ligands in the sender celltypes and their cognate receptors in the receiver cell types.
#'
#' @examples
#' \dontrun{
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
#' DE_sender = calculate_niche_de(seurat_obj, niches, "sender")
#' DE_receiver = calculate_niche_de(seurat_obj, niches, "receiver")
#' expression_pct = 0.10
#' DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
#' DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
#' specificity_score_LR_pairs = "min_lfc"
#' DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
#' }
#'
#' @export
#'
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
#' @title Processing differential expression output of the receiver cell types
#'
#' @description \code{process_receiver_target_de} Processing differential expression output of the receiver cell types -- used before ligand-activity and ligand-target inference.
#'
#' @usage
#' process_receiver_target_de(DE_receiver_targets, niches, expression_pct, specificity_score = "min_lfc")
#'
#' @param DE_receiver_targets Output of `calculate_niche_de` with `type = receiver`
#' @inheritParams process_niche_de
#' @inheritParams combine_sender_receiver_de
#'
#' @return A tibble containing the processed DE information of the receiver cell types -- used before ligand-activity and ligand-target inference.
#'
#' @examples
#' \dontrun{
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
#' DE_receiver = calculate_niche_de(seurat_obj, niches, "receiver")
#' expression_pct = 0.10
#' DE_receiver_processed = process_receiver_target_de(DE_receiver_targets = DE_receiver, niches = niches, expression_pct = expression_pct)
#' }
#'
#' @export
#'
process_receiver_target_de = function(DE_receiver_targets, niches, expression_pct, specificity_score = "min_lfc"){

  DE_receiver = DE_receiver_targets %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2))
  DE_receiver = DE_receiver %>% mutate_cond(is.na(present), present = FALSE) %>% mutate_cond(is.na(diff) | is.nan(diff), diff = 1) %>% mutate(score = diff*avg_log2FC) %>% arrange(-score)

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
#' @title Calculate the ligand activities and infer ligand-target links based on a list of niche-specific genes per receiver cell type
#'
#' @description \code{get_ligand_activities_targets} Calculate the ligand activities and infer ligand-target links based on a list of niche-specific genes per receiver cell type.
#'
#' @usage
#' get_ligand_activities_targets(niche_geneset_list, ligand_target_matrix, top_n_target)
#'
#' @param niche_geneset_list List of lists/niches giving the geneset of interest for the receiver cell type in each niche.
#' @inheritParams nichenet_seuratobj_aggregate
#' @param top_n_target To predict active, affected targets of the prioritized ligands, consider only DE genes if they also belong to the a priori top n ("top_n_targets") targets of a ligand. Default = 200.
#'
#' @return A tibble of ligands, their activities and targets in each receiver cell type
#'
#' @examples
#' \dontrun{
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
# )
#' DE_receiver = calculate_niche_de(seurat_obj, niches, "receiver")
#' expression_pct = 0.10
#' lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff.
#' specificity_score_targets = "min_lfc"
#' DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
#' background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
#' geneset_KC = DE_receiver_processed_targets %>% filter(receiver == niches$KC_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
#' geneset_MoMac2 = DE_receiver_processed_targets %>% filter(receiver == niches$MoMac2_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
#' geneset_MoMac1 = DE_receiver_processed_targets %>% filter(receiver == niches$MoMac1_niche$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
#'top_n_target = 250
#'niche_geneset_list = list(
#'  "KC_niche" = list(
#'    "receiver" = "KCs",
#'    "geneset" = geneset_KC,
#'    "background" = background),
#'  "MoMac1_niche" = list(
#'    "receiver" = "MoMac1",
#'    "geneset" = geneset_MoMac1 ,
#'    "background" = background),
#'  "MoMac2_niche" = list(
#'    "receiver" = "MoMac2",
#'    "geneset" = geneset_MoMac2 ,
#'    "background" = background)
#')
#'ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
#' }
#'
#' @export
#'
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
#' @title Calculate differential expression between spatially different subpopulations of the same cell type
#'
#' @description \code{calculate_spatial_DE} Calculate differential expression between spatially different subpopulations of the same cell type
#'
#' @usage
#' calculate_spatial_DE(seurat_obj, spatial_info, assay_oi = "SCT")
#'
#' @param seurat_obj Seurat object
#' @param spatial_info Tibble giving information about which celltypes should be compared to each other for defining spatial differential expression. Contains the columns "celltype_region_oi", "celltype_other_region", "niche", "celltype_type".
#' @param assay_oi Assay for the DE analysis: RNA, SCT, ...
#'
#' @return A tibble with DE output
#'
#' @examples
#' \dontrun{
#' seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
#' spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
#'celltype_other_region = c("LSECs_central","Hepatocytes_central","Stellate cells_central")
#') %>%
#'  mutate(niche =  "KC_niche", celltype_type = "sender")
#' calculate_spatial_DE(seurat_obj, spatial_info)
#' }
#'
#' @export
#'
calculate_spatial_DE = function(seurat_obj, spatial_info, assay_oi = "SCT"){
  spatial_info$celltype_region_oi %>% lapply(function(celltype_oi, seurat_obj, spatial_info){
    other_celltype = spatial_info %>% filter(celltype_region_oi == celltype_oi) %>% pull(celltype_other_region) %>% unique()
    niche_oi = spatial_info %>% filter(celltype_region_oi == celltype_oi) %>% pull(niche) %>% unique()

    print(paste0("Calculate Spatial DE between: ",celltype_oi, " and ", other_celltype))

    DE_table = FindMarkers(object = seurat_obj, ident.1 = celltype_oi, ident.2 = other_celltype, min.pct = 0, logfc.threshold = 0, only.pos = FALSE, assay = assay_oi) %>% rownames_to_column("gene") %>% as_tibble()
    DE_table = DE_table %>% mutate(celltype = celltype_oi, niche = niche_oi) %>% arrange(-avg_log2FC)
  }, seurat_obj, spatial_info) %>% bind_rows()
}
#' @title Process the spatialDE output
#'
#' @description \code{process_spatial_de} Process the spatialDE output
#'
#' @usage
#' process_spatial_de(DE_table, type, lr_network, expression_pct, specificity_score = "lfc")
#'
#' @param DE_table Output of `calculate_spatial_DE`
#' @inheritParams process_niche_de
#' @inheritParams combine_sender_receiver_de
#'
#' @return A tibble of processed spatial DE information
#'
#' @examples
#' \dontrun{
#' seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
#' spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
#'celltype_other_region = c("LSECs_central","Hepatocytes_central","Stellate cells_central")
#') %>%
#'  mutate(niche =  "KC_niche", celltype_type = "sender")
#' DE_table= calculate_spatial_DE(seurat_obj, spatial_info)
#' processed_spatialDE = process_spatial_de(DE_table, type = "sender", lr_network, expression_pct = 0.10, specificity_score = "lfc")
#' }
#'
#' @export
#'
process_spatial_de = function(DE_table, type, lr_network, expression_pct, specificity_score = "lfc"){
  ### process the DE_tables of senders and receiver

  if(type == "sender"){
    DE_sender = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC)
    if(specificity_score == "lfc"){
      DE_sender = DE_sender %>% rename(ligand = gene, ligand_score_spatial = avg_log2FC, sender = celltype)
    }
    if(specificity_score == "score"){
      DE_sender = DE_sender %>% rename(ligand = gene, ligand_score_spatial = score, sender = celltype)
    }
    DE_sender = DE_sender %>% select(niche, sender, ligand, ligand_score_spatial)  %>% arrange(-ligand_score_spatial) %>% filter(ligand %in% lr_network$ligand)
    return(DE_sender)
  }
  if(type == "receiver"){
    DE_receiver = DE_table %>% mutate(significant = p_val_adj <= 0.05, present = pct.1 >= expression_pct) %>% mutate(pct.1 = pct.1+0.0001, pct.2 = pct.2 + 0.0001) %>% mutate(diff = (pct.1/pct.2)) %>% mutate(score = diff*avg_log2FC)
    if(specificity_score == "lfc"){
      DE_receiver = DE_receiver %>% rename(receptor = gene, receptor_score_spatial = avg_log2FC, receiver = celltype)
    }
    if(specificity_score == "score"){
      DE_receiver = DE_receiver %>% rename(receptor = gene, receptor_score_spatial = score, receiver = celltype)
    }
    DE_receiver = DE_receiver %>% select(niche, receiver, receptor, receptor_score_spatial)  %>% arrange(-receptor_score_spatial) %>% filter(receptor %in% lr_network$receptor)
    return(DE_receiver)
  }
}
#' @title Makes a table similar to the output of `calculate_spatial_DE` and `process_spatial_de`, but now in case you don't have spatial information for the sender and/or receiver celltype. This is needed for comparability reasons.
#'
#' @description \code{get_non_spatial_de} Makes a table similar to the output of `calculate_spatial_DE` and `process_spatial_de`, but now in case you don't have spatial information for the sender and/or receiver celltype. This is needed for comparability reasons.
#'
#' @usage
#' get_non_spatial_de(niches, spatial_info, type, lr_network)
#'
#' @inheritParams calculate_spatial_DE
#' @inheritParams calculate_niche_de
#' @inheritParams process_niche_de
#' @inheritParams combine_sender_receiver_de
#'
#' @return A tibble of mock processed spatial DE information in case you don't have spatial information for the sender and/or receiver celltype.
#'
#' @examples
#' \dontrun{
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
# )
#' seurat_obj = readRDS(url("https://zenodo.org/record/5840787/files/seurat_obj_subset_integrated_zonation.rds"))
#' spatial_info = tibble(celltype_region_oi = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
#'celltype_other_region = c("LSECs_central","Hepatocytes_central","Stellate cells_central")
#') %>%
#'  mutate(niche =  "KC_niche", celltype_type = "sender")
#' get_non_spatial_de(niches, spatial_info, type = "receiver", lr_network)
#' }
#'
#' @export
#'
get_non_spatial_de = function(niches, spatial_info, type, lr_network){

  if(type == "sender"){

    niche_df = niches %>% names() %>% lapply(function(niche_oi, niches){
      tibble(niche = niche_oi, sender = niches[[niche_oi]]$sender)
    }, niches) %>% bind_rows()
    spatial_df = niches %>% purrr::map("sender") %>% unlist() %>% unique() %>% setdiff(spatial_info %>% filter(celltype_type == "sender") %>% pull(celltype_region_oi)) %>% lapply(function(sender_oi, niches, lr_network){
      spatial_df =  tibble(ligand = lr_network$ligand %>% unique(), ligand_score_spatial = 0)
      spatial_df = spatial_df %>% mutate(sender = sender_oi)
    },  niches, lr_network) %>% bind_rows() %>% select(sender, ligand, ligand_score_spatial)  %>% inner_join(niche_df, by = "sender")

  }
  if(type == "receiver"){

    niche_df = niches %>% names() %>% lapply(function(niche_oi, niches){
      tibble(niche = niche_oi, receiver = niches[[niche_oi]]$receiver)
    }, niches) %>% bind_rows()
    spatial_df = niches %>% purrr::map("receiver") %>% unlist() %>% unique() %>% setdiff(spatial_info %>% filter(celltype_type == "receiver") %>% pull(celltype_region_oi)) %>% lapply(function(receiver_oi, niches, lr_network){
      spatial_df =  tibble(receptor = lr_network$receptor %>% unique(), receptor_score_spatial = 0)
      spatial_df = spatial_df %>% mutate(receiver = receiver_oi)
    },  niches, lr_network) %>% bind_rows() %>% select(receiver, receptor, receptor_score_spatial)  %>% inner_join(niche_df, by = "receiver")


  }
  return(spatial_df)
}
#' @title Use the information from the niche- and spatial differential expression analysis of ligand-senders and receptor-receivers pairs, in addition to the ligand activity prediction and ligand-target inferernce, in order to make a final ligand-receptor and ligand-target prioritization table.
#'
#' @description \code{get_prioritization_tables} Use the information from the niche- and spatial differential expression analysis of ligand-senders and receptor-receivers pairs, in addition to the ligand activity prediction and ligand-target inferernce, in order to make a final ligand-receptor and ligand-target prioritization table.
#'
#' @usage
#' get_prioritization_tables(output_nichenet_analysis, prioritizing_weights)
#'
#' @param output_nichenet_analysis List containing following data frames: DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed,
#' receiver_spatial_DE_processed, ligand_activities_targets, DE_receiver_processed_targets, exprs_tbl_ligand, exprs_tbl_receptor, exprs_tbl_target
#' @param prioritizing_weights Named numeric vector in the form of:
#' #' prioritizing_weights = c("scaled_ligand_score" = 5,
# "scaled_ligand_expression_scaled" = 1,
# "ligand_fraction" = 1,
# "scaled_ligand_score_spatial" = 2,
# "scaled_receptor_score" = 0.5,
# "scaled_receptor_expression_scaled" = 0.5,
# "receptor_fraction" = 1,
# "ligand_scaled_receptor_expression_fraction" = 1,
# "scaled_receptor_score_spatial" = 0,
# "scaled_activity" = 0,
# "scaled_activity_normalized" = 1,
# "bona_fide" = 1)
#'
#' @return A list containing a prioritization table for ligand-receptor interactions, and one for ligand-target interactions
#'
#' @examples
#' \dontrun{
#' prioritizing_weights = c("scaled_ligand_score" = 5,
# "scaled_ligand_expression_scaled" = 1,
# "ligand_fraction" = 1,
# "scaled_ligand_score_spatial" = 2,
# "scaled_receptor_score" = 0.5,
# "scaled_receptor_expression_scaled" = 0.5,
# "receptor_fraction" = 1,
# "ligand_scaled_receptor_expression_fraction" = 1,
# "scaled_receptor_score_spatial" = 0,
# "scaled_activity" = 0,
# "scaled_activity_normalized" = 1,
# "bona_fide" = 1)
#' output_nichenet_analysis = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
# ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
#' prioritization_tables = get_prioritization_tables(output_nichenet_analysis, prioritizing_weights)
#' }
#'
#' @export
#'
get_prioritization_tables = function(output_nichenet_analysis, prioritizing_weights){
  combined_information = output_nichenet_analysis$DE_sender_receiver %>%
    inner_join(output_nichenet_analysis$ligand_scaled_receptor_expression_fraction_df, by = c("receiver", "ligand", "receptor")) %>%
    inner_join(output_nichenet_analysis$sender_spatial_DE_processed, by = c("niche", "sender", "ligand")) %>%
    inner_join(output_nichenet_analysis$receiver_spatial_DE_processed, by = c("niche", "receiver", "receptor")) %>%
    inner_join(output_nichenet_analysis$ligand_activities_targets, by = c("receiver", "ligand")) %>%
    left_join(output_nichenet_analysis$DE_receiver_processed_targets, by = c("niche", "receiver", "target")) %>% ## if ligand has no target genes --> it gets NA as target value --> these ligands should not be removed --> therefore left_join instead of inner_join
    inner_join(output_nichenet_analysis$exprs_tbl_ligand, by = c("sender", "ligand")) %>%
    inner_join(output_nichenet_analysis$exprs_tbl_receptor, by = c("receiver", "receptor"))  %>%
    left_join(output_nichenet_analysis$exprs_tbl_target, by = c("receiver", "target"))  ## if ligand has no target genes --> it gets NA as target value --> these ligands should not be removed --> therefore left_join instead of inner_join

  # reorder the columns

  combined_information = combined_information %>% mutate(ligand_receptor = paste(ligand, receptor, sep = "--"))  %>%  mutate(bonafide_score = 1) %>%  mutate_cond(bonafide == FALSE, bonafide_score = 0.5)

  combined_information = combined_information %>% select(
    niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
    ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_spatial,
    receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_spatial,
    ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor, bonafide_score,
    target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
    activity, activity_normalized,
    scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
    scaled_ligand_score_spatial, scaled_receptor_score_spatial,
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
                       (prioritizing_weights["scaled_ligand_score_spatial"] * scaled_ligand_score_spatial) +
                       (prioritizing_weights["scaled_receptor_score_spatial"] * scaled_receptor_score_spatial) +
                       (prioritizing_weights["ligand_scaled_receptor_expression_fraction"] * ligand_scaled_receptor_expression_fraction) +
                       (prioritizing_weights["scaled_activity"] * scaled_activity) +
                       (prioritizing_weights["scaled_activity_normalized"] * scaled_activity_normalized) +
                       (prioritizing_weights["ligand_fraction"] * scaled_ligand_fraction_adapted ) +
                       (prioritizing_weights["receptor_fraction"] * scaled_receptor_fraction_adapted  ) +
                       (prioritizing_weights["bona_fide"] * bonafide_score)
                    )* (1/length(prioritizing_weights))) %>% dplyr::arrange(-prioritization_score)

  prioritization_tbl_ligand_receptor = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide,
                                                                                   ligand_score,ligand_significant, ligand_present, ligand_expression, ligand_expression_scaled, ligand_fraction, ligand_score_spatial,
                                                                                   receptor_score, receptor_significant, receptor_present, receptor_expression, receptor_expression_scaled, receptor_fraction, receptor_score_spatial,
                                                                                   ligand_scaled_receptor_expression_fraction, avg_score_ligand_receptor,
                                                                                   activity, activity_normalized,
                                                                                   scaled_ligand_score, scaled_ligand_expression_scaled, scaled_receptor_score, scaled_receptor_expression_scaled, scaled_avg_score_ligand_receptor,
                                                                                   scaled_ligand_score_spatial, scaled_receptor_score_spatial,
                                                                                   scaled_ligand_fraction_adapted, scaled_receptor_fraction_adapted,
                                                                                   scaled_activity, scaled_activity_normalized,
                                                                                   prioritization_score) %>% distinct()
  prioritization_tbl_ligand_target = combined_information_prioritized %>% select(niche, receiver, sender, ligand_receptor, ligand, receptor, bonafide, target,
                                                                                 target_score, target_significant, target_present, target_expression, target_expression_scaled, target_fraction, ligand_target_weight,
                                                                                 activity, activity_normalized, scaled_activity, scaled_activity_normalized, prioritization_score) %>% distinct()

  return(list(prioritization_tbl_ligand_receptor = prioritization_tbl_ligand_receptor, prioritization_tbl_ligand_target = prioritization_tbl_ligand_target))

}
