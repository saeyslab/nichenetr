scale_quantile_adapted = function(x, outlier_cutoff = 0){
  y = nichenetr::scale_quantile(x,outlier_cutoff  = outlier_cutoff)
  y = y + 0.001
  return(y)
}

check_names <- function(column, seurat_obj = NA){
  # If seurat_obj is NA
  if (typeof(seurat_obj) == "logical") {
    if (column != make.names(column)) {
      stop(paste0("'", column, "' is not a syntactically valid R name - check make.names"))
    }
  } else {
    if (!all(unique(seurat_obj[[column, drop=TRUE]]) == make.names(unique(seurat_obj[[column, drop=TRUE]])))){
      stop(paste0("'", column, "' column should have syntactically valid R names - see make.names"))
    }
  }
}

#' @title Calculate differential expression of one cell type versus all other cell types
#'
#' @description \code{calculate_de} Calculate differential expression of one cell type versus all other cell types. If condition_oi is provided, only consider cells from that condition.
#'
#' @usage
#' calculate_de(seurat_obj, celltype_id, condition_oi = NA, group_id = NA, min.pct = 0.1, assay_oi = "RNA")
#'
#' @param seurat_obj Seurat object
#' @param celltype_id Name of the meta data column that indicates the cell type of a cell
#' @param condition_oi If provided, subset seurat_obj so DE is only calculated for cells belonging to condition_oi
#' @param group_id group_id Name of the meta data column that indicates from which group/condition a cell comes from
#' @param assay_oi Which assay need to be used for DE calculation via `FindMarkers`. Default RNA, alternatives: SCT.
#'
#' @return A dataframe containing the DE results
#'
#' @examples
#' \dontrun{
#' seurat_obj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
#' seurat_obj$celltype <- make.names(seurat_obj$celltype)
#' # Calculate cell-type specific markers across conditions
#' calculate_de(seurat_obj, "celltype")
#' # Calculate LCMV-specific cell-type markers
#' calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#' }
#'
#' @export
#'
calculate_de = function(seurat_obj, celltype_id,
                        condition_oi = NA, group_id = NA,
                        min.pct = 0.1, assay_oi = "RNA"){

  if (any(!is.na(group_id), !is.na(condition_oi)) & !all(!is.na(group_id), !is.na(condition_oi))){
    stop("Please input both group_id and condition_oi")
  }

  # Check names
  sapply(c(celltype_id, condition_oi, group_id) %>% .[!is.na(.)], check_names)
  sapply(celltype_id, check_names, seurat_obj)

  # Subset seurat obj to condition of interest
  if (!is.na(condition_oi)) {
    seurat_obj = seurat_obj[,seurat_obj[[group_id]] == condition_oi]
  }

  DE_table = FindAllMarkers(seurat_obj, min.pct = min.pct) %>%
                rename(cluster_id = cluster)

  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table)
  if(!SeuratV4){
    DE_table = DE_table %>% dplyr::rename(avg_log2FC = avg_logFC)
  }

  return(DE_table)

}
#' @title Calculate average of gene expression per cell type.
#'
#' @description \code{get_exprs_avg}  Calculate average of gene expression per cell type. If condition_oi is provided, only consider cells from that condition.
#' @usage get_exprs_avg(seurat_obj, celltype_id, condition_oi = NA, group_id = NA)
#'
#' @return Data frame with average gene expression per cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom Seurat NormalizeData
#'
#' @examples
#' \dontrun{
#' seurat_obj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
#' seurat_obj$celltype <- make.names(seuratObj$celltype)
#' # Calculate average expression across conditions
#' expression_info = get_exprs_avg(seurat_obj, "celltype")
#' # Calculate LCMV-specific average expression
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#' }
#'
#' @export
#'
get_exprs_avg = function(seurat_obj, celltype_id,
                         condition_oi = NA, group_id = NA){

  requireNamespace("dplyr")

  if (any(!is.na(group_id), !is.na(condition_oi)) & !all(!is.na(group_id), !is.na(condition_oi))){
    stop("Please input both group_id and condition_oi")
  }

  # Check names
  sapply(c(celltype_id, condition_oi, group_id) %>% .[!is.na(.)], check_names)
  sapply(celltype_id, check_names, seurat_obj)

  # Subset seurat object
  if (!is.na(condition_oi)) {
    seurat_obj = seurat_obj[,seurat_obj[[group_id]] == condition_oi]
  }

  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  avg_celltype <- AverageExpression(seurat_obj, assay = "RNA", slot = "data", group.by = celltype_id) %>%
                    .$RNA %>% data.frame() %>% rownames_to_column("gene") %>%
                    pivot_longer(!gene, names_to = "cluster_id", values_to = "avg_expr")

  return (avg_celltype)



}
#' @title Process DE or expression information into intercellular communication focused information.
#'
#' @description \code{process_info_to_ic} First, only keep information of ligands for senders_oi, and information of receptors for receivers_oi.
#' Then, combine information for senders and receivers by linking ligands to receptors based on the prior knowledge ligand-receptor network.
#' @usage process_info_to_ic(info_object, ic_type = "sender", lr_network)
#' @param table_object Output of `get_exprs_avg` or `calculate_de`
#' @param table_type "expression" or "DE": indicates whether the table contains expression or DE information
#' @param lr_network
#' @param senders_oi
#' @param receivers_oi
#' @return Dataframe combining sender and receiver information linked to each other through joining by the ligand-receptor network.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' seurat_obj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
#' seurat_obj$celltype <- make.names(seuratObj$celltype)
#' # Calculate LCMV-specific average expression
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#' # Calculate LCMV-specific cell-type markers
#' DE_table = calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#' processed_expr_info = process_table_to_ic(expression_info, table_type = "expression", lr_network)
#' processed_DE_table <- process_table_to_ic(DE_table, table_type = "DE", lr_network,
#' senders_oi = c("CD4.T", "Treg", "Mono", "NK", "B", "DC"), receivers_oi = "CD8.T")
#' }
#'
#' @export
#'
process_table_to_ic = function(table_object, table_type = "expression",
                               lr_network, senders_oi = NA, receivers_oi = NA){

  ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
  receptors = lr_network %>% dplyr::pull(receptor) %>% unique()

  sender_table <- table_object %>% dplyr::rename(sender = cluster_id, ligand = gene)
  receiver_table <- table_object %>% dplyr::rename(receiver = cluster_id, receptor = gene)

  if (table_type == "expression"){
    if (!any(is.na(senders_oi))) warning("senders_oi is given. The expression data will be scaled with all remaining cell types, so it is recommended that senders_oi = NA")
    if (!any(is.na(receivers_oi))) warning("receivers_oi is given. The expression data will be scaled with all remaining cell types, so it is recommended that receivers_oi = NA")

    sender_table <- sender_table %>% dplyr::rename(avg_ligand = avg_expr)
    receiver_table <- receiver_table %>% dplyr::rename(avg_receptor = avg_expr)
    columns_select <- c("sender", "receiver", "ligand", "receptor", "avg_ligand", "avg_receptor", "ligand_receptor_prod")

  } else if (table_type == "DE"){
    if (any(is.na(senders_oi))) warning("senders_oi is NA. For DE filtering, it is best if this parameter is given.")
    if (any(is.na(receivers_oi))) warning("receivers_oi is NA. For DE filtering, it is best if this parameter is given.")

    sender_table = sender_table %>% dplyr::rename(avg_ligand = avg_log2FC, p_val_ligand = p_val,  p_adj_ligand = p_val_adj)
    receiver_table =  receiver_table %>% dplyr::rename(avg_receptor = avg_log2FC, p_val_receptor = p_val,  p_adj_receptor = p_val_adj)
    columns_select <- c("sender", "receiver", "ligand", "receptor", "lfc_ligand", "lfc_receptor", "ligand_receptor_lfc_avg", "p_val_ligand", "p_adj_ligand", "p_val_receptor", "p_adj_receptor")
  }

  # Filter senders and receivers if it is not NA
  sender_table <- sender_table %>% {if (all(!is.na(senders_oi))) filter(., sender %in% senders_oi) else (.)}
  receiver_table <- receiver_table %>% {if (all(!is.na(receivers_oi))) filter(., receiver %in% receivers_oi) else (.)}

  # Join sender-ligand-receptor-receiver
  sender_receiver_table <- sender_table %>% dplyr::inner_join(lr_network, by = "ligand") %>%
    dplyr::inner_join(receiver_table, by = c("receptor"))

  # Calculate average expression
  sender_receiver_table <- sender_receiver_table %>%
    mutate(ligand_receptor_avg = case_when(
                table_type == "expression" ~ avg_ligand * avg_receptor,
                table_type == "DE" ~ (avg_ligand + avg_ligand)/2
          )
    )  %>% arrange(-ligand_receptor_avg) %>%
    # Rename columns appropriately
    {if (table_type == "expression") rename(., "ligand_receptor_prod" = "ligand_receptor_avg")
      else rename(., "ligand_receptor_lfc_avg" = "ligand_receptor_avg", "lfc_ligand" = "avg_ligand", "lfc_receptor" = "avg_receptor")} %>%
    select(all_of(columns_select)) %>% dplyr::distinct()

  return(sender_receiver_table)

}


#' @title generate_prioritization_tables
#'
#' @description \code{generate_prioritization_tables}  Perform a prioritization of cell-cell interactions (similar to MultiNicheNet).
#' User can choose the importance attached to each of the following prioritization criteria: differential expression of ligand and receptor, cell-type specificity of expression of ligand and receptor, NicheNet ligand activity
#' @usage generate_prioritization_tables(sender_receiver_info, sender_receiver_de, ligand_activities, prioritizing_weights = c("de_ligand" = 1,"de_receptor" = 1,"activity_scaled" = 2,"exprs_ligand" = 2,"exprs_receptor" = 2))
#'
#' @param sender_receiver_info Output of `get_exprs_avg` -> `process_table_to_ic`
#' @param sender_receiver_de Output of`calculate_de` -> `process_table_to_ic`
#' @param ligand_activities Output of `predict_ligand_activities`
#' @param prioritizing_weights Named vector indicating the relative weights of each prioritization criterion
#'
#' @return Data frames of prioritized sender-ligand-receiver-receptor interactions.
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
#' lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% dplyr::distinct(ligand, receptor)
#' ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
#' seurat_obj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
#' seurat_obj$celltype <- make.names(seuratObj$celltype)
#' sender_celltypes = c("CD4.T","Treg", "Mono", "NK", "B", "DC")
#' receiver = "CD8.T"
#'
#' # Convert lr_network from mouse to human
#' lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
#' colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
#' rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
#' ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
#'
#' # Ligand activity analysis
#' seurat_obj_receiver = subset(seurat_obj, idents = receiver) %>% SetIdent(value = .[["aggregate"]])
#' geneset_oi = FindMarkers(object = seurat_obj_receiver, ident.1 = "LCMV, ident.2 = "SS, min.pct = 0.10) %>% rownames_to_column("gene") %>%
#'      filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)]
#' expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.10) %>% unlist() %>% unique()
#' expressed_genes_receiver = get_expressed_genes(receiver, seurat_obj, pct = 0.10)
#' expressed_ligands = intersect(lr_network %>% pull(ligand) %>% unique(), expressed_genes_sender)
#' expressed_receptors = intersect(lr_network %>% pull(receiver) %>% unique(), expressed_genes_receiver)
#' potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>% pull(from) %>% unique()
#' ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)],
#'                                             ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
#'
#' # Calculate LCMV-specific average expression
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#'
#' # Calculate LCMV-specific cell-type markers
#' DE_table = calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", group_id = "aggregate")
#'
#' # Process tables
#' processed_expr_info = process_table_to_ic(expression_info, table_type = "expression", lr_network)
#' processed_DE_table <- process_table_to_ic(DE_table, table_type = "DE", lr_network,
#'                                           senders_oi = sender_celltypes, receivers_oi = receiver)
#' }
#'
#' # Generate prioritization tables
#' prioritizing_weights = c("de_ligand" = 1, "de_receptor" = 1, "activity_scaled" = 2, "exprs_ligand" = 1, "exprs_receptor" = 1)
#' generate_prioritization_tables(processed_expr_info,
#'                                processed_DE_table,
#'                                ligand_activities,
#'                                prioritizing_weights)
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


  # Cell-type specificity of expression of ligand:  per ligand: score each sender combination based on expression
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
                                               dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)

  # Cell-type specificity of expression of receptor:  per receptor: score each receiver combination based on expression
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
    dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)

  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    #dplyr::inner_join(ligand_activities %>% rename(ligand=test_ligand, activity=pearson) %>% select(ligand, activity) %>% dplyr::distinct()) %>%
    #dplyr::mutate(lr_interaction = paste(ligand, receptor, sep = "_")) %>% dplyr::mutate(id = paste(lr_interaction, sender, receiver, sep = "_")) %>%
    dplyr::inner_join(sender_receiver_info) %>%
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
