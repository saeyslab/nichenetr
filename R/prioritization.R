check_names <- function(column, seurat_obj = NULL){
  # If seurat_obj is NA
  if (is.null(seurat_obj)) {
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
#' @description \code{calculate_de} Calculate differential expression of one cell type versus all other cell types using Seurat::FindAllMarkers. If condition_oi is provided, only consider cells from that condition.
#'
#' @usage
#' calculate_de(seurat_obj, celltype_colname, condition_oi = NA, condition_colname = NA, assay_oi = "RNA", ...)
#'
#' @param seurat_obj Seurat object
#' @param celltype_colname Name of the meta data column that indicates the cell type of a cell
#' @param condition_oi If provided, subset seurat_obj so DE is only calculated for cells belonging to condition_oi
#' @param condition_colname Name of the meta data column that indicates from which group/condition a cell comes from
#' @param assay_oi Which assay need to be used for DE calculation. If NULL, will use \code{DefaultAssay}
#' @param ... Arguments passed to Seurat::FindAllMarkers(by default: features = NULL, min.pct = 0, logfc.threshold = 0, return.thresh = 1)
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
#' calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#' }
#'
#' @export
#'
calculate_de = function(seurat_obj, celltype_colname,
                        condition_oi = NULL, condition_colname = NULL,
                        assay_oi = "RNA",
                        ...){

  # Default settings to return all genes with their p-val and LFC
  FindAllMarkers_args = list(assay = assay_oi,
                             features = NULL, min.pct = 0,
                             logfc.threshold = 0,
                             return.thresh = 1)

  # Replace this with user arguments
  FindAllMarkers_args[names(list(...))] =  list(...)

  if (any(!is.null(condition_colname), !is.null(condition_oi)) & !all(!is.null(condition_colname), !is.null(condition_oi))){
    stop("Please input both condition_colname and condition_oi")
  }

  # Subset seurat obj to condition of interest
  if (!is.null(condition_oi)) {
    seurat_obj = seurat_obj[,seurat_obj[[condition_colname]] == condition_oi]
  }

  # Set celltype as identity class
  Idents(seurat_obj) <- seurat_obj[[celltype_colname, drop=TRUE]]

  FindAllMarkers_args$object <- seurat_obj
  DE_table <- do.call(FindAllMarkers, FindAllMarkers_args) %>%
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
#'
#' @inheritParams calculate_de
#' @param condition_oi If provided, subset seurat_obj so average expression is only calculated for cells belonging to condition_oi
#' @param ... Arguments passed to Seurat::AverageExpression, usually for slot/layer to use (default: data)
#'
#' @return Data frame with average gene expression per cell type.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' seurat_obj = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj.rds"))
#' seurat_obj$celltype <- make.names(seuratObj$celltype)
#' # Calculate average expression across conditions
#' expression_info = get_exprs_avg(seurat_obj, "celltype")
#' # Calculate LCMV-specific average expression
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#' }
#'
#' @export
#'
get_exprs_avg = function(seurat_obj, celltype_colname,
                         condition_oi = NULL, condition_colname = NULL,
                         assay_oi = NULL, ...){

  requireNamespace("dplyr")

  # Check that assay_oi is in seurat_obj
  if (!is.null(assay_oi)){
    if (!assay_oi %in% names(seurat_obj@assays)) {
      stop(paste0("assay_oi '", assay_oi, "' does not exist in seurat_obj"))
    }
  # Otherwise, use DefaultAssay
  } else {
    assay_oi <- DefaultAssay(seurat_obj)
  }

  if (any(!is.null(condition_colname), !is.null(condition_oi)) & !all(!is.null(condition_colname), !is.null(condition_oi))){
    stop("Please input both condition_colname and condition_oi")
  }

  # Subset seurat object
  if (!is.null(condition_oi)) {
    seurat_obj <- seurat_obj[,seurat_obj[[condition_colname]] == condition_oi]
  }

  celltypes <- unique(seurat_obj[[celltype_colname, drop=TRUE]])

  avg_celltype <- AverageExpression(seurat_obj, group.by = celltype_colname, assays = assay_oi, ...) %>%
                    .[[assay_oi]] %>% data.frame(check.names=FALSE) %>% rownames_to_column("gene") %>%
                    pivot_longer(!gene, names_to = "cluster_id", values_to = "avg_expr")

  # If any celltypes had an underscore in their name
  if (any(grepl("_", celltypes))){
    # Map the new names and the original names
    # This is so it works in the case the original name also already has a hyphen in in
    mapping <- data.frame(orig_name = sort(celltypes), cluster_id = sort(unique(avg_celltype$cluster_id)))
    avg_celltype <- avg_celltype %>% left_join(mapping, by = "cluster_id") %>%
                    dplyr::select(gene, cluster_id = orig_name, avg_expr)


  }

  return (avg_celltype)

}
#' @title Process DE or expression information into intercellular communication focused information.
#'
#' @description \code{process_table_to_ic} First, only keep information of ligands for senders_oi, and information of receptors for receivers_oi.
#' Then, combine information for senders and receivers by linking ligands to receptors based on the prior knowledge ligand-receptor network.
#' @usage process_table_to_ic(table_object, table_type = "expression", lr_network, senders_oi = NULL, receivers_oi = NULL)
#' @param table_object Output of \code{get_exprs_avg}, \code{calculate_de}, or \code{FindMarkers}
#' @param table_type "expression", "celltype_DE", or "group_DE": indicates whether the table contains expression, celltype markers, or condition-specific information
#' @param lr_network Prior knowledge Ligand-Receptor network (columns: ligand, receptor)
#' @param senders_oi Default NULL: all celltypes will be considered as senders. If you want to select specific senders of interest: you can add this here as character vector.
#' @param receivers_oi Default NULL: all celltypes will be considered as receivers If you want to select specific receivers of interest: you can add this here as character vector.
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
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#' # Calculate LCMV-specific cell-type markers
#' DE_table = calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#' # Calculate LCMV-specific genes across cell types
#' condition_markers <- FindMarkers(object = seuratObj, ident.1 = "LCMV", ident.2 = "SS",
#'                                 group.by = "aggregate", min.pct = 0, logfc.threshold = 0) %>% rownames_to_column("gene")
#' processed_expr_info = process_table_to_ic(expression_info, table_type = "expression", lr_network)
#' processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network,
#' senders_oi = c("CD4.T", "Treg", "Mono", "NK", "B", "DC"), receivers_oi = "CD8.T")
#' processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "condition_DE", lr_network)
#' }
#'
#' @export
#'
process_table_to_ic = function(table_object, table_type = "expression",
                               lr_network, senders_oi = NULL, receivers_oi = NULL){

  # Rename lr_network columns to "ligand", "receptor" if the columns "from","to" exist
  if (all(c("from", "to") %in% colnames(lr_network)) &
      !all(c("ligand", "receptor") %in% colnames(lr_network))){
    lr_network <- lr_network %>% dplyr::rename(ligand = from, receptor = to)
  }

  ligands = lr_network %>% dplyr::pull(ligand) %>% unique()
  receptors = lr_network %>% dplyr::pull(receptor) %>% unique()

  if (table_type == "expression"){
    if (!is.null(senders_oi)) warning("senders_oi is given. The expression data will be scaled with all remaining cell types, so it is recommended that senders_oi = NULL")
    if (!is.null(receivers_oi)) warning("receivers_oi is given. The expression data will be scaled with all remaining cell types, so it is recommended that receivers_oi = NULL")

    sender_table <- table_object %>% dplyr::rename(sender = cluster_id, ligand = gene, avg_ligand = avg_expr)
    receiver_table <- table_object %>% dplyr::rename(receiver = cluster_id, receptor = gene, avg_receptor = avg_expr)
    columns_select <- c("sender", "receiver", "ligand", "receptor", "avg_ligand", "avg_receptor", "ligand_receptor_prod")

  } else if (table_type == "celltype_DE"){
    if (is.null(senders_oi)) warning("senders_oi is NULL For DE filtering, it is best if this parameter is given.")
    if (is.null(receivers_oi)) warning("receivers_oi is NULL For DE filtering, it is best if this parameter is given.")

    sender_table <- table_object %>% dplyr::rename(sender = tidyr::starts_with("cluster"), ligand = gene, avg_ligand = avg_log2FC, p_val_ligand = p_val,  p_adj_ligand = p_val_adj, pct_expressed_sender = pct.1)
    receiver_table <-  table_object %>% dplyr::rename(receiver = tidyr::starts_with("cluster"), receptor = gene, avg_receptor = avg_log2FC, p_val_receptor = p_val, p_adj_receptor = p_val_adj, pct_expressed_receiver = pct.1)
    columns_select <- c("sender", "receiver", "ligand", "receptor", "lfc_ligand", "lfc_receptor", "ligand_receptor_lfc_avg", "p_val_ligand", "p_adj_ligand", "p_val_receptor", "p_adj_receptor", "pct_expressed_sender", "pct_expressed_receiver")

  } else if (table_type == "group_DE") {
    if (!is.null(senders_oi)) stop("senders_oi is given. Since we do not consider cell type specificity with group DE, please change this to NULL")
    if (!is.null(receivers_oi)) stop("receivers_oi is given. Since we do not consider cell type specificity with group DE, please change this to NULL")

    sender_table = table_object %>% dplyr::rename(ligand = gene, avg_ligand = avg_log2FC, p_val_ligand = p_val,  p_adj_ligand = p_val_adj)
    receiver_table = table_object %>% dplyr::rename(receptor = gene, avg_receptor = avg_log2FC, p_val_receptor = p_val, p_adj_receptor = p_val_adj)
    columns_select <- c("ligand", "receptor", "lfc_ligand", "lfc_receptor", "ligand_receptor_lfc_avg", "p_val_ligand", "p_adj_ligand", "p_val_receptor", "p_adj_receptor")

  }

  # Filter senders and receivers if it is not NULL
  sender_table <- sender_table %>% {if (!is.null(senders_oi)) filter(., sender %in% senders_oi) else (.)}
  receiver_table <- receiver_table %>% {if (!is.null(receivers_oi)) filter(., receiver %in% receivers_oi) else (.)}

  # Join sender-ligand-receptor-receiver
  sender_receiver_table <- sender_table %>% dplyr::inner_join(lr_network, by = "ligand") %>%
    dplyr::inner_join(receiver_table, by = "receptor")

  # Calculate average expression
  sender_receiver_table <- sender_receiver_table %>%
    mutate(ligand_receptor_avg = case_when(
                table_type == "expression" ~ avg_ligand * avg_receptor,
                grepl("DE$", table_type) ~ (avg_ligand + avg_receptor)/2
          )
    )  %>% arrange(-ligand_receptor_avg) %>%
    # Rename columns appropriately
    {if (table_type == "expression") rename(., "ligand_receptor_prod" = "ligand_receptor_avg")
      else rename(., "ligand_receptor_lfc_avg" = "ligand_receptor_avg", "lfc_ligand" = "avg_ligand", "lfc_receptor" = "avg_receptor")} %>%
    select(all_of(columns_select)) %>% dplyr::distinct()

  return(sender_receiver_table)

}

#' @title Generate tables used for \code{generate_prioritization_tables}
#' @description Calculate differential expression, average expression, and condition specificity of ligands and receptors.
#' @param seuratObj Seurat object
#' @param lr_network_filtered Ligand-receptor network that has been filtered to only contain ligands and receptors that are expressed
#' @param condition_reference Reference condition for condition specificity calculation
#' @param scenario "case_control" or "one_condition": if "case_control", calculate condition specificity. If "one_condition", only calculate cell type specificity.
#' @param ... Arguments passed to \code{FindAllMarkers}, \code{FindMarkers}, and \code{AverageExpression}
#' @inheritParams calculate_de
#' @inheritParams process_table_to_ic
#' @return List of dataframes containing sender-receiver DE, sender-receiver expression, and condition DE
#' @export
generate_info_tables <- function(seuratObj,
                                 celltype_colname,
                                 senders_oi,
                                 receivers_oi,
                                 lr_network_filtered,
                                 condition_colname = NULL,
                                 condition_oi = NULL,
                                 condition_reference = NULL,
                                 scenario = "case_control",
                                 assay_oi = NULL,
                                 ...){
  requireNamespace("dplyr")

  # Check that celltype_colname exists in Seurat object
  if (!celltype_colname %in% names(seuratObj@meta.data)){
    stop(paste0("celltype_colname '", celltype_colname, "' does not exist in seuratObj"))
  }

  if (any(!is.null(condition_colname), !is.null(condition_oi), !is.null(condition_reference)) &
      !all(!is.null(condition_colname), !is.null(condition_oi), !is.null(condition_reference))){
    stop("condition_* arguments must be either all NULL or all provided")
  }

  # If condition is not provided
  if (all(is.null(condition_colname), is.null(condition_oi), is.null(condition_reference)) &
      scenario == "case_control"){
    stop("condition_* arguments are not provided but the 'case_control' scenario is selected. Please either provide condition_* arguments or change scenario to 'one_condition', as condition specificity cannot be calculated.")
  }

  # If condition is provided but scenario is one_condition
  if (!is.null(condition_colname))
    # Check that condition_colname exists in Seurat object
    if (!condition_colname %in% names(seuratObj@meta.data)){
      stop(paste0("condition_colname '", condition_colname, "' does not exist in seuratObj"))
    }
    if(scenario == "one_condition"){
    warning("condition_* arguments are provided but the 'one_condition' scenario is selected. Only cells from the condition_oi will be used to calculate cell type specificity, condition specificity will not be calculated.")
  }

  # Check that senders_oi is in Seurat object
  if (!all(senders_oi %in% unique(seuratObj[[celltype_colname, drop=TRUE]]))){
    stop(paste0("Not all senders_oi exist in the seurat object"))
  }

  # Check that receivers_oi is in Seurat object
  if (!all(receivers_oi %in% unique(seuratObj[[celltype_colname, drop=TRUE]]))){
    stop(paste0("Not all receivers_oi exist in the seurat object"))
  }

  # Check that scenario is either case_control or one_condition
  if (!scenario %in% c("case_control", "one_condition")){
    stop(paste0("scenario must be either 'case_control' or 'one_condition'"))
  }

  # Check that assay_oi is in seurat_obj
  if (!is.null(assay_oi)){
    if (!assay_oi %in% names(seuratObj@assays)) {
      stop(paste0("assay_oi '", assay_oi, "' does not exist in Seurat object"))
    }
    # Otherwise, use DefaultAssay
  } else {
    assay_oi <- DefaultAssay(seuratObj)
  }

  # If condition_colname is given, give message
  if (!is.null(condition_colname)) {
    message("condition_* is given. Only cells from that condition will be considered in cell type specificity calculation.")
  }

  # Calculate DE of ligands and receptors
  DE_table <- calculate_de(seuratObj,
                           celltype_colname = celltype_colname,
                           condition_colname = condition_colname,
                           condition_oi = condition_oi,
                           features = unique(unlist(lr_network_filtered)),
                           ...)

  # Average expression information
  expression_info <- get_exprs_avg(seuratObj,
                                   celltype_colname = celltype_colname,
                                   condition_colname = condition_colname,
                                   condition_oi = condition_oi,
                                   features = unique(unlist(lr_network_filtered)),
                                   assay_oi = assay_oi,
                                   ...)

  if (scenario == "case_control"){
    # Calculate condition specificity
    Idents(seuratObj) <- seuratObj[[condition_colname, drop=TRUE]]

    # Default settings to return all genes with their p-val and LFC
    FindMarkers_args <- list(object = seuratObj,
                             ident.1 = condition_oi,
                             ident.2 = condition_reference,
                             group.by = condition_colname,
                             assay = assay_oi,
                             features = unique(unlist(lr_network_filtered)),
                             min.pct = 0,
                             logfc.threshold = 0)

    # Replace this with user arguments
    FindMarkers_args[names(list(...))] <-   list(...)

    condition_markers <- do.call(FindMarkers, FindMarkers_args) %>%
                            rownames_to_column("gene")

    processed_condition_markers <- process_table_to_ic(condition_markers,
                                                       table_type = "group_DE",
                                                       lr_network_filtered)
  } else {
    processed_condition_markers <- NULL
  }

  # Combine DE of senders and receivers -> used for prioritization
  processed_DE_table <- process_table_to_ic(DE_table,
                                            table_type = "celltype_DE",
                                            lr_network_filtered,
                                            senders_oi = senders_oi,
                                            receivers_oi = receivers_oi)

  processed_expr_table <- process_table_to_ic(expression_info,
                                              table_type = "expression",
                                              lr_network_filtered)

  return (list(sender_receiver_de = processed_DE_table,
               sender_receiver_info = processed_expr_table,
               lr_condition_de = processed_condition_markers))

}

#' @title  Perform a prioritization of cell-cell interactions (similar to MultiNicheNet).
#' @description User can choose the importance attached to each of the following prioritization criteria: differential expression of ligand and receptor, cell-type specificity of expression of ligand and receptor, NicheNet ligand activity
#' @param sender_receiver_info Output of \code{generate_info_tables} OR \code{get_exprs_avg} -> \code{process_table_to_ic}
#' @param sender_receiver_de Output of \code{generate_info_tables} OR \code{calculate_de} -> \code{process_table_to_ic}
#' @param ligand_activities Output of \code{predict_ligand_activities}
#' @param lr_condition_de Output of \code{generate_info_tables} OR \code{FindMarkers} -> \code{process_table_to_ic}
#' @param prioritizing_weights Named vector indicating the relative weights of each prioritization criterion (default: NULL). If NULL, the weights are determined by the chosen "scenario". If provided, the vector must contain the following names: "de_ligand", "de_receptor", "activity_scaled", "exprs_ligand", "exprs_receptor", "ligand_condition_specificity", "receptor_condition_specificity"
#' @param scenario "case_control" or "one_condition": if "case_control", all weights are set to 1. If "one_condition", the weights are set to 0 for condition specificity and 1 for the remaining criteria.
#'
#' @return Data frames of prioritized sender-ligand-receiver-receptor interactions.
#' The resulting dataframe contains columns from the input dataframes, but columns from \code{lr_condition_de} are suffixed with \code{_group} (some columns from \code{lr_condition_de} are also not present).
#' Additionally, the following columns are added:
#' \itemize{
#' \item \code{lfc_pval_*}: product of -log10(pval) and the LFC of the ligand/receptor
#' \item \code{p_val_adapted_*}: p-value adapted to the sign of the LFC to only consider interactions where the ligand/receptor is upregulated in the sender/receiver
#' \item \code{activity_zscore}: z-score of the ligand activity
#' \item \code{prioritization_score}: The prioritization score for each interaction, calculated as a weighted sum of the prioritization criteria.
#' }
#' Moreover, \code{scaled_*} columns are scaled using the corresponding column's ranking or the \code{scale_quantile_adapted} function. The columns used for prioritization are scaled_p_val_adapted_ligand, scaled_p_val_adapted_receptor, scaled_activity, scaled_avg_exprs_ligand, scaled_avg_exprs_receptor, scaled_p_val_adapted_ligand_group, scaled_p_val_adapted_receptor_group
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'
#' # Calculate tables with generate_info_tables
#' info_tables <- generate_info_tables(seurat_obj, "celltype", sender_celltypes, receiver, expressed_ligands, expressed_receptors, lr_network, "aggregate", "LCMV", "SS")
#'
#' # Generate prioritization tables
#' generate_prioritization_tables(info_tables$sender_receiver_info,
#'                                info_tables$sender_receiver_de,
#'                                ligand_activities,
#'                                info_tables$lr_condition_de,
#'                                scenario = "case_control")
#'
#' # Alternatively, these tables can also be calculated manually:
#' # Calculate LCMV-specific average expression
#' expression_info = get_exprs_avg(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#'
#' # Calculate LCMV-specific cell-type markers
#' DE_table = calculate_de(seurat_obj, "celltype", condition_oi = "LCMV", condition_colname = "aggregate")
#'
#' # Calculate condition-specific markers
#' condition_markers <- FindMarkers(object = seuratObj, ident.1 = "LCMV", ident.2 = "SS",
#'                                  group.by = "aggregate", min.pct = 0, logfc.threshold = 0) %>% rownames_to_column("gene")
#'
#' # Process tables
#' processed_expr_info = process_table_to_ic(expression_info, table_type = "expression", lr_network)
#' processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network,
#'                                           senders_oi = sender_celltypes, receivers_oi = receiver)
#' processed_condition_DE_table <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network)
#'
#' # Custom weights
#' prioritizing_weights = c("de_ligand" = 1, "de_receptor" = 1, "activity_scaled" = 2, "exprs_ligand" = 1, "exprs_receptor" = 1, "ligand_condition_specificity" = 0.5, "receptor_condition_specificity" = 0.5)
#'
#' generate_prioritization_tables(processed_expr_info,processed_DE_table, ligand_activities, processed_condition_DE_table, prioritizing_weights)
#'
#'}
#' @export
generate_prioritization_tables = function(sender_receiver_info, sender_receiver_de, ligand_activities, lr_condition_de = NULL,
                                          prioritizing_weights = NULL,
                                          scenario = "case_control"){
  requireNamespace("dplyr")

  # If prioritizing weights is NULL, check that scenario is either case_control or one_condition
  if (is.null(prioritizing_weights)){
    if (!scenario %in% c("case_control", "one_condition")){
      stop(paste0("scenario must be either 'case_control' or 'one_condition'"))
    }
    if(scenario == "case_control"){
      # Check if lr_condition_de is NULL
      if (is.null(lr_condition_de)){
        stop("lr_condition_de is NULL. Please  provide lr_condition_de or change scenario to 'one_condition'")
      }

      prioritizing_weights = c("de_ligand" = 1,
                               "de_receptor" = 1,
                               "activity_scaled" = 1,
                               "exprs_ligand" = 1,
                               "exprs_receptor" = 1,
                               "ligand_condition_specificity" = 1,
                               "receptor_condition_specificity" = 1)
    }
    else if(scenario == "one_condition"){
      if (!is.null(lr_condition_de)){
        warning("lr_condition_de is provided but will not be used for prioritization. Change scenario to 'case_control' if condition specificity should also be considered.")
      }

      prioritizing_weights = c("de_ligand" = 1,
                               "de_receptor" = 1,
                               "activity_scaled" = 1,
                               "exprs_ligand" = 1,
                               "exprs_receptor" = 1,
                               "ligand_condition_specificity" = 0,
                               "receptor_condition_specificity" = 0)
    }

  } else {
    # Check that prioritizing weights has correct names
    weight_names <-  c("de_ligand", "de_receptor", "activity_scaled", "exprs_ligand", "exprs_receptor", "ligand_condition_specificity", "receptor_condition_specificity")
    if (!all(weight_names %in% names(prioritizing_weights))){
      stop(paste0("prioritizing_weights must have the following names:", weight_names))
    }
    # Check that prioritizing_weights is a named vector
    if (!is.vector(prioritizing_weights) | !is.numeric(prioritizing_weights)){
        stop("prioritizing_weights must be a named numeric vector")
    }

    # message using custom weights
    message("Using custom weights for prioritization")
    scenario <- NULL
  }

  # If "rank" column doesn't exist for ligand_activities, calculate rank
  if (!"rank" %in% colnames(ligand_activities)){
    ligand_activities <- ligand_activities %>% dplyr::mutate(rank = rank(desc(aupr_corrected)))
  }

  sender_receiver_tbl = sender_receiver_de %>% dplyr::distinct(sender, receiver)

  # Ligand DE prioritization
  sender_ligand_prioritization = sender_receiver_de %>% dplyr::ungroup() %>%
    dplyr::select(sender, ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
                  dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                                p_val_adapted_ligand = -log10(p_val_ligand)*sign(lfc_ligand))
  sender_ligand_prioritization = sender_ligand_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                scaled_p_val_adapted_ligand = rank(p_val_adapted_ligand, ties.method = "average", na.last = FALSE)/max(rank(p_val_adapted_ligand, ties.method = "average", na.last = FALSE))) %>%
                                                                  dplyr::arrange(-lfc_pval_ligand)

  # Receptor DE prioritization
  receiver_receptor_prioritization = sender_receiver_de %>% dplyr::ungroup() %>% dplyr::select(receiver, receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
                                      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                                                    p_val_adapted_receptor = -log10(p_val_receptor)*sign(lfc_receptor) )
  receiver_receptor_prioritization = receiver_receptor_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_adapted_receptor = rank(p_val_adapted_receptor, ties.method = "average", na.last = FALSE)/max(rank(p_val_adapted_receptor, ties.method = "average", na.last = FALSE))) %>% dplyr::arrange(-lfc_pval_receptor)

  # Ligand activity prioritization
  ligand_activity_prioritization = ligand_activities %>%
                                      {if ('receiver' %in% colnames(.)) select(., test_ligand, aupr_corrected, rank, receiver) else select(., test_ligand, aupr_corrected, rank)} %>%
                                      rename(activity=aupr_corrected, ligand=test_ligand) %>%
                                      dplyr::mutate(activity_zscore = nichenetr::scaling_zscore(activity),
                                                    scaled_activity = scale_quantile_adapted(activity, outlier_cutoff = 0.01)) %>% dplyr::arrange(-activity_zscore)


  # Cell-type specificity of expression of ligand:  per ligand scale across cell types
  ligand_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(sender, ligand, avg_ligand) %>% dplyr::distinct() %>% dplyr::group_by(ligand) %>%
                                               dplyr::mutate(scaled_avg_exprs_ligand = scale_quantile_adapted(avg_ligand)) %>% dplyr::arrange(-scaled_avg_exprs_ligand)

  # Cell-type specificity of expression of receptor:  per receptor scale across cell types
  receptor_celltype_specificity_prioritization = sender_receiver_info %>% dplyr::select(receiver, receptor, avg_receptor) %>% dplyr::distinct() %>% dplyr::group_by(receptor) %>%
                                               dplyr::mutate(scaled_avg_exprs_receptor = scale_quantile_adapted(avg_receptor)) %>% dplyr::arrange(-scaled_avg_exprs_receptor)

  if (!is.null(lr_condition_de)){
    # Condition specificity of ligand (upregulation)
    ligand_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(ligand, lfc_ligand, p_val_ligand) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_ligand = -log10(p_val_ligand)*lfc_ligand,
                    p_val_adapted_ligand = -log10(p_val_ligand)*sign(lfc_ligand))
    ligand_condition_prioritization = ligand_condition_prioritization %>% dplyr::mutate(scaled_lfc_ligand = rank(lfc_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_ligand = rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_ligand), ties.method = "average", na.last = FALSE)),
                                                                                        scaled_lfc_pval_ligand = rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_ligand, ties.method = "average", na.last = FALSE)),
                                                                                        scaled_p_val_adapted_ligand = rank(p_val_adapted_ligand, ties.method = "average", na.last = FALSE)/max(rank(p_val_adapted_ligand, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_ligand) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -ligand)

    # Condition specificity of receptor (upregulation)
    receptor_condition_prioritization = lr_condition_de %>% dplyr::ungroup() %>% dplyr::select(receptor, lfc_receptor, p_val_receptor) %>% dplyr::distinct() %>%
      dplyr::mutate(lfc_pval_receptor = -log10(p_val_receptor)*lfc_receptor,
                    p_val_adapted_receptor = -log10(p_val_receptor)*sign(lfc_receptor))
    receptor_condition_prioritization = receptor_condition_prioritization %>% dplyr::mutate(scaled_lfc_receptor = rank(lfc_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_receptor = rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)/max(rank(desc(p_val_receptor), ties.method = "average", na.last = FALSE)),
                                                                                            scaled_lfc_pval_receptor = rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)/max(rank(lfc_pval_receptor, ties.method = "average", na.last = FALSE)),
                                                                                            scaled_p_val_adapted_receptor = rank(p_val_adapted_receptor, ties.method = "average", na.last = FALSE)/max(rank(p_val_adapted_receptor, ties.method = "average", na.last = FALSE))) %>%
      dplyr::arrange(-lfc_pval_receptor) %>% rename_with(.fn = function(column_name) paste0(column_name, "_group"), .cols = -receptor)

  } else {
    if (any(prioritizing_weights[grep("specificity", names(prioritizing_weights))] > 0)) {
      stop("No lr_condition_de table given, yet the relevant weights are nonzero.\nEither set weights of 'ligand_condition_specificity' and 'receptor_condition_specificity' to zero or provide lr_condition_de.")
    }
  }

  weights <- prioritizing_weights
  # final group-based prioritization
  group_prioritization_tbl = sender_receiver_de %>%
    dplyr::inner_join(sender_receiver_info) %>%
    {if (weights["de_ligand"] > 0) dplyr::inner_join(., sender_ligand_prioritization) else (.)} %>%
    {if (weights["activity_scaled"] > 0) dplyr::inner_join(., ligand_activity_prioritization) else (.)} %>%
    {if (weights["de_receptor"] > 0) dplyr::inner_join(., receiver_receptor_prioritization) else (.)} %>%
    {if (weights["exprs_ligand"] > 0) dplyr::inner_join(., ligand_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["exprs_receptor"] > 0) dplyr::inner_join(., receptor_celltype_specificity_prioritization) else (.)} %>%
    {if (weights["ligand_condition_specificity"] > 0) dplyr::inner_join(., ligand_condition_prioritization) else (.)} %>%
    {if (weights["receptor_condition_specificity"] > 0) dplyr::inner_join(., receptor_condition_prioritization) else (.)}


  # have a weighted average the final score (no product!!)
  sum_prioritization_weights = weights["de_ligand"] + weights["de_receptor"] + weights["activity_scaled"] + weights["exprs_ligand"] + weights["exprs_receptor"] + weights["ligand_condition_specificity"] + weights["receptor_condition_specificity"]
  group_prioritization_tbl = group_prioritization_tbl %>% rowwise() %>%
    dplyr::mutate(prioritization_score =
                    (
                        (prioritizing_weights["de_ligand"] * ifelse("scaled_p_val_adapted_ligand" %in% names(group_prioritization_tbl), scaled_p_val_adapted_ligand, 0)) +
                        (prioritizing_weights["de_receptor"] * ifelse("scaled_p_val_adapted_receptor" %in% names(group_prioritization_tbl), scaled_p_val_adapted_receptor, 0)) +
                        (prioritizing_weights["activity_scaled"] * ifelse("scaled_activity" %in% names(group_prioritization_tbl), scaled_activity, 0)) +
                        (prioritizing_weights["exprs_ligand"] * ifelse("scaled_avg_exprs_ligand" %in% names(group_prioritization_tbl), scaled_avg_exprs_ligand, 0)) +
                        (prioritizing_weights["exprs_receptor"] * ifelse("scaled_avg_exprs_receptor" %in% names(group_prioritization_tbl), scaled_avg_exprs_receptor, 0)) +
                        (prioritizing_weights["ligand_condition_specificity"] * ifelse("scaled_p_val_adapted_ligand_group" %in% names(group_prioritization_tbl), scaled_p_val_adapted_ligand_group, 0)) +
                        (prioritizing_weights["receptor_condition_specificity"] * ifelse("scaled_p_val_adapted_receptor_group" %in% names(group_prioritization_tbl), scaled_p_val_adapted_receptor_group, 0))
                    )* (1/sum_prioritization_weights)) %>% dplyr::arrange(-prioritization_score) %>%
    ungroup()

  return (group_prioritization_tbl)

}

get_top_n_lr_pairs = function(prioritization_tables, top_n, groups_oi = NULL, senders_oi = NULL, receivers_oi = NULL, rank_per_group = TRUE){
  prioritization_tbl_oi = prioritization_tables$group_prioritization_tbl %>% dplyr::filter(group == top_group & fraction_expressing_ligand_receptor > 0) %>% dplyr::distinct(group, sender, receiver, ligand, receptor, receiver, id, prioritization_score)
  if(!is.null(groups_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(group %in% groups_oi)
  }
  if(!is.null(senders_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(sender %in% senders_oi)
  }
  if(!is.null(receivers_oi)){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::filter(receiver %in% receivers_oi)
  }
  if(rank_per_group == TRUE){
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::group_by(group) %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% dplyr::filter(prioritization_rank <= top_n)
  } else {
    prioritization_tbl_oi = prioritization_tbl_oi %>% dplyr::mutate(prioritization_rank = rank(desc(prioritization_score))) %>% dplyr::filter(prioritization_rank <= top_n)
  }
  return(prioritization_tbl_oi)
}
