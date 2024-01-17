context("Prioritization scheme")
test_that("Prioritization scheme works", {
    options(timeout = 3600)
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
    seurat_obj_test = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
    seurat_obj_test = Seurat::UpdateSeuratObject(seurat_obj_test)

    lr_network = lr_network %>% distinct(from, to)
    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
    seurat_obj_test = alias_to_symbol_seurat(seurat_obj_test, "mouse")

    celltypes <- unique(seurat_obj_test$celltype)
    lr_network_renamed <- lr_network %>% rename(ligand=from, receptor=to)
    condition_oi <- "LCMV"
    condition_reference <- "SS"
    sender_celltypes <- "Mono"
    receiver <- "CD8 T"

    nichenet_output = suppressWarnings(nichenet_seuratobj_aggregate(seurat_obj = seurat_obj_test, receiver = receiver,
                condition_oi = condition_oi, condition_reference = condition_reference,  condition_colname = "aggregate",
                sender = sender_celltypes,
                ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network))

    expect_type(nichenet_output,"list")
    ligand_activities <- nichenet_output$ligand_activities

    feature_list <- intersect(rownames(seurat_obj_test), unique(union(lr_network_renamed$ligand, lr_network_renamed$receptor)))
    # Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
    DE_table <- calculate_de(seurat_obj_test, celltype_colname = "celltype",
                            condition_colname = "aggregate", condition_oi = condition_oi,
                            features = feature_list)

    # Average expression information - only for LCMV condition
    expression_info <- get_exprs_avg(seurat_obj_test, "celltype", condition_colname = "aggregate", condition_oi = condition_oi)

    # Calculate condition specificity - only for datasets with two conditions!
    condition_markers <- FindMarkers(object = seurat_obj_test, ident.1 = condition_oi, ident.2 = condition_reference,
                                    group.by = "aggregate", min.pct = 0, logfc.threshold = 0,
                                    features = feature_list) %>% rownames_to_column("gene")

    # Combine DE of senders and receivers -> used for prioritization
    processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network_renamed,
                                            senders_oi = sender_celltypes, receivers_oi = receiver)

    processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network_renamed)

    processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network_renamed)

    # Check processed tables
    expect_type(processed_DE_table,"list")
    expect_type(processed_expr_table,"list")
    expect_type(processed_condition_markers,"list")

    # Check that processed_DE_table has been processed correctly
    expect_equal(length(unique(processed_DE_table$sender)), length(sender_celltypes))
    expect_equal(length(unique(processed_DE_table$receiver)), length(receiver))

    expect_equal(processed_DE_table %>% filter(ligand == "Il1rn", sender == "Mono") %>% select(p_val_ligand, lfc_ligand, pct_expressed_sender, p_adj_ligand, sender, ligand) %>% distinct(across(everything())),
                 DE_table %>% filter(gene == "Il1rn", cluster_id == "Mono") %>% select(-pct.2),
                 check.attributes = FALSE)
    expect_equal(processed_DE_table %>% filter(ligand == "Il1rn", sender == "Mono", receptor == "Il1r2") %>% select(p_val_receptor, lfc_receptor, pct_expressed_receiver, p_adj_receptor, receiver, receptor) ,
                 DE_table %>% filter(gene == "Il1r2", cluster_id == "CD8 T") %>% select(-pct.2),
                 check.attributes = FALSE)
    temp_row <- processed_DE_table %>% filter(ligand == "Il1rn", sender == "Mono", receiver == "CD8 T", receptor == "Il1r2")
    expect_equal((temp_row$lfc_ligand + temp_row$lfc_receptor)/2, temp_row$ligand_receptor_lfc_avg)


    # Check that processed_expr_table has been processed correctly
    expect_equal(length(unique(processed_expr_table$sender)), length(celltypes))
    expect_equal(length(unique(processed_expr_table$receiver)), length(celltypes))

    temp_row <- processed_expr_table %>% filter(ligand == "Il1rn", sender == "Mono", receiver == "CD8 T", receptor == "Il1r2")
    expect_equal(temp_row$avg_ligand * temp_row$avg_receptor, temp_row$ligand_receptor_prod)

    # Check that processed_condition_markers has been processed correctly
    expect_equal(processed_condition_markers %>% filter(ligand == "Il1rn") %>% select(ligand, p_val_ligand, lfc_ligand, p_adj_ligand) %>% distinct(across(everything())),
                 condition_markers %>% filter(gene == "Il1rn") %>% select(-pct.1, -pct.2),
                 check.attributes = FALSE)
    expect_equal(processed_condition_markers %>% filter(receptor == "Il1r2") %>% select(receptor, p_val_receptor, lfc_receptor, p_adj_receptor) %>% distinct(across(everything())),
                condition_markers %>% filter(gene == "Il1r2") %>% select(-pct.1, -pct.2),
                check.attributes = FALSE)
    temp_row <- processed_condition_markers %>% filter(ligand == "Il1rn", receptor == "Il1r2")
    expect_equal((temp_row$lfc_ligand + temp_row$lfc_receptor)/2, temp_row$ligand_receptor_lfc_avg)

    # Check errors and warnings in case of improper usage
    expect_error(process_table_to_ic(condition_markers, table_type = "group_DE", lr_network = lr_network_renamed, receivers_oi = receiver))
    expect_error(process_table_to_ic(condition_markers, table_type = "group_DE", lr_network = lr_network_renamed, senders_oi = sender_celltypes))
    expect_error(process_table_to_ic(DE_table, table_type = "random", lr_network = lr_network_renamed))
    expect_warning(process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network = lr_network_renamed))
    expect_warning(process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network = lr_network_renamed))
    expect_warning(process_table_to_ic(expression_info, lr_network = lr_network_renamed, receivers_oi = receiver))
    expect_warning(process_table_to_ic(expression_info, lr_network = lr_network_renamed, senders_oi = sender_celltypes))

    # Default weights
    prioritizing_weights = c("de_ligand" = 1,
                            "de_receptor" = 1,
                            "activity_scaled" = 2,
                            "exprs_ligand" = 1,
                            "exprs_receptor" = 1,
                            "ligand_condition_specificity" = 0.5,
                            "receptor_condition_specificity" = 0.5)

    prior_table <- generate_prioritization_tables(processed_expr_table,
                               processed_DE_table,
                               ligand_activities,
                               processed_condition_markers,
                               prioritizing_weights)

    expect_type(prior_table,"list")

    # Check that columns contain columns from processed_DE_table, processed_expr_table, ligand_activities, and processed_condition_markers
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% select(colnames(processed_DE_table)),
                 processed_DE_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% mutate(sender = as.character(sender), receiver = as.character(receiver)),
                 check.attributes = FALSE)
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% select(colnames(processed_expr_table)),
                    processed_expr_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69", receiver == "CD8 T"),
                    check.attributes = FALSE)
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% pull(activity),
                    ligand_activities %>% filter(test_ligand == "Lgals1") %>% pull(aupr_corrected))
    temp_cols <- c("lfc_ligand", "lfc_receptor", "p_val_ligand", "p_val_receptor")
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% select(all_of(paste0(temp_cols, "_group"))),
                    processed_condition_markers %>% filter(ligand == "Lgals1", receptor == "Cd69") %>% select(all_of(temp_cols)),
                    check.attributes = FALSE)

    # Check that prioritization score is the same as the sum of the weighted scores
    temp_weights <- c(prioritizing_weights, prioritizing_weights["de_ligand"], prioritizing_weights["de_receptor"])
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>%
        mutate(prioritization_score = rowSums(across(c(scaled_lfc_ligand, scaled_lfc_receptor, scaled_activity, scaled_avg_exprs_ligand, scaled_avg_exprs_receptor, scaled_lfc_ligand_group, scaled_lfc_receptor_group, scaled_p_val_ligand_adapted, scaled_p_val_receptor_adapted)) * temp_weights) / sum(temp_weights)) %>% pull(prioritization_score),
                prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% pull(prioritization_score),
        check.attributes=FALSE)

    # Do not pass condition markers
    expect_error(generate_prioritization_tables(processed_expr_table,
                            processed_DE_table,
                            ligand_activities,
                            lr_condition_de = NULL,
                            prioritizing_weights))

    # Define priorization weights to 0 except for activity scaled
    prioritizing_weights = c("de_ligand" = 0,
                            "de_receptor" = 0,
                            "activity_scaled" = 1,
                            "exprs_ligand" = 0,
                            "exprs_receptor" = 0,
                            "ligand_condition_specificity" = 0,
                            "receptor_condition_specificity" = 0)

    prior_table <- generate_prioritization_tables(processed_expr_table,
                               processed_DE_table,
                               ligand_activities,
                               processed_condition_markers,
                               prioritizing_weights)

    # Check colnames
    expect_equal(colnames(prior_table),
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table), "activity", "rank", "activity_zscore", "scaled_activity", "prioritization_score")))


    prior_table_top10 <- prior_table %>% distinct(ligand, prioritization_score) %>% mutate(rank = rank(desc(prioritization_score), ties.method = "average")) %>% arrange(rank, ligand) %>% pull(ligand) %>% .[1:10]
    ligands_top10 <- ligand_activities %>% arrange(rank, test_ligand) %>% pull(test_ligand) %>% .[1:10]

    # When using only activity scaled, the top 10 ligands should be the same as the top 10 ligands from the ligand activity table
    expect_equal(prior_table_top10, ligands_top10)

    # All weights zero
    prioritizing_weights = c("de_ligand" = 0,
                             "de_receptor" = 0,
                             "activity_scaled" = 0,
                             "exprs_ligand" = 0,
                             "exprs_receptor" = 0,
                             "ligand_condition_specificity" = 0,
                             "receptor_condition_specificity" = 0)

    prior_table <- generate_prioritization_tables(processed_expr_table,
                                                  processed_DE_table,
                                                  ligand_activities,
                                                  processed_condition_markers,
                                                  prioritizing_weights)

    # Check colnames
    expect_equal(colnames(prior_table),
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table), "prioritization_score")))

    prioritizing_weights["de_ligand"] = 1
    prior_table <- generate_prioritization_tables(processed_expr_table,
                                                  processed_DE_table,
                                                  ligand_activities,
                                                  processed_condition_markers,
                                                  prioritizing_weights)

    expect_equal(colnames(prior_table),
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table),
                          "lfc_pval_ligand", "p_val_ligand_adapted", "scaled_lfc_ligand", "scaled_p_val_ligand", "scaled_lfc_pval_ligand", "scaled_p_val_ligand_adapted", "prioritization_score")))



})
