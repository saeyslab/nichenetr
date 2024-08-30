context("Prioritization scheme")
test_that("Wrapper function for seurat", {
  options(timeout = 3600)

  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  lr_network = lr_network %>% distinct(from, to)

  seurat_obj_test = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
  seurat_obj_test = Seurat::UpdateSeuratObject(seurat_obj_test)

  if (as.numeric(substr(packageVersion("Seurat"), 1, 1)) < 5){
    seurat_obj_test = alias_to_symbol_seurat(seurat_obj_test, "mouse")
  } else if (grepl("^5", packageVersion("Seurat")) & grepl("^5", seurat_obj_test@version)){
    expect_warning(alias_to_symbol_seurat(seurat_obj_test, "mouse"))
  }

  lr_network_filtered <- lr_network %>% filter(from %in% rownames(seurat_obj_test), to %in% rownames(seurat_obj_test))

  generate_info_tables_args <- list(
    seuratObj = seurat_obj_test,
    celltype_colname = "celltype",
    senders_oi = "Mono",
    receivers_oi = "CD8 T",
    lr_network = lr_network_filtered,
    condition_colname = "aggregate",
    condition_oi = "LCMV",
    condition_reference = "SS",
    scenario = "case_control",
    assay_oi = "RNA"
  )

  tmp <- do.call(generate_info_tables, generate_info_tables_args)

  expect_type(tmp, "list")
  expect_equal(length(tmp), 3)
  expect_named(tmp, c("sender_receiver_de", "sender_receiver_info", "lr_condition_de"))

  # Cell type, senders, receivers, condition doesn't exist
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "celltype_colname", "cell.type")), "does not exist")
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "senders_oi", "Mono1")), "Not all senders_oi exist")
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "receivers_oi", "CD8 T1")), "Not all receivers_oi exist")
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "condition_colname", "condition")), "does not exist")

  # Not all conditions are given
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "condition_colname", NULL)), "arguments must be either all NULL or all provided")
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args, "condition_oi", NULL)), "arguments must be either all NULL or all provided")
  expect_error(do.call(generate_info_tables, replace(generate_info_tables_args,
                                                     c("condition_colname", "condition_oi", "condition_reference"),
                                                      c(NULL, NULL, NULL))),
               "arguments are not provided but the 'case_control' scenario is selected")

  # Conditions given
  expect_warning(do.call(generate_info_tables, replace(generate_info_tables_args, "scenario", "one_condition")),
                 "arguments are provided but the 'one_condition' scenario is selected")

  tmp2 <- do.call(generate_info_tables, replace(generate_info_tables_args, "scenario", "one_condition"))
  expect_null(tmp2$lr_condition_de)

  # Give additional arguments - change DE test
  tmp2 <- do.call(generate_info_tables, replace(generate_info_tables_args, "test.use", "t"))

  # pval should be different but LFC should remain the same
  expect_equal(nrow(tmp$sender_receiver_de), nrow(tmp2$sender_receiver_de))
  expect_false(isTRUE(all.equal(tmp$sender_receiver_de$p_val_ligand, tmp2$sender_receiver_de$p_val_ligand)))
  expect_identical(tmp2$sender_receiver_de$lfc_ligand, tmp$sender_receiver_de$lfc_ligand)

  # Now for lr_condition_de
  expect_equal(nrow(tmp$lr_condition_de), nrow(tmp2$lr_condition_de))
  expect_false(isTRUE(all.equal(tmp$lr_condition_de$p_val_ligand, tmp2$lr_condition_de$p_val_ligand)))
  expect_identical(tmp2$lr_condition_de$lfc_ligand, tmp$lr_condition_de$lfc_ligand)

  # Averae expression should remain the same
  expect_identical(tmp2$sender_receiver_info, tmp$sender_receiver_info)

  # Change logfc.threshold and min.pct
  tmp2 <- do.call(generate_info_tables, replace(generate_info_tables_args, c("logfc.threshold", "min.pct"), c(0.25, 0.2)))
  expect_true(nrow(tmp2$sender_receiver_de) < nrow(tmp$sender_receiver_de))
  expect_true(nrow(tmp2$lr_condition_de) < nrow(tmp$lr_condition_de))
  expect_identical(tmp2$sender_receiver_info, tmp$sender_receiver_info)

  # Change slot if object is seurat obj v5
  # slot_name <- ifelse(as.numeric(substr(seurat_obj_test@version, 1, 1)) >= 5 & inherits(seurat_obj_test[["RNA"]], "Assay5"),
  #     "layer", "slot")
  slot_name <- "slot"

  tmp3 <- do.call(generate_info_tables, replace(generate_info_tables_args, slot_name, "counts"))

  # tmp3 should have different values for p_val_ligand, lfc_ligand, average values, and lr_condition_de
  expect_false(isTRUE(all.equal(tmp$sender_receiver_de$p_val_receptor, tmp3$sender_receiver_de$p_val_receptor)))
  # LFC remains the same for v5
  if (grepl("^5", packageVersion("Seurat")) & grepl("^5", seurat_obj_test@version)){
    expect_true(isTRUE(all.equal(tmp$sender_receiver_de$lfc_receptor, tmp3$sender_receiver_de$lfc_receptor)))
  } else {
    expect_false(isTRUE(all.equal(tmp$sender_receiver_de$lfc_receptor, tmp3$sender_receiver_de$lfc_receptor)))
  }
  expect_false(isTRUE(all.equal(tmp$sender_receiver_info$avg_ligand, tmp3$sender_receiver_info$avg_ligand)))
  expect_false(isTRUE(all.equal(tmp$lr_condition_de$p_val_receptor, tmp3$lr_condition_de$p_val_receptor)))

  # It should have the same length for everything though
  expect_equal(nrow(tmp$sender_receiver_info), nrow(tmp3$sender_receiver_info))
  expect_equal(nrow(tmp$sender_receiver_de), nrow(tmp3$sender_receiver_de))
  expect_equal(nrow(tmp$lr_condition_de), nrow(tmp3$lr_condition_de))

  # Results should be different for one_condition if we give NULL conditions vs not
  tmp2 <- do.call(generate_info_tables, replace(generate_info_tables_args, "scenario", "one_condition"))
  tmp3 <- do.call(generate_info_tables, generate_info_tables_args %>% replace(., c("condition_colname", "condition_oi", "condition_reference"), c(NULL, NULL, NULL)) %>% replace("scenario", "one_condition"))

  expect_false(isTRUE(all.equal(tmp3$sender_receiver_de, tmp2$sender_receiver_de)))
  expect_false(isTRUE(all.equal(tmp3$sender_receiver_info, tmp2$sender_receiver_info)))


})

test_that("Prioritization scheme works", {
    options(timeout = 3600)

    lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
    lr_network = lr_network %>% distinct(from, to)
    ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
    weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))

    seurat_obj_test = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))

    seurat_objs = list(Seurat::UpdateSeuratObject(seurat_obj_test))
    pcts <- c(0.1, 0.05)

    # Test for v5
    if (grepl("^5", packageVersion("Seurat"))){
      seurat_objs[[2]] <- CreateSeuratObject(counts = Seurat::GetAssayData(seurat_obj_test, layer = "counts"),
                                             meta.data = seurat_obj_test@meta.data) %>% SetIdent(value = .$celltype) %>%
        NormalizeData()
    }


    # Test both newly created object and updated object from v3
    for (i in 1:length(seurat_objs)){
      seurat_obj_test <- seurat_objs[[i]]
      pct <- pcts[i]

      if (grepl("^5", packageVersion("Seurat"))){
        expect_warning(alias_to_symbol_seurat(seurat_obj_test, "mouse"))
      } else {
        seurat_obj_test <- alias_to_symbol_seurat(seurat_obj_test, "mouse")
      }


    weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))

    celltypes <- unique(seurat_obj_test$celltype)
    condition_oi <- "LCMV"
    condition_reference <- "SS"
    sender_celltypes <- "Mono"
    receiver <- "CD8 T"

      nichenet_output = suppressWarnings(nichenet_seuratobj_aggregate(seurat_obj = seurat_obj_test, receiver = receiver,
                                                                      condition_oi = condition_oi, condition_reference = condition_reference,  condition_colname = "aggregate",
                                                                      sender = sender_celltypes, expression_pct = pct,
                                                                      ligand_target_matrix = ligand_target_matrix, weighted_networks = weighted_networks, lr_network = lr_network))

      expect_type(nichenet_output,"list")
      ligand_activities <- nichenet_output$ligand_activities

    feature_list <- unique(unlist(lr_network)) %>% .[. %in% rownames(seurat_obj_test)]
    # Only calculate DE for LCMV condition, with genes that are in the ligand-receptor network
    DE_table <- calculate_de(seurat_obj_test, celltype_colname = "celltype",
                            condition_colname = "aggregate", condition_oi = condition_oi,
                            features = feature_list)


    # Average expression information - only for LCMV condition
    expression_info <- get_exprs_avg(seurat_obj_test, "celltype",
                                     condition_colname = "aggregate", condition_oi = condition_oi,
                                     features = feature_list)

    # Test cell type names conversion for Seurat object
    # Replace space with underscore
    seurat_obj_test$celltype2 <- gsub(" ", "_", seurat_obj_test$celltype)
    new_celltypes <- suppressWarnings(get_exprs_avg(seurat_obj_test, "celltype2") %>% pull(cluster_id) %>% unique())
    expect_equal(new_celltypes, sort(unique(seurat_obj_test$celltype2)))

    # Replace CD8 T with CD8_T-test & replace Mono with Mono-test
    seurat_obj_test$celltype2 <- gsub("CD8 T", "CD8_T-test", seurat_obj_test$celltype) %>% gsub("Mono", "Mono-test", .)
    new_celltypes <- suppressWarnings(get_exprs_avg(seurat_obj_test, "celltype2") %>% pull(cluster_id) %>% unique())
    expect_equal(new_celltypes, sort(unique(seurat_obj_test$celltype2)))

    # Calculate condition specificity - only for datasets with two conditions!
    condition_markers <- FindMarkers(object = seurat_obj_test, ident.1 = condition_oi, ident.2 = condition_reference,
                                       group.by = "aggregate", min.pct = 0, logfc.threshold = 0,
                                       features = feature_list) %>% rownames_to_column("gene")

    # Combine DE of senders and receivers -> used for prioritization
    processed_DE_table <- process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network,
                                            senders_oi = sender_celltypes, receivers_oi = receiver)

    processed_expr_table <- process_table_to_ic(expression_info, table_type = "expression", lr_network)

    processed_condition_markers <- process_table_to_ic(condition_markers, table_type = "group_DE", lr_network)

    lr_network_filterd <- lr_network %>% filter(from %in% feature_list, to %in% feature_list)
    # Check that results from wrapper function generate_info_tables and the separate ones are the same
    info_tables <- generate_info_tables(seurat_obj_test, "celltype", sender_celltypes, receiver,
                                        lr_network_filterd, "aggregate",
                                        condition_oi, condition_reference, "case_control", "RNA")

    expect_identical(info_tables$sender_receiver_info, processed_expr_table)
    expect_identical(info_tables$sender_receiver_de, processed_DE_table)
    expect_identical(info_tables$lr_condition_de, processed_condition_markers)

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
    expect_error(process_table_to_ic(condition_markers, table_type = "group_DE", lr_network = lr_network, receivers_oi = receiver))
    expect_error(process_table_to_ic(condition_markers, table_type = "group_DE", lr_network = lr_network, senders_oi = sender_celltypes))
    expect_error(process_table_to_ic(DE_table, table_type = "random", lr_network = lr_network))
    expect_warning(process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network = lr_network))
    expect_warning(process_table_to_ic(DE_table, table_type = "celltype_DE", lr_network = lr_network))
    expect_warning(process_table_to_ic(expression_info, lr_network = lr_network, receivers_oi = receiver))
    expect_warning(process_table_to_ic(expression_info, lr_network = lr_network, senders_oi = sender_celltypes))

    # Default weights
    prior_table <- generate_prioritization_tables(processed_expr_table,
                                                  processed_DE_table,
                                                  ligand_activities,
                                                  processed_condition_markers)
    default_prior_table <- prior_table
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
    temp_weights = c("de_ligand" = 0.5,
                     "de_receptor" = 0.5,
                     "activity_scaled" = 1,
                     "exprs_ligand" = 0.5,
                     "exprs_receptor" = 0.5,
                     "ligand_condition_specificity" = 1,
                     "receptor_condition_specificity" = 1)
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>%
        mutate(prioritization_score = rowSums(across(c(scaled_p_val_adapted_ligand, scaled_p_val_adapted_receptor, scaled_activity, scaled_avg_exprs_ligand, scaled_avg_exprs_receptor, scaled_p_val_adapted_ligand_group, scaled_p_val_adapted_receptor_group)) * temp_weights) / sum(temp_weights)) %>% pull(prioritization_score),
                prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% pull(prioritization_score),
        check.attributes=FALSE)

    # Custom weights
    prioritizing_weights <- c("de_ligand" = 1,
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
                                                  prioritizing_weights = prioritizing_weights)
    temp_weights <- prioritizing_weights
    temp_weights[c("de_ligand", "de_receptor", "exprs_ligand", "exprs_receptor")] <- 0.5
    expect_equal(prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>%
                   mutate(prioritization_score = rowSums(across(c(scaled_p_val_adapted_ligand, scaled_p_val_adapted_receptor, scaled_activity, scaled_avg_exprs_ligand, scaled_avg_exprs_receptor, scaled_p_val_adapted_ligand_group, scaled_p_val_adapted_receptor_group)) * temp_weights) / sum(temp_weights)) %>% pull(prioritization_score),
                 prior_table %>% filter(ligand == "Lgals1", sender == "Mono", receptor == "Cd69") %>% pull(prioritization_score),
                 check.attributes=FALSE)

    # Check "one_condition" scenario
    prior_table <- generate_prioritization_tables(processed_expr_table,
                                                  processed_DE_table,
                                                  ligand_activities,
                                                  scenario = "one_condition")
    # Check that colnames don't have "_group
    expect_false(any(grepl("_group", colnames(prior_table))))

    # Try giving lr_condition_de
    expect_warning(generate_prioritization_tables(processed_expr_table,
                                  processed_DE_table,
                                  ligand_activities,
                                  lr_condition_de = processed_condition_markers,
                                  scenario = "one_condition"))

    # Do not pass condition markers
    expect_error(generate_prioritization_tables(processed_expr_table,
                                                processed_DE_table,
                                                ligand_activities,
                                                lr_condition_de = NULL,
                                                scenario = "case_control"))

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
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table),
                          "activity", "rank", "activity_zscore", "scaled_activity", "prioritization_score", "prioritization_rank")))


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
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table), "prioritization_score", "prioritization_rank")))

    prioritizing_weights["de_ligand"] = 1
    prior_table <- generate_prioritization_tables(processed_expr_table,
                                                  processed_DE_table,
                                                  ligand_activities,
                                                  processed_condition_markers,
                                                  prioritizing_weights)


    expect_equal(colnames(prior_table),
                 unique(c(colnames(processed_DE_table), colnames(processed_expr_table),
                          "lfc_pval_ligand", "p_val_adapted_ligand", "scaled_lfc_ligand", "scaled_p_val_ligand", "scaled_lfc_pval_ligand", "scaled_p_val_adapted_ligand", "prioritization_score", "prioritization_rank")))


    # Check mushroom plot
    usable_colnames <- default_prior_table %>% select(ends_with(c("_ligand", "_receptor", "_sender", "_receiver"))) %>% colnames() %>%
       stringr::str_remove("_ligand|_receptor|_sender|_receiver") %>% unique

    mushroom_plots <- list()
    for (colname in usable_colnames){
      # If range of values is between 0 and 0.001
      ext <- c("_ligand", "_receptor")
      if (colname == "pct_expressed") ext <- c("_sender", "_receiver")
      values_range <- range(default_prior_table[, c(paste0(colname, ext[1]), paste0(colname, ext[2]))])

      # Only scaled values can occupy both size and color
      if (values_range[1] >= 0 & values_range[2] <= 1.001){
        mushroom_plots[[colname]] <- make_mushroom_plot(default_prior_table, top_n = 30, size = colname, color = colname)
      } else {
        # If range of values is not between 0 and 0.001, it should throw an error when used as size
        expect_error(make_mushroom_plot(default_prior_table, top_n = 30, size = colname, color = colname))

        # But ok when used as color
        mushroom_plots[[colname]] <- make_mushroom_plot(default_prior_table, top_n = 30, size = "scaled_lfc", color = colname, true_color_range = TRUE)
      }
    }

    # Expect all to be ggplot objects
    expect_true(all(sapply(mushroom_plots, inherits, "gg")))

    # If a column name that doesn't exist is passed, it should throw an error
    expect_error(make_mushroom_plot(default_prior_table, top_n = 30, size = "non_existent_colname", color = usable_colnames[1]))
    expect_error(make_mushroom_plot(default_prior_table, top_n = 30, size = usable_colnames[1], color = "non_existent_colname"))
    }





})
