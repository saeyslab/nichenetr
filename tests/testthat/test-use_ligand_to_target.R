context("Model usage functions")

test_that("Get top n or percentage targets or ligands: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","BMP2",c("IL4","IL13"))
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)

  expect_type(extract_top_fraction_targets("BMP2",0.01,ligand_target_matrix),"double")
  expect_type(extract_top_n_targets("BMP2",100,ligand_target_matrix),"double")
  expect_error(extract_top_n_targets("BMP2",100000,ligand_target_matrix))
  expect_error(extract_top_n_targets("BMP4",100000,ligand_target_matrix))

  expect_type(extract_top_fraction_ligands("ID3",0.01,ligand_target_matrix),"double")
  expect_error(extract_top_n_ligands("ID3",100000,ligand_target_matrix))
  expect_type(extract_top_n_ligands("ID3",1,ligand_target_matrix),"double")
})
test_that("Get targets genes of a ligand: no error", {
  weighted_networks = construct_weighted_networks(lr_network, sig_network, gr_network,source_weights_df)
  ligands = list("TNF","BMP2",c("IL4","IL13"))
  ligand_target_matrix = construct_ligand_target_matrix(weighted_networks, ligands)

  expect_type(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "distribution", fdr_method = "global", output = "logical",ligands_position = "cols"),"logical")
  expect_type(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "distribution", fdr_method = "global", output = "logical",ligands_position = "rows"),"logical")

  expect_type(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "quantile", fdr_method = "local", output = "logical",ligands_position = "cols"),"logical")
  expect_type(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "quantile", fdr_method = "local", output = "gene_symbols",ligands_position = "cols"),"character")
  expect_equal(length(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "quantile", fdr_method = "local", output = "logical",ligands_position = "cols")),nrow(ligand_target_matrix))
  expect_equal(length(get_target_genes_ligand_oi("BMP2", ligand_target_matrix %>% t(), cutoff_method = "quantile", fdr_method = "local", output = "logical",ligands_position = "rows")), nrow(ligand_target_matrix))

  expect_error(get_target_genes_ligand_oi("BMP2", ligand_target_matrix, cutoff_method = "decoy", fdr_method = "local", output = "logical",ligands_position = "cols"))
})

