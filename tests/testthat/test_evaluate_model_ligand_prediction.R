context("Model evaluation functions: ligand prediction")
test_that("Convert expression settings to settings", {
  settings = lapply(expression_settings_validation,convert_expression_settings_evaluation)
  ligands = extract_ligands_from_settings(expression_settings_validation, combination = FALSE)
  ligands = unlist(ligands)

  settings_ligand_pred = convert_settings_ligand_prediction(settings = settings, all_ligands = ligands, validation = TRUE, single = TRUE)
  expect_type(settings_ligand_pred,"list")

  expect_type(settings_ligand_pred[[1]]$ligand,"character")
  expect_equal(length(settings_ligand_pred[[1]]$ligand),1)


  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")


  expect_equal(length(settings_ligand_pred),length(ligands)*length(settings))


  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = TRUE, single = FALSE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(settings))
  expect_type(settings_ligand_pred[[1]]$ligand,"character")
  expect_equal(length(settings_ligand_pred[[1]]$ligand),1)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = FALSE, single = TRUE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(ligands)*length(settings))
  expect_equal(settings_ligand_pred[[1]]$ligand,NULL)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

  settings_ligand_pred = convert_settings_ligand_prediction(settings, ligands, validation = FALSE, single = FALSE)
  expect_type(settings_ligand_pred,"list")
  expect_equal(length(settings_ligand_pred),length(settings))
  expect_equal(settings_ligand_pred[[1]]$ligand,NULL)
  expect_type(settings_ligand_pred[[1]]$name,"character")
  expect_type(settings_ligand_pred[[1]]$from,"character")
  expect_type(settings_ligand_pred[[1]]$response,"logical")

})

