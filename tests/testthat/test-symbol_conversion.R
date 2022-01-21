context("Symbol conversion functions")

test_that("human-mouse and symbol-alias conversion works", {
  human_symbols = c("TNF","IFNG","IL-15","CTGF")
  aliases = human_symbols %>% convert_alias_to_symbols(organism = "human")
  expect_type(aliases,"character")

  mouse_symbols_v1 = human_symbols %>% convert_human_to_mouse_symbols() %>% .[!is.na(.)]
  expect_type(mouse_symbols_v1,"character")

  mouse_symbols_v2 = human_symbols %>% convert_human_to_mouse_symbols(version = 2) %>% .[!is.na(.)]
  expect_type(mouse_symbols_v2,"character")

  aliases_mouse = mouse_symbols_v1 %>% convert_alias_to_symbols(organism = "mouse")
  expect_type(aliases_mouse,"character")

  aliases_mouse_v2 = mouse_symbols_v2 %>% convert_alias_to_symbols(organism = "mouse")
  expect_type(aliases_mouse_v2,"character")

  human_symbols_back_v1 = mouse_symbols_v1 %>% convert_mouse_to_human_symbols(version = 1)
  expect_type(human_symbols_back_v1,"character")

  human_symbols_back_v2 = mouse_symbols_v2 %>% convert_mouse_to_human_symbols(version = 2)
  expect_type(human_symbols_back_v2,"character")

  aliases = human_symbols_back_v2 %>% convert_alias_to_symbols(organism = "human")
  expect_type(aliases,"character")
})
test_that("Seurat alias conversion works", {
  seurat_object_lite = readRDS(url("https://zenodo.org/record/3531889/files/seuratObj_test.rds"))
  seurat_object_lite2 = seurat_object_lite %>% alias_to_symbol_seurat(organism = "mouse")
  testthat::expect_equal(typeof(seurat_object_lite2), "S4")
})
