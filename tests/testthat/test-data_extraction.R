context("Extract large data functions")
test_that("Complete networks and expression validation settings can be loaded", {
  expression_settings_validation = get_expression_settings_validation()
  expect_type(expression_settings_validation, "list")

  all_networks = get_all_networks()
  expect_type(all_networks, "list")

})
