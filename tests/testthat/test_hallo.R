context("Hallo function")

test_that("string is truely outputted",{
  expect_equal(hello("Hadley"),"Hello, Hadley")
  expect_identical(hello("Hadley"),"Hello, Hadley")

})
