library(dplyr)
library(Seurat)
library(testthat)

test_check("nichenetr", filter = "application_visualization")
test_check("nichenetr", filter = "application_networks")
test_check("nichenetr", filter = "application_prediction")
