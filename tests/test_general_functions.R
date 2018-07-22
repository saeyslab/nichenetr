library(dplyr)
library(testthat)

test_check("nichenetr",filter = "use_ligand_to_target")
test_check("nichenetr",filter = "calculate_popularity_bias")
test_check("nichenetr",filter = "random_networks")
# test_check("nichenetr",filter = "data_extraction") # in previous versions, networks and expression settings were delivered together with the package - no not anymore to reduce installation time, so no test required to load this data.
