# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Reload all code:              'Ctrl + Shift + L'!!!! #   Reload all code:              'Ctrl + Shift + L'!!!!
#   Reload Tall code:              'Ctrl + Shift + D'!!!! #   Document all code  devtools::document()



# install.packages("devtools")

devtools::use_package("tidyverse")
#' @import tidyverse




x = sample(1000)
y = sample(500)
devtools::use_data(x,mtcars,overwrite = T)
devtools::use_data(y,internal = TRUE,overwrite = T) # you can add all datasets in same dataset here!

