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
#   Reload all code:              'Ctrl + Shift + D'!!!! #   Document all code  devtools::document()

# devtools::use_package("tidyverse")

#' Say hello to whoever you want!
#'
#' return hello to everyone, everything you want
#'
#' blababjlballdfdfjdlfjdfkdfjdmfkldjfkdqjfmsldqjfmkldqfjsdqmlfjsdqmlfjsdqmfjqmfjsdqmfjqmj
#' dfdfdfdfddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddddfdfdf
#'
#' @param a Character vector. The name of whoever you want to greet. must follow following format
#' \describe{
#'  \item{grn}{Gene regulatory network}
#'  \item{sig}{signaling network}
#'  }
#'
#' @return a string
#' @examples
#' hello("Hadley")
#' hello("World")
#' @family test functions
hello <- function(a = "world") {
  # print(paste0(c("Hello,",a)))
  paste0(c("Hello,",a))
}
# hello("Hadley")



