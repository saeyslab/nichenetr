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
#' a = hello("Hadley")
#' print(a)
#' hello("World")
#' @family test functions #### does not work yet
#'
#' @export
#'
hello = function(a = "world") {
  # print(paste0(c("Hello,",a)))
  paste0("Hello, ",a)
}



