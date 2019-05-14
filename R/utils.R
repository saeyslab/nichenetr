#' @title Normalize values in a vector by the z-score method
#'
#' @description \code{scaling_zscore} Normalize values in a vector by the z-score method.
#'
#' @usage
#' scaling_zscore(x)
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' scaling_zscore(rnorm(5))
#'
#' @export
#'
scaling_zscore = function(x){
  if (typeof(x) == "double"){
    if(sd(x, na.rm = TRUE) > 0){
      return((x - mean(x))/sd(x))
    } else{
      return((x - mean(x)))
    }
  } else {return(x)}
}
#' @title Normalize values in a vector by the modified z-score method.
#'
#' @description \code{scaling_modified_zscore} Normalize values by the modified z-score method (uses median and median absolute deviation instead meean)
#'
#' @usage
#' scaling_modified_zscore(x)
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector.
#'
#' @examples
#' scaling_modified_zscore(rnorm(5))
#'
#' @export
#'
scaling_modified_zscore = function (x) {
  if (typeof(x) == "double") {
    if (mad(x, na.rm = TRUE) != 0) {
      return(0.6745 * (x - median(x))/mad(x))
    }
    else {
      return(0.6745 * (x - median(x)))
    }
  }
  else {
    return(x)
  }
}

mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
