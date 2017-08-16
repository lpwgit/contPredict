#' Check for scalar argument
#' This function checks for user's input is a scalar
#'
#' @details user's entry must meet the following criteria:
#'     character,
#'     length 1,
#'     not NA,
#'     non-zero length
#'
#' @param x character (1)
#'
#' @return boolean TRUE or FALSE
#'
#' @export
is_scalar_character <- function(x){

  is.character(x) && length (x) == 1 && !is.na(x) && nzchar(x)
}
