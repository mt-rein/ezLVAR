#' Extract the data with factor scores
#'
#' @param step2output An object obtained with the [step2] function.
#'
#' @returns data The data as they were supplied to the [step1] function, with appended factor scores.
#' @export
get_data <- function(step2output) {
  data <- step2output$data
  return(data)
}
