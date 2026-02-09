#' Extract the data with factor scores
#'
#' @param step2output An object obtained with the [step2] function.
#'
#' @returns data The data as they were supplied to the [step1] function, with
#'   appended factor scores.
#'
#' @export
get_data <- function(step2output) {
  data <- step2output$data
  return(data)
}


#' Print an estimate data frame with rounded numeric values (internal only used
#' in summary functions)
#'
#' @param df A data frame
#' @param digits The number of digits to round numeric values to
#'
#' @returns df The data frame with numeric values rounded to the specified
#'   number of digits, printed without row names.
print_rounded <- function(df, digits) {
  num_cols <- sapply(df, is.numeric)
  df[num_cols] <- lapply(df[num_cols], round, digits = digits)
  print(df, row.names = FALSE)
}
