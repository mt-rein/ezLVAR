#' Step 1: Estimate the Measurement Model
#'
#' @description
#' This function performs step 1 of 3S-LVAR: the estimation of the measurement model. It is a wrapper function around `lavaan::cfa()`.
#' Its output will be passed on to [step2()].
#'
#'
#' @param data A data frame.
#' @param measurementmodel A string describing the measurement model using the `lavaan` syntax.
#'                         Can be a list of strings determining the measurement blocks (e.g., one string for the MM of factor 1, and a second string for the MM of factor 2).
#' @param id_var String containing the name of the ID variable.
#' @param group_var String containing the name of the grouping variable. Optional.
#' @param invariances Character vector containing the parameters constrained to equality across groups (cf. the `group.equal` argument in `lavaan`). Optional.
#' @param partial_noninvariances Character vector containing the non-invariances of the model (cf. the `group.partial` argument in `lavaan`). Optional. Must be a list of strings if `measurementmodel` is a list (e.g., one string for block 1, and a second string for block 2).
#' @param ... Arguments to be passed on to [lavaan::cfa()].
#'
#' @returns An object of class `3slvar_step1`, which is a list comprising the output of the measurement model estimation (`mm_output`), the original data set (`data`), the measurement model (`measurementmodel`), and the ID variable (`id_var`).
#' @export

step1 <- function(data, measurementmodel, id_var,
                  group_var = NULL,
                  invariances = NULL,
                  partial_noninvariances = NULL, ...) {

  #### Errors and Warnings ####
  # ensure that both partial_noninvariances and measurementmodel are a list:
  if (!is.null(partial_noninvariances) &&
      is.list(measurementmodel) &&
      !is.list(partial_noninvariances)) {
    stop("If the measurement model uses blocks, then non-invariances also need to be defined in a list.")
  }

  # ensure that both partial_noninvariances and measurementmodel are of the same length
  if (is.list(measurementmodel) &&
      is.list(partial_noninvariances) &&
      length(measurementmodel) != length(partial_noninvariances)) {
    stop("The number of elements in the lists for measurement model and partial_noninvariances must be equal.")
  }

  #### Confirmatory Factor Analysis per block #####
  # if there are no MM blocks defined, turn the MM into a list with a single MM block (to facilitate coding)
  if (!is.list(measurementmodel)) {
    measurementmodel <- list(measurementmodel)
    partial_noninvariances <- list(partial_noninvariances)
  }

  # how many blocks?
  n_blocks <- length(measurementmodel)

  mm_output <- vector("list", length = n_blocks)
  # estimate MM in each block:
  for (m in 1:n_blocks) {
    mm_output[[m]] <- lavaan::cfa(measurementmodel[[m]],
                                  data = data,
                                  group = group_var,
                                  meanstructure = TRUE,
                                  group.equal = invariances,
                                  group.partial = partial_noninvariances[[m]],
                                  orthogonal = TRUE,
                                  ...)
  }

  # assemble output
  output <- list("mm_output" = mm_output,
                 "data" = data,
                 "measurementmodel" = measurementmodel,
                 "id_var" = id_var)
  class(output) <- "3slvar_step1"
  return(output)
}
