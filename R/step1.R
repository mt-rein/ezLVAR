#### Step 1 core function ####

#' Step 1: Estimate the Measurement Model
#'
#' @description This function performs step 1 of 3S-LVAR: the estimation of the
#'   measurement model. It is a wrapper function around [lavaan::cfa]. Its
#'   output will be passed on to [step2].
#'
#'
#' @param data A data frame.
#' @param measurementmodel A string describing the measurement model using the
#'   `lavaan` syntax. Can be a list of strings determining the measurement
#'   blocks (e.g., one string for the MM of factor 1, and a second string for
#'   the MM of factor 2).
#' @param id_var String containing the name of the ID variable.
#' @param group_var String containing the name of the grouping variable.
#'   Optional.
#' @param group.equal Character vector containing the parameters constrained to
#'   equality across groups (cf. the `group.equal` argument in `lavaan`).
#'   Optional. Must be a list of strings if `measurementmodel` is a list (e.g.,
#'   one string for block 1, and a second string for block 2).
#' @param group.partial Character vector containing the non-invariances of the
#'   model (cf. the `group.partial` argument in `lavaan`). Optional. Must be a
#'   list of strings if `measurementmodel` is a list (e.g., one string for block
#'   1, and a second string for block 2).
#' @param ... Arguments to be passed on to [lavaan::cfa].
#'
#' @returns An object of class `3slvar_step1`, which is a list comprising the
#'   following elements: \item{mm_output}{the output of the measurement model
#'   estimation.} \item{data}{the original data set as supplied to the
#'   function.} \item{measurementmodel}{the measurement model as supplied to the
#'   function.} \item{id_var}{the name of the ID variable as supplied to the
#'   function.}
#'
#' @export
step1 <- function(
  data,
  measurementmodel,
  id_var,
  group_var = NULL,
  group.equal = NULL,
  group.partial = NULL,
  ...
) {
  #### Errors and Warnings ####
  # ensure that both group.partial and measurementmodel are a list:
  if (
    !is.null(group.partial) &&
      is.list(measurementmodel) &&
      !is.list(group.partial)
  ) {
    stop(
      "If the measurement model uses blocks, then non-invariances also need to be defined in a list."
    )
  }

  # ensure that both group.partial and measurementmodel are of the same
  # length
  if (
    is.list(measurementmodel) &&
      is.list(group.partial) &&
      length(measurementmodel) != length(group.partial)
  ) {
    stop(
      "The number of elements in the lists for measurement model and group.partial must be equal."
    )
  }

  #### Confirmatory Factor Analysis per block #####

  # if there are no MM blocks defined, turn the MM into a list with a single MM
  # block (to facilitate coding)
  if (!is.list(measurementmodel)) {
    measurementmodel <- list(measurementmodel)
    group.partial <- list(group.partial)
  }

  # how many blocks?
  n_blocks <- length(measurementmodel)

  mm_output <- vector("list", length = n_blocks)
  # estimate MM in each block:
  for (m in 1:n_blocks) {
    mm_output[[m]] <- lavaan::cfa(
      measurementmodel[[m]],
      data = data,
      group = group_var,
      meanstructure = TRUE,
      group.equal = group.equal[[m]],
      group.partial = group.partial[[m]],
      orthogonal = TRUE,
      ...
    )
  }

  # assemble output
  output <- list(
    "mm_output" = mm_output,
    "data" = data,
    "measurementmodel" = measurementmodel,
    "id_var" = id_var
  )
  class(output) <- "3slvar_step1"
  return(output)
}

#### summary function ####

#' @param step1output An object obtained with the [step1] function.
#' @param print_MM_summary Should the summary of the estimated measurement
#'   models be printed (i.e, call `lavaan::summary` on all blocks?
#' @param ... Arguments that are forwarded to `lavaan::summary()`.
#'
#' @rdname step1
#'
#' @export
summary.3slvar_step1 <- function(step1output, print_MM_summary = TRUE, ...) {
  fit_step1 <- step1output$mm_output
  n_blocks <- length(fit_step1)
  group_var <- lavaan::lavInspect(fit_step1[[1]], "group") |> as.character()
  group_names <- lavaan::lavInspect(fit_step1[[1]], "group.label")
  if (rlang::is_empty(group_var)) {
    n_groups <- 1
  } else {
    n_groups <- length(group_names)
  }

  cat("Number of measurement blocks: ", n_blocks, "\n")
  if (!rlang::is_empty(group_var)) {
    cat(paste0(
      "Grouping variable (number of groups): ",
      group_var,
      " (",
      n_groups,
      ")\n"
    ))
  }

  if (print_MM_summary) {
    for (block in 1:n_blocks) {
      cat("============================================================\n")
      cat("Summary of measurement block ", block, ":\n")
      print(lavaan::summary(fit_step1[[block]], ...))
      cat("\n\n")
    }
  }

  invisible(step1output)
}
