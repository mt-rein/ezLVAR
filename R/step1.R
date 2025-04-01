#' Step 1: Estimate the Measurement Model
#'
#' @description
#' This function performs step 1 of 3S-LVAR: the estimation of the measurement model. It is a wrapper function around `lavaan::cfa()`.
#' Its output will be passed on to `step2()`.
#'
#'
#' @param data A data frame.
#' @param measurementmodel A string describing the measurement model using the lavaan syntax.
#'                         Can be a list of strings determining the number of measurement blocks (e.g., one string for the MM of factor 1, and a second string for the MM of factor 2).
#' @param group String containing the name of the grouping variable. Optional.
#' @param invariances String containing the parameters constrained to equality across groups (cf. the `group.equal` argument in lavaan)
#' @param partial_noninvariances String containing the non-invariances of the model (cf. the `group.partial` argument in lavaan). Must be a list of strings if `measurementmodel` is a list (e.g., one string for block 1, and a second string for block 2).
#' @param ... Arguments to be passed on to `lavaan::cfa()`.
#'
#' @return A list comprising the output of the measurement model estimation (`MMoutput`), the original data set (`data`), and the measurement model (`measurement model`).
#' @export

step1 <- function(data, measurementmodel, group = NULL,
                  invariances = NULL,
                  partial_noninvariances = NULL, ...){

  if(is.list(measurementmodel) & !is.list(partial_noninvariances) & !is.null(partial_noninvariances)){
    stop("If the measurement model uses blocks, then non-invariances also need to be defined in a list.")
  }

  if(is.list(measurementmodel) & is.list(partial_noninvariances) & !is.null(partial_noninvariances)){
    if(length(measurementmodel) != length(partial_noninvariances)){
      stop("The number of elements in the lists for measurement model and partial_noninvariances must be equal.")
    }
  }

  data <- data |> as.data.frame()

  # Are there measurement blocks in the MM?
  if(is.list(measurementmodel)){
    # If there are measurement blocks, how many do we have?
    M <- length(measurementmodel)

    # create list, with elements equal to number of measurement blocks
    MMoutput <- vector("list", length = M)
    for(m in 1:M){
      # estimate MM in each block
      MMoutput[[m]] <- lavaan::cfa(measurementmodel[[m]],
                                   data = data,
                                   group = group,
                                   group.equal = invariances,
                                   group.partial = partial_noninvariances[[m]],
                                   ...)
    }
  } else {
    # if there are no measurement blocks, estimate full MM
    MMoutput <- lavaan::cfa(measurementmodel,
                            data = data,
                            group = group,
                            group.equal = invariances,
                            group.partial = partial_noninvariances,
                            ...)
  }

  # assemble output
  output <- list("MMoutput" = MMoutput,
                 "data" = data,
                 "measurementmodel" = measurementmodel)
  return(output)
}
