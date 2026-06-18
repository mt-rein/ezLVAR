#' Creates an Innovation Covariance Matrix
#'
#' @description This is a helper function to create the Q matrix, which is
#'   required to specify the structural model for [step3]. The Q matrix
#'   describes the innovation (co)variances (i.e., the residuals of the latent
#'   variables).
#'
#' @param step2output An object obtained with the [step2] function.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`,
#'   `labels`, `lbound`, and `ubound` are expanded to accommodate the
#'   specification of a random intercept.
#' @param startvalues A square matrix of numeric values that represents the
#'   starting values for each parameter. The number of rows/columns must be
#'   equal to the number of latent factors in the model. Optional. If `NULL`,
#'   the starting values for all innovation variances (on the diagonal) are set
#'   to the variances of the respective factor score variables, and the
#'   innovation covariances (on the off-diagonal) are set to 0.
#' @param free A square matrix of TRUE and FALSE values that indicates which
#'   parameters are freely estimated. The number of rows/columns must be equal
#'   to the number of latent factors in the model. Optional. If `NULL`, all
#'   elements of the innovation covariance matrix are freely estimated.
#' @param labels A square matrix of strings that indicates the labels for each
#'   parameter. The number of rows/columns must be equal to the number of latent
#'   factors in the model. Optional. If `NULL`, labels will be automatically
#'   generated.
#' @param lbound A square matrix of numeric values that indicates the lower
#'   bounds for each parameter (if a value is NA, no bounds are imposed on that
#'   parameter). The number of rows/columns must be equal to the number of
#'   latent factors in the model. Optional. If `NULL`, no bounds are imposed.
#' @param ubound A square matrix of numeric values that indicates the upper
#'   bounds for each parameter (if a value is NA, no bounds are imposed on that
#'   parameter). The number of rows/columns must be equal to the number of
#'   latent factors in the model. Optional. If `NULL`, no bounds are imposed.
#'
#' @returns Q An `OpenMx` matrix object that is entered into [step3()].
#'
#' @export
create_Q <- function(
  step2output,
  random_intercept = FALSE,
  startvalues = NULL,
  free = NULL,
  labels = NULL,
  lbound = NULL,
  ubound = NULL
) {
  factors <- step2output$other$factors
  n_factors <- length(factors)

  #### Errors ####
  matrices <- list(
    "startvalues" = startvalues,
    "free" = free,
    "labels" = labels,
    "lbound" = lbound,
    "ubound" = ubound
  )

  if (!is.logical(random_intercept)) {
    stop("The random_intercept argument must be TRUE or FALSE.")
  }

  # all matrices must have n_factors rows and columns
  for (mat in names(matrices)) {
    if (
      !is.null(matrices[[mat]]) &&
        (nrow(matrices[[mat]]) != n_factors ||
          ncol(matrices[[mat]]) != n_factors)
    ) {
      stop(glue::glue(
        "The '{mat}' element must have the same number of rows and columns as the number of latent factors in the model."
      ))
    }
  }

  #### create elements that have not been provided by the user ####
  # create matrix with start values (the (co)variances of the factor scores
  # variables)
  if (is.null(startvalues)) {
    startvalues <- var(step2output$data[, factors], na.rm = TRUE)
  }

  # create matrix that indicates free parameters
  if (is.null(free)) {
    free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
  }

  # create label matrix
  if (is.null(labels)) {
    labels <- matrix(NA, n_factors, n_factors)

    # Fill the matrix with symmetric labels
    for (i in 1:n_factors) {
      for (j in 1:n_factors) {
        labels[i, j] <- paste0(
          "zeta_",
          factors[min(i, j)],
          "_",
          factors[max(i, j)]
        )
      }
    }
  }

  # create lbound and ubound objects if not specified by user
  if (is.null(lbound)) {
    lbound <- matrix(NA, n_factors, n_factors)
  }
  if (is.null(ubound)) {
    ubound <- matrix(NA, n_factors, n_factors)
  }

  #### if random intercept has been requested, expand the matrices ####
  if (random_intercept) {
    zero_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
    false_matrix <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)

    # expand startvalue matrix
    startvalues <- rbind(
      cbind(startvalues, zero_matrix),
      cbind(zero_matrix, zero_matrix)
    )

    # expand free matrix
    free <- rbind(cbind(free, false_matrix), cbind(false_matrix, false_matrix))

    # expand label matrix
    labels <- rbind(cbind(labels, na_matrix), cbind(na_matrix, na_matrix))

    # expand lbound and ubound
    lbound <- rbind(cbind(lbound, na_matrix), cbind(na_matrix, na_matrix))
    ubound <- rbind(cbind(ubound, na_matrix), cbind(na_matrix, na_matrix))
  }

  #### create OpenMx model object ####
  Q <- OpenMx::mxMatrix(
    type = "Full",
    name = "Q",
    nrow = nrow(startvalues),
    ncol = nrow(startvalues),
    free = free,
    values = startvalues,
    labels = labels,
    lbound = lbound,
    ubound = ubound,
    byrow = TRUE
  )

  return(Q)
}
