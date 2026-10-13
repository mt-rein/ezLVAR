#' Creates a Regression Effects Matrix
#'
#' @description This is a helper function to create the P0 matrix. The P0 matrix
#'   describes the initial (co)variances in the State Space model.
#'
#' @param step2output An object obtained with the [step2] function.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`,
#'   `labels`, `lbound`, and `ubound` are expanded to accommodate the
#'   specification of a random intercept.
#' @param startvalues A square matrix of numeric values that represents the
#'   starting values for each parameter. The number of rows/columns must be
#'   equal to the number of latent factors in the model. Optional. If `NULL`,
#'   the starting values are 1000 on the diagonal and 100 on the off-diagonal.
#'
#' @returns P0 An `OpenMx` matrix object that is entered into [step3()].
#'
#' @export
create_P0 <- function(
  step2output,
  random_intercept = FALSE,
  startvalues = NULL
) {
  factors <- step2output$other$factors
  n_factors <- length(factors)

  #### Errors ####
  if (!is.logical(random_intercept)) {
    stop("The random_intercept argument must be TRUE or FALSE.")
  }

  #### create elements that have not been provided by the user ####
  # create matrix that indicates free parameters
  free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
  # create label matrix
  if (is.null(labels)) {
    labels <- outer(factors, factors, FUN = function(i, j) {
      paste0("phi_", i, "_", j)
    })
  }

  # create label matrix
  labels <- matrix(NA, n_factors, n_factors)
  # Fill the matrix with symmetric labels
  for (i in 1:n_factors) {
    for (j in 1:n_factors) {
      labels[i, j] <- paste0(
        "P0_",
        factors[min(i, j)],
        "_",
        factors[max(i, j)]
      )
    }
  }

  # create matrix with start values (1000 on diagonal, 100 on off-diagonal)
  if (is.null(startvalues)) {
    startvalues <- matrix(100, nrow = n_factors, ncol = n_factors)
    diag(startvalues) <- 1000
  }

  #### if random intercept has been requested, expand the matrices ####
  if (random_intercept) {
    others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, others), cbind(others, free))

    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    intercepts <- matrix(NA, n_factors, n_factors)
    for (i in 1:n_factors) {
      for (j in 1:n_factors) {
        intercepts[i, j] <- paste0(
          "P0_icp_",
          factors[min(i, j)],
          "_",
          factors[max(i, j)]
        )
      }
    }
    labels <- rbind(cbind(labels, others), cbind(others, intercepts))

    others <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(
      cbind(startvalues, others),
      cbind(others, startvalues)
    )
  }

  #### create OpenMx model object ####
  A <- OpenMx::mxMatrix(
    type = "Full",
    name = "A",
    nrow = nrow(startvalues),
    ncol = nrow(startvalues),
    free = free,
    values = startvalues,
    labels = labels,
    lbound = lbound,
    ubound = ubound,
    byrow = TRUE
  )

  return(A)
}
