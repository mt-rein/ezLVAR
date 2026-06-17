#' Creates a Regression Effects Matrix
#'
#' @description This is a helper function to create the A matrix, which is
#'   required to specify the structural model for [step3]. The A matrix
#'   describes the regression coefficients in the State Space model. Diagonal
#'   entries represent autoregressive effects, and off-diagonal entries
#'   represent cross-lagged effects.
#'
#' @param step2output An object obtained with the [step2] function.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`,
#'   `labels`, `lbound`, and `ubound` are expanded to accommodate the
#'   specification of a random intercept.
#' @param startvalues A square matrix of numeric values that represents the
#'   starting values for each parameter. The number of rows/columns must be
#'   equal to the number of latent factors in the model. Optional. If `NULL`,
#'   the starting values for all regression effects are set to 0.
#' @param free A square matrix of TRUE and FALSE values that indicates which
#'   parameters are freely estimated. The number of rows/columns must be equal
#'   to the number of latent factors in the model. Optional. If `NULL`, all
#'   regression effects are freely estimated.
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
#' @returns A An `OpenMx` matrix object that is entered into [step3()].
#'
#' @export
create_A <- function(
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

  #### If there is no random intercept ####
  if (!random_intercept) {
    # create matrix with start values (all values equal to 0)
    if (is.null(startvalues)) {
      startvalues <- matrix(0, nrow = n_factors, ncol = n_factors)
    }

    # create matrix that indicates free parameters
    if (is.null(free)) {
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }

    # create label matrix
    if (is.null(labels)) {
      labels <- outer(factors, factors, FUN = function(i, j) {
        paste0("phi_", i, "_", j)
      })
    }

    # create lbound and ubound objects if not specified by user
    if (is.null(lbound)) {
      lbound <- NA
    }
    if (is.null(ubound)) {
      ubound <- NA
    }
  }

  #### If there is a random intercept ####
  if (random_intercept) {
    # create matrix with start values (all values equal to 0)
    if (is.null(startvalues)) {
      startvalues <- matrix(0, nrow = n_factors, ncol = n_factors)
    }

    # expand startvalue matrix with random intercept specification
    zero_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(
      cbind(startvalues, zero_matrix),
      cbind(zero_matrix, diag(n_factors))
    )

    # create matrix that indicates free parameters
    if (is.null(free)) {
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }
    # expand free matrix with random intercept specification
    fixed <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, fixed), cbind(fixed, fixed))

    # create label matrix
    if (is.null(labels)) {
      labels <- outer(factors, factors, FUN = function(i, j) {
        paste0("phi_", i, "_", j)
      })
    }
    # expand label matrix with random intercept specification
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, na_matrix), cbind(na_matrix, na_matrix))

    # expand lbound and ubound matrices with random intercept specification
    if (is.null(lbound)) {
      lbound <- NA
    } else {
      lbound <- rbind(cbind(lbound, na_matrix), cbind(na_matrix, na_matrix))
    }
    if (is.null(ubound)) {
      ubound <- NA
    } else {
      ubound <- rbind(cbind(ubound, na_matrix), cbind(na_matrix, na_matrix))
    }

    # create OpenMx model object
    A <- OpenMx::mxMatrix(
      type = "Full",
      name = "A",
      nrow = 2 * n_factors,
      ncol = 2 * n_factors,
      free = free,
      values = startvalues,
      labels = labels,
      lbound = lbound,
      ubound = ubound,
      byrow = TRUE
    )
  }

  # create OpenMx model object
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
