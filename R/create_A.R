#' Creates a Regression Effects Matrix
#'
#' @description
#' This is a helper function to create the A matrix, which is required to specify the structural model for [step3()].
#' The A matrix describes the regression coefficients in the State Space model. Diagonal entries represent autoregressive effects, and off-diagonal entries represent cross-lagged effects.
#'
#' @param step2output The output obtained with the [step2()] function.
#' @param startvalues A square matrix that represents the starting values for each parameter. The number of rows/columns must be equal to the number of latent factors in the model.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`, `labels`, `lbound`, and `ubound` are expanded to accommodate the specification of a random intercept.
#' @param free A matrix of TRUE and FALSE values that indicates which parameters are freely estimated. Optional. If `NULL`, all regression effects are freely estimated. If not `NULL`, the matrix must have the same dimensions as `startvalues`.
#' @param labels A matrix of strings that indicates the labels for each parameter. Optional. If `NULL`, labels will be automatically generated. If not `NULL`, the matrix must have the same dimensions as `startvalues`.
#' @param lbound A matrix of numeric values that indicates the lower bounds for each parameter (if a value is NA, no bounds are imposed on that parameter). Optional. If `NULL`, no bounds are imposed. If not `NULL`, the matrix must have the same dimensions as `startvalues`.
#' @param ubound A matrix of numeric values that indicates the upper bounds for each parameter (if a value is NA, no bounds are imposed on that parameter). Optional. If `NULL`, no bounds are imposed. If not `NULL`, the matrix must have the same dimensions as `startvalues`.
#'
#' @returns A An `OpenMx` matrix object that is entered into [step3()].
#' @export

create_A <- function(step2output, startvalues,
                     random_intercept = FALSE, free = NULL,
                     labels = NULL, lbound = NULL, ubound = NULL) {

  # checks:
  matrices <- list(free = free, labels = labels, lbound = lbound, ubound = ubound)

  if(nrow(startvalues) != ncol(startvalues)) {
    stop("The 'startvalues' matrix must be a square matrix.")
  }

  # Loop through each and check dimensions
  for (name in names(matrices)) {
    if (!is.null(matrices[[name]]) && any(dim(startvalues) != dim(matrices[[name]]))) {
      stop(glue::glue("The '{name}' matrix must have the same dimensions as the 'startvalues' matrix."))
    }
  }


  if(!is.logical(random_intercept)){
    stop("The random_intercept argument must be TRUE or FALSE.")
  }


  # extract information on factors:
  factors <- step2output$other$factors
  n_factors <- length(factors)

  if (nrow(startvalues) != n_factors) {
    stop("The number of rows/columns of the startvalues matrix must be equal to the number of latent factors in the model.")
  }

  #### If there is no random intercept:
  if(!random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }

    # create label matrix
    if(is.null(labels)){
      labels <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
    }

    # create lbound and ubound objects if not specified by user
    if(is.null(lbound)){
      lbound <- NA
    }
    if(is.null(ubound)){
      ubound <- NA
    }
  }

  #### If there is a random intercept:
  if(random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }
    # expand free matrix with random intercept specification
    fixed <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, fixed),
                  cbind(fixed, fixed))

    # create label matrix
    if(is.null(labels)){
      labels <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
    }
    # expand label matrix with random intercept specification
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, na_matrix),
                    cbind(na_matrix, na_matrix))

    # expand lbound and ubound matrices with random intercept specification
    if(is.null(lbound)){
      lbound <- NA
    } else {
      lbound <- rbind(cbind(lbound, na_matrix),
                      cbind(na_matrix, na_matrix))
    }
    if(is.null(ubound)){
      ubound <- NA
    } else {
      ubound <- rbind(cbind(ubound, na_matrix),
                      cbind(na_matrix, na_matrix))
    }

    # expand startvalue matrix with random intercept specification
    zero_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(cbind(startvalues, zero_matrix),
                         cbind(zero_matrix, diag(n_factors)))


    # create OpenMx model object
    A <- OpenMx::mxMatrix(type = "Full", name = "A",
                          nrow = 2*n_factors, ncol = 2*n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

  # create OpenMx model object
  A <- OpenMx::mxMatrix(type = "Full", name = "A",
                        nrow = nrow(startvalues), ncol = nrow(startvalues),
                        free = free,
                        values = startvalues,
                        labels = labels,
                        lbound = lbound,
                        ubound = ubound,
                        byrow = TRUE)

  return(A)
}
