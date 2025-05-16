#' Creates a Regression Effects Matrix
#'
#' @description
#' This is a helper function to create the A matrix, which is required to specify the structural model for [step3()].
#' The A matrix describes the regression coefficients in the State Space model. Diagonal entries represent autoregressive effects, and off-diagonal entries represent cross-lagged effects.
#'
#' @param step2output The output obtained with the [step2()] function.
#' @param startvalues A matrix that represents the starting values. Must have the same dimensions as `free`, `labels`, `lbound` and `ubound` if these are not `NULL`.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`, `labels`, `lbound`, and `ubound` are expanded to accommodate the specification a random intercept.
#' @param free A matrix of TRUE and FALSE values that indicates which parameters are freely estimated. Optional. If NULL, all regression effects are freely estimated. If not NULL, the matrix must have the same dimensions as `startvalues`.
#' @param labels A matrix of strings that indicates the labels for each parameter. Optional. If NULL, labels will be automatically generated. If not NULL, the matrix must have the same dimensions as `startvalues`.
#' @param lbound A matrix of numeric values that indicates the lower bounds for each parameter. Optional. If NULL, no bounds are imposed. If not NULL, the matrix must have the same dimensions as `startvalues`.
#' @param ubound A matrix of numeric values that indicates the upper bounds for each parameter. Optional. If NULL, no bounds are imposed. If not NULL, the matrix must have the same dimensions as `startvalues`.
#'
#' @return A An `OpenMx` matrix object that is enter into [step3()].
#' @export

create_A <- function(step2output, startvalues,
                     random_intercept = FALSE, free = NULL,
                     labels = NULL, lbound = NULL, ubound = NULL){

  # Checks:
  if(!is.null(free)){
    if(dim(startvalues) != dim(free)){
      stop("The free matrix must have the same dimensions as the startvalues matrix.")
    }
  }

  if(!is.null(labels)){
    if(dim(startvalues) != dim(labels)){
      stop("The labels matrix must have the same dimensions as the startvalues matrix.")
    }
  }

  if(!is.null(lbound)){
    if(dim(startvalues) != dim(lbound)){
      stop("The lbound matrix must have the same dimensions as the startvalues matrix.")
    }
  }

  if(!is.null(ubound)){
    if(dim(startvalues) != dim(ubound)){
      stop("The ubound matrix must have the same dimensions as the startvalues matrix.")
    }
  }

  if(!is.logical(random_intercept)){
    stop("The random_intercept argument must be TRUE or FALSE.")
  }


  # extract information on factors:
  factors <- step2output$other$factors
  n_factors <- length(factors)

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

      # create OpenMx model object
      A <- OpenMx::mxMatrix(type = "Full", name = "A",
                            nrow = n_factors, ncol = n_factors,
                            free = free,
                            values = startvalues,
                            labels = labels,
                            lbound = lbound,
                            ubound = ubound,
                            byrow = TRUE)
  }

  #### If there is a random intercept:
  if(random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }
    # expand free matrix with random intercept specification
    others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, others),
                  cbind(others, others))

    # create label matrix
    if(is.null(labels)){
      labels <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
    }
    # expand label matrix with random intercept specification
    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, others),
                    cbind(others, others))

    # expand lbound and ubound matrices with random intercept specification
    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    if(is.null(lbound)){
      lbound <- NA
    } else {
      lbound <- rbind(cbind(lbound, others),
                      cbind(others, others))
    }
    if(is.null(ubound)){
      ubound <- NA
    } else {
      ubound <- rbind(cbind(ubound, others),
                      cbind(others, others))
    }

    # expand startvalue matrix with random intercept specification
    startvalues <- rbind(cbind(startvalues, matrix(0, nrow = n_factors, ncol = n_factors)),
                         cbind(matrix(0, nrow = n_factors, ncol = n_factors), diag(n_factors)))


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

  return(A)
}
