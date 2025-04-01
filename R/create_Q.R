#' Creates an Innovation Covariance Matrix
#'
#' @description
#' This is a helper function to create the Q matrix, which is required to specify the structural model for `step3()`.
#' The Q matrix describes the innovation (co)variances (i.e., the residuals of the latent variables).
#'
#'
#' @param step2output The output obtained with the `step2()` function.
#' @param startvalues A matrix that represents the starting values. Must have the same dimensions as `free`, `labels`, `lbound` and `ubound` if these are not `NULL`.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`, `labels`, `lbound`, and `ubound` are automatically expanded to accommodate the specification a random intercept.
#' @param free A matrix of TRUE and FALSE values that indicates which parameters are freely estimated. Optional. If NULL, all variances and covariances are freely estimated. If not NULL, the matrix must have the same dimensions as `startvalues`, as well as `labels`, `lbound` and `ubound` if these are not `NULL`.
#' @param labels A matrix of strings that indicates the labels for each parameter. Optional. If NULL, labels will be automatically generated. If not NULL, the matrix must have the same dimensions as `startvalues`, as well as `free`, `lbound` and `ubound` if these are not `NULL`.
#' @param lbound A matrix of numeric values that indicates the lower bounds for each parameter. Optional. If NA, no bounds are imposed.
#' @param ubound A matrix of numeric values that indicates the upper bounds for each parameter. Optional. If NA, no bounds are imposed.
#'
#' @return Q An `OpenMx` matrix object that is passed on to `step3()`.
#' @export
#'
#' @examples
create_Q <- function(step2output, startvalues,
                     random_intercept = FALSE, free = NULL,
                     labels = NULL, lbound = NA, ubound = NA){

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

  if(!random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }

    # create label matrix
    if(is.null(labels)){
      labels <- matrix(NA, n_factors, n_factors)

      # Fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels[i, j] <- paste0("zeta_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
    }

    # create lbound and ubound objects if not specified by user
    if(is.null(lbound)){
      lbound <- NA
    }
    if(is.null(ubound)){
      ubound <- NA
    }

    # create OpenMx model object
    Q <- OpenMx::mxMatrix(type = "Full", name = "Q",
                          nrow = n_factors, ncol = n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

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
      labels <- matrix(NA, n_factors, n_factors)

      # Fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels[i, j] <- paste0("zeta_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
    }
    # expand label matrix with random intercept specification
    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, others),
                    cbind(others, others))

    # expand lbound and ubound matrices with random intercept specification:
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
    others <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(cbind(startvalues, others),
                    cbind(others, others))


    # create OpenMx model object
    Q <- OpenMx::mxMatrix(type = "Full", name = "Q",
                          nrow = 2*n_factors, ncol = 2*n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

  return(Q)
}
