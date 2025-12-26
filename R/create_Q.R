#' Creates an Innovation Covariance Matrix
#'
#' @description
#' This is a helper function to create the Q matrix, which is required to specify the structural model for [step3()].
#' The Q matrix describes the innovation (co)variances (i.e., the residuals of the latent variables).
#'
#'
#' @param step2output The output obtained with the [step2()] function.
#' @param startvalues A matrix that represents the starting values for each parameter. Must have the same dimensions as `free`, `labels`, `lbound` and `ubound` if these are not `NULL`.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`, `labels`, `lbound`, and `ubound` are automatically expanded to accommodate the specification a random intercept.
#' @param free A matrix of TRUE and FALSE values that indicates which parameters are freely estimated. Optional. If NULL, all variances and covariances are freely estimated. If not NULL, the matrix must have the same dimensions as `startvalues`, as well as `labels`, `lbound` and `ubound` if these are not `NULL`.
#' @param labels A matrix of strings that indicates the labels for each parameter. Optional. If NULL, labels will be automatically generated. If not NULL, the matrix must have the same dimensions as `startvalues`, as well as `free`, `lbound` and `ubound` if these are not `NULL`.
#' @param lbound A matrix of numeric values that indicates the lower bounds for each parameter. Optional. If NA, no bounds are imposed.
#' @param ubound A matrix of numeric values that indicates the upper bounds for each parameter. Optional. If NA, no bounds are imposed.
#'
#' @returns Q An `OpenMx` matrix object that is enter into [step3()].
#' @export

create_Q <- function(step2output, startvalues,
                     random_intercept = FALSE, free = NULL,
                     labels = NULL, lbound = NULL, ubound = NULL){

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
    fixed <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, fixed),
                  cbind(fixed, fixed))

    # create label matrix
    if(is.null(labels)){
      labels <- matrix(NA, n_factors, n_factors)

      # fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels[i, j] <- paste0("zeta_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
    }
    # expand label matrix with random intercept specification
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, na_matrix),
                    cbind(na_matrix, na_matrix))

    # expand lbound and ubound matrices with random intercept specification:
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
                    cbind(zero_matrix, zero_matrix))


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
