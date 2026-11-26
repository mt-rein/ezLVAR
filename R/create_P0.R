#' Creates a P0 Matrix for OpenMx
#'
#' @description This is a helper function to create the P0 matrix. The P0 matrix
#'   describes the initial (co)variances in the State Space model.
#'
#' @param step2output An object obtained with the [step2] function.
#' @param random_intercept Logical. If TRUE, the matrices `startvalues`, `free`,
#'   `labels`, `lbound`, and `ubound` are automatically expanded to accommodate
#'   the specification of a random intercept.
#' @param startvalues A square matrix of numeric values that represents the
#'   starting values for each parameter. The number of rows/columns must be
#'   equal to the number of latent factors in the model. If
#'   `random_intercept = TRUE`, `startvalues` must be a list of two square
#'   matrices: the first corresponds to the within-person (co)variances and the
#'   second corresponds to the between-person (co)variances. Optional. If
#'   `NULL`, start values will be computed as the (within-person and
#'   between-person) (co)variances of the factor score variables.
#' @param free A square matrix of TRUE and FALSE values that indicates which
#'   parameters are freely estimated. The number of rows/columns must be equal
#'   to the number of latent factors in the model. If `random_intercept = TRUE`,
#'   `free` must be a list of two square matrices: the first corresponds to the
#'   within-person (co)variances and the second corresponds to the
#'   between-person (co)variances. Optional. If `NULL`, all within- and between-person (co)variances are freely estimated.
#' @param labels A square matrix of strings that indicates the labels for each
#'   parameter (if a value is NA, no label is assigned). The number of
#'   rows/columns must be equal to the number of latent factors in the model. If
#'   `random_intercept = TRUE`, `labels` must be a list of two square matrices:
#'   the first corresponds to the within-person (co)variances and the second
#'   corresponds to the between-person (co)variances. Optional. If `NULL`,
#'   labels will be automatically generated.
#' @param lbound A square matrix of numeric values that indicates the lower
#'   bounds for each parameter (if a value is NA, no bounds are imposed on that
#'   parameter). The number of rows/columns must be equal to the number of
#'   latent factors in the model. If `random_intercept = TRUE`, `lbound` must be
#'   a list of two square matrices: the first corresponds to the within-person
#'   (co)variances and the second corresponds to the between-person
#'   (co)variances. Optional. If `NULL`, no lower bounds are set.
#' @param ubound A square matrix of numeric values that indicates the upper
#'   bounds for each parameter (if a value is NA, no bounds are imposed on that
#'   parameter). The number of rows/columns must be equal to the number of
#'   latent factors in the model. If `random_intercept = TRUE`, `ubound` must be
#'   a list of two square matrices: the first corresponds to the within-person
#'   (co)variances and the second corresponds to the between-person
#'   (co)variances. Optional. If `NULL`, no upper bounds are set.
#'
#' @returns P0 An `OpenMx` matrix object that is entered into [step3()].
#'
#' @export
create_P0 <- function(
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
  if (!is.logical(random_intercept) || is.na(random_intercept)) {
    stop("The random_intercept argument must be TRUE or FALSE.")
  }

  #### no random intercept specification ####
  if (!random_intercept) {
    #### create elements that have not been provided by the user ####
    # create start values matrix ((co)variances of factor scores)
    if (is.null(startvalues)) {
      startvalues <- stats::var(step2output$data[, factors, drop = FALSE], na.rm = TRUE)
    }

    # create matrix that indicates free parameters (all freely estimated)
    if (is.null(free)) {
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }

    # create label matrix
    if (is.null(labels)) {
      labels <- matrix(NA, n_factors, n_factors)
      # fill the matrix with symmetric labels
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
    }

    # create lbound and ubound objects (no bounds)
    if (is.null(lbound)) {
      lbound <- matrix(NA, n_factors, n_factors)
    }
    if (is.null(ubound)) {
      ubound <- matrix(NA, n_factors, n_factors)
    }
  }

  #### with random intercept specification ####
  if (random_intercept) {
    # matrices that will be used to expand start values/labels/free values etc
    # to accommodate random intercept
    zero_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
    false_matrix <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)

    #### create elements that have not been provided by the user ####
    # create start values (pooled within-person covariance matrix and
    # between-person covariance matrix)
    if (is.null(startvalues)) {
      # extract ID variable and factor scores
      id_name <- step2output$other$id_var
      id <- step2output$data[[id_name]]

      X <- step2output$data[, factors, drop = FALSE]

      # use only rows with complete factor scores and non-missing ID
      keep <- stats::complete.cases(X) & !is.na(id)

      id <- id[keep]
      X <- X[keep, , drop = FALSE]

      # pooled within-person covariance matrix:
      # person-mean-center the factor scores
      X_centered <- X

      for (f in factors) {
        person_means <- stats::ave(X[[f]], id, FUN = mean)
        X_centered[[f]] <- X[[f]] - person_means
      }

      X_centered <- as.matrix(X_centered)

      # compute the pooled within-person covariance matrix
      n_obs_total <- nrow(X_centered)
      n_persons <- length(unique(id))

      within_cov <- crossprod(X_centered) / (n_obs_total - n_persons)
      startvalues_w <- within_cov

      # between-person covariance matrix
      # compute person means for each factor:
      person_means <- stats::aggregate(
        X,
        by = list(id = id),
        FUN = mean
      )

      # compute the between-person covariance matrix:
      between_cov <- stats::cov(person_means[, factors, drop = FALSE])
      startvalues_b <- between_cov
    } else {
      startvalues_w <- startvalues[[1]]
      startvalues_b <- startvalues[[2]]
    }

    # create complete start values matrix
    startvalues <- rbind(
      cbind(startvalues_w, zero_matrix),
      cbind(zero_matrix, startvalues_b))

    # create matrix that indicates free parameters (all freely estimated)
    if (is.null(free)) {
      free_w <- free_b <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    } else {
      free_w <- free[[1]]
      free_b <- free[[2]]
    }
    # create complete free matrix:
    free <- rbind(cbind(free_w, false_matrix), cbind(false_matrix, free_b))

    # create label matrix
    if (is.null(labels)) {
      labels_w <- matrix(NA, n_factors, n_factors)
      # fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels_w[i, j] <- paste0(
            "P0_w_",
            factors[min(i, j)],
            "_",
            factors[max(i, j)]
          )
        }
      }

      labels_b <- matrix(NA, n_factors, n_factors)
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels_b[i, j] <- paste0(
            "P0_b_",
            factors[min(i, j)],
            "_",
            factors[max(i, j)]
          )
        }
      }
    } else {
      labels_w <- labels[[1]]
      labels_b <- labels[[2]]
    }

    # create complete labels matrix:
    labels <- rbind(cbind(labels_w, na_matrix), cbind(na_matrix, labels_b))

    # create lbound and ubound objects (no bounds)
    if (is.null(lbound)) {
      lbound_w <- lbound_b <- matrix(NA, n_factors, n_factors)
    } else {
      lbound_w <- lbound[[1]]
      lbound_b <- lbound[[2]]
    }
    if (is.null(ubound)) {
      ubound_w <- ubound_b <- matrix(NA, n_factors, n_factors)
    } else {
      ubound_w <- ubound[[1]]
      ubound_b <- ubound[[2]]
    }

    # create complete lbound and ubound matrices:
    lbound <- rbind(cbind(lbound_w, na_matrix), cbind(na_matrix, lbound_b))
    ubound <- rbind(cbind(ubound_w, na_matrix), cbind(na_matrix, ubound_b))
  }


  #### create OpenMx model object ####
  P0 <- OpenMx::mxMatrix(
    type = "Full",
    name = "P0",
    nrow = nrow(startvalues),
    ncol = nrow(startvalues),
    free = free,
    values = startvalues,
    labels = labels,
    lbound = lbound,
    ubound = ubound
  )

  return(P0)
}
