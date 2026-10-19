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
#'   equal to the number of latent factors in the model. If `random_intercept =
#'   TRUE`, the start values must be a list of square matrices. The first
#'   element of the list corresponds to the upper left matrix of P0, the second
#'   element to the lower left matrix of P0. Optional. If `NULL`, start values
#'   will be computed as the (co)variances of the factor score variables. If
#'   `random_intercept = TRUE`, the starting values are set to the pooled
#'   within-person (co)variances and the between-person (co)variances of the
#'   factor score variables.
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
  labels <- outer(factors, factors, FUN = function(i, j) {
    paste0("phi_", i, "_", j)
  })

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
    if (!random_intercept) {
      startvalues <- stats::var(step2output$data[, factors], na.rm = TRUE)
    } else {
      # Extract ID and factor scores
      id_name <- step2output$other$id_var
      id <- step2output$data[[id_name]]

      X <- step2output$data[, factors, drop = FALSE]

      # Use only rows with complete factor scores and non-missing ID
      keep <- complete.cases(X) & !is.na(id)

      id <- id[keep]
      X <- X[keep, , drop = FALSE]

      # Number of observations per person
      n_per_person <- table(id)

      # ------------------------------------------------------------------
      #
      # 1. Within-person covariance/variance
      #
      # ------------------------------------------------------------------

      # For within-person covariance, persons need at least 2 observations
      valid_within_ids <- names(n_per_person)[n_per_person >= 2L]
      keep_within <- id %in% valid_within_ids

      id_within <- id[keep_within]
      X_within <- X[keep_within, , drop = FALSE]

      # Person-center the factor scores
      X_centered <- X_within

      for (v in factors) {
        person_means <- ave(X_within[[v]], id_within, FUN = mean)
        X_centered[[v]] <- X_within[[v]] - person_means
      }

      X_centered <- as.matrix(X_centered)

      # Pooled within-person covariance matrix
      n_obs_within <- nrow(X_centered)
      n_persons_within <- length(unique(id_within))

      within_cov <- crossprod(X_centered) / (n_obs_within - n_persons_within)
      startvalues_w <- within_cov

      #
      # ------------------------------------------------------------------
      #
      # 2. Between-person covariance/variance
      #
      # ------------------------------------------------------------------

      # Compute person means for each factor
      person_means <- aggregate(
        X,
        by = list(id = id),
        FUN = mean
      )

      # Remove the ID column
      X_means <- person_means[, factors, drop = FALSE]

      # Between-person covariance matrix
      between_cov <- stats::cov(X_means)

      startvalues_b <- between_cov

      startvalues <- list(startvalues_w, startvalues_b)
    }
  }

  #### if random intercept has been requested, expand the matrices ####
  if (random_intercept) {
    zero_matrix <- matrix(0, nrow = n_factors, ncol = n_factors)
    false_matrix <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    na_matrix <- matrix(NA, nrow = n_factors, ncol = n_factors)

    free <- rbind(cbind(free, false_matrix), cbind(false_matrix, free))

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
    labels <- rbind(cbind(labels, na_matrix), cbind(na_matrix, intercepts))

    others <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(
      cbind(startvalues[[1]], zero_matrix),
      cbind(zero_matrix, startvalues[[2]])
    )
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
    byrow = TRUE
  )

  return(P0)
}
