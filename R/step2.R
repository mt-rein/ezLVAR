#### Step 2 Core function ####

#' Step 2: Obtain Factor Scores
#'
#' @description This function performs step 2 of 3S-LVAR: compute factor scores
#'   and their uncertainty. Its output will be passed on to [step3].
#'
#' @param step1output An object obtained with the [step1] function.
#'
#' @returns An object of class `3slvar_step2`, which is a list comprising the
#'   following elements: \item{data}{The original data set with the appended
#'   factor scores.} \item{lambda_star}{A data frame containing the values of
#'   lambda* per individual.} \item{theta_star}{A data frame containing the
#'   values of theta* per individual.} \item{other}{A list containing
#'   information on the names of the `factors`, `indicators`, and the ID
#'   variable (`id_var`).}
#'
#' @export
step2 <- function(step1output) {
  ## Preparations
  fit_step1 <- step1output$mm_output
  data <- step1output$data
  measurementmodel <- step1output$measurementmodel
  id_var <- step1output$id_var
  unique_ids <- unique(data[[id_var]])
  n_persons <- length(unique_ids)
  indicators <- lavaan::lavaanify(measurementmodel) |> lavaan::lavNames("ov")
  factors <- lavaan::lavaanify(measurementmodel) |> lavaan::lavNames("lv")
  n_factors <- length(factors)
  group_var <- lavaan::lavInspect(fit_step1[[1]], "group") |> as.character()
  group_names <- lavaan::lavInspect(fit_step1[[1]], "group.label")
  if (rlang::is_empty(group_var)) {
    n_groups <- 1
  } else {
    n_groups <- length(group_names)
  }
  n_blocks <- length(measurementmodel)

  #### Preparations ####

  # create lists to store MM parameter values per block (lists with n_blocks
  # elements, where each element is a list with n_groups elements)
  psi_block <-
    alpha_block <-
      lambda_block <-
        theta_block <-
          tau_block <-
            vector(mode = "list", length = n_blocks)

  for (m in 1:n_blocks) {
    # get the estimates of the current block
    est_block <- lavaan::lavInspect(fit_step1[[m]], "est")
    # if there are multiple groups, est_block is a list.
    # if single-group, wrap into a list of length 1:
    if (n_groups == 1) {
      est_block <- list(est_block)
    }

    # extract the parameter estimates per group by subsetting across all groups
    psi_block[[m]] <- lapply(est_block, function(x) x[["psi"]])
    alpha_block[[m]] <- lapply(est_block, function(x) x[["alpha"]])
    lambda_block[[m]] <- lapply(est_block, function(x) x[["lambda"]])
    theta_block[[m]] <- lapply(est_block, function(x) x[["theta"]])
    tau_block[[m]] <- lapply(est_block, function(x) x[["nu"]])
  }

  # create lists to store MM parameter values per group (i.e., combine the
  # block-wise MM parameters per group) -> Lists with n_groups elements, where
  # each element is a list with n_blocks elements
  psi_group <-
    alpha_group <-
      lambda_group <-
        theta_group <-
          tau_group <-
            vector(mode = "list", length = n_groups)
  for (g in 1:n_groups) {
    # put MM matrices of each group in the respective list
    for (m in 1:n_blocks) {
      psi_group[[g]][[m]] <- psi_block[[m]][[g]]
      alpha_group[[g]][[m]] <- alpha_block[[m]][[g]]
      lambda_group[[g]][[m]] <- lambda_block[[m]][[g]]
      theta_group[[g]][[m]] <- theta_block[[m]][[g]]
      tau_group[[g]][[m]] <- tau_block[[m]][[g]]
    }
    # combine the matrices of the different measurement blocks into a single
    # matrix per group
    psi_group[[g]] <- lavaan::lav_matrix_bdiag(psi_group[[g]])
    alpha_group[[g]] <- do.call(rbind, alpha_group[[g]])
    lambda_group[[g]] <- lavaan::lav_matrix_bdiag(lambda_group[[g]])
    theta_group[[g]] <- lavaan::lav_matrix_bdiag(theta_group[[g]])
    tau_group[[g]] <- do.call(rbind, tau_group[[g]])

    # name the matrices' rows and columns
    rownames(psi_group[[g]]) <- colnames(psi_group[[g]]) <- factors
    rownames(lambda_group[[g]]) <- indicators
    colnames(lambda_group[[g]]) <- factors
    rownames(theta_group[[g]]) <- colnames(theta_group[[g]]) <- indicators
  }

  # if there are multiple groups, obtain the correct group names, otherwise just
  # assign "group1".
  if (n_groups > 1) {
    names(psi_group) <-
      names(alpha_group) <-
        names(lambda_group) <-
          names(theta_group) <-
            names(tau_group) <-
              group_names
  } else {
    names(psi_group) <-
      names(alpha_group) <-
        names(lambda_group) <-
          names(theta_group) <-
            names(tau_group) <-
              "group1"
  }

  #### compute factor scores, lambda_star, theta_star ####
  # add factor score variables to data frame:
  full_data <- data
  for (f in factors) {
    full_data[[f]] <- NA_real_
  }

  # create empty matrices with n_persons rows and n_factors columns
  lambda_star <-
    theta_star <-
      matrix(NA, nrow = n_persons, ncol = n_factors + 1) |>
      as.data.frame()
  colnames(lambda_star) <- colnames(theta_star) <- c(id_var, factors)

  for (i in 1:n_persons) {
    id_i <- unique_ids[i]
    data_i <- data[data[[id_var]] == id_i, , drop = FALSE]
    if (n_groups > 1) {
      # get group value of the individual:
      group_i <- stats::na.omit(data_i[[group_var]]) |> unique()
      # throw errors if the grouping variable for an individual is empty or has
      # more than 1 value
      if (length(group_i) > 1) {
        stop(paste0(
          "Individual ",
          id_i,
          " has more than 1 unique value on the grouping variable '",
          group_var,
          "'."
        ))
      }
      if (rlang::is_empty(group_i)) {
        stop(paste0(
          "Individual ",
          id_i,
          " does not have a valid value on the grouping variable '",
          group_var,
          "'."
        ))
      }
    } else {
      # if there is no grouping variable, just assign the value "group1"
      group_i <- "group1"
    }

    # get corresponding MM values:
    psi_i <- psi_group[[group_i]]
    alpha_i <- alpha_group[[group_i]]
    lambda_i <- lambda_group[[group_i]]
    theta_i <- theta_group[[group_i]]
    tau_i <- tau_group[[group_i]]

    # compute sigma, A, and centered indicators
    sigma_i <- lambda_i %*% psi_i %*% t(lambda_i) + theta_i
    A_i <- psi_i %*% t(lambda_i) %*% solve(sigma_i)
    indicators_cent <- sweep(
      as.matrix(data_i[, indicators]),
      2,
      (tau_i + lambda_i %*% alpha_i),
      "-"
    )

    # compute factor scores and add to factor scores list
    fs_i <- sweep(t(A_i %*% t(indicators_cent)), 2, alpha_i, "+")
    rows_i <- which(data[[id_var]] == id_i)
    full_data[rows_i, factors] <- fs_i

    # compute lambda_star and theta_star and add to respective matrix
    lambda_star[i, id_var] <- theta_star[i, id_var] <- id_i
    lambda_star_i <- (A_i %*% lambda_i) |> diag() |> as.numeric()
    lambda_star[i, factors] <- lambda_star_i
    theta_star_i <- (A_i %*% theta_i %*% t(A_i)) |> diag() |> as.numeric()
    theta_star[i, factors] <- theta_star_i
  }

  # assemble output
  output <- list(
    "data" = full_data,
    "lambda_star" = lambda_star,
    "theta_star" = theta_star,
    "other" = list(
      "factors" = factors,
      "indicators" = indicators,
      "id_var" = id_var
    )
  )
  class(output) <- "3slvar_step2"
  return(output)
}

#### summary function ####

#' @param step2output An object obtained with the [step2] function.
#' @param round_digits Numerical. The number of digits to round numeric values
#'   in the output. Default is 3.
#'
#' @rdname step2
#'
#' @export
summary.3slvar_step2 <- function(step2output, round_digits = 3) {
  factors <- step2output$other$factors

  descriptives_lambdastar <- data.frame(
    factor = character(),
    average = numeric(),
    min = numeric(),
    max = numeric()
  )
  for (i in 1:length(factors)) {
    fac <- factors[i]
    average <- mean(step2output$lambda_star[[fac]])
    minimum <- min(step2output$lambda_star[[fac]])
    maximum <- max(step2output$lambda_star[[fac]])
    descriptives_lambdastar[i, 1] <- fac
    descriptives_lambdastar[i, 2:4] <- c(average, minimum, maximum)
  }

  descriptives_thetastar <- data.frame(
    factor = character(),
    average = numeric(),
    min = numeric(),
    max = numeric()
  )
  for (i in 1:length(factors)) {
    fac <- factors[i]
    average <- mean(step2output$theta_star[[fac]])
    minimum <- min(step2output$theta_star[[fac]])
    maximum <- max(step2output$theta_star[[fac]])
    descriptives_thetastar[i, 1] <- fac
    descriptives_thetastar[i, 2:4] <- c(average, minimum, maximum)
  }

  cat("Descriptives for Lambda* (model-based reliability)\n")
  print_rounded(descriptives_lambdastar, digits = round_digits)
  cat("\n")
  cat("Descriptives for Theta* (residual variance in the factor scores):\n")
  print_rounded(descriptives_thetastar, digits = round_digits)

  invisible(step2output)
}
