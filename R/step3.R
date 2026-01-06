#' Step 3: Estimate the Structural Model
#'
#' @description
#' This function performs Step 3 of 3S-LVAR: Estimate the Structural Model using the factor scores as single indicators while accounting for their inherent uncertainty.
#' The estimation is implemented in the State Space Model framework in `OpenMx`.
#'
#' @param step2output The output obtained with the `step2()` function.
#' @param A An `OpenMx` matrix object that describes the regression coefficients in the State Space model. Diagonal entries represent autoregressive effects, and off-diagonal entries represent cross-lagged effects. See [OpenMx::mxExpectationStateSpace]. The helper function [create_A()] is available to create this object.
#' @param Q An `OpenMx` matrix object that describes the innovation covariance matrix. See [OpenMx::mxExpectationStateSpace]. The helper function [create_Q()] is available to create this object.
#' @param B An `OpenMx` matrix object that describes the covariate effects on the latent variables. Optional. If `NULL`, a zero matrix (i.e., no covariate effects) is automatically created.
#' @param D An `OpenMx` matrix object that describes the covariate effects on the observed items. Optional. If `NULL`, a zero matrix (i.e., no covariate effects) is automatically created.
#' @param x0 An `OpenMx` matrix object that describes the initial latent variable vector that initiates the Kalman Filter. Optional. If `NULL`, an object is automatically created where all values are fixed to 0.
#' @param P0 An `OpenMx` matrix object that describes the initial error covariance matrix that initiates the Kalman Filter. Optional. If `NULL`, a matrix will be automatically created. The parameters are freely estimated, with large starting values (100 for the diagonal and 10 for the off-diagonal).
#' @param u An `OpenMx` matrix object that indicates the covariates. Required if `B` or `D` are included. If `NULL`, an object with no covariates is automatically created.
#' @param group_var A string containing the name of the grouping variable for the structural model.
#' @param mixture Logical. If TRUE, a mixture model is applied.
#' @param n_clusters Numerical. The number of clusters (required if `mixture = TRUE`).
#' @param n_starts Numerical. The total number of random starts. Only used if `mixture = TRUE`.
#' @param n_best_starts Numerical. The number of starts that are completed until convergence in Phase 2. Only used if `mixture = TRUE`.
#' @param maxit_phase1 Numerical. The maximum number of iterations in Phase 1. Only used if `mixture = TRUE`.
#' @param maxit_phase2 Numerical. The maximum number of iterations in Phase 2. Only used if `mixture = TRUE`.
#' @param stablecluster_criterion Numerical. If the average change in the posterior probabilities across persons after the E-Step in Phase 1 is smaller than this value, the clustering is considered stable.  Only used if `mixture = TRUE`.
#' @param convergence_criterion Numerical. If the change in the observed data log likelihood after the M-Step is smaller than this value, convergence is achieved. Only used if `mixture = TRUE`.
#' @param GEM_iterations The maximum number of iterations in OpenMx during the M-Step (when updating the model parameters) in Phase 1. Only used if `mixture = TRUE`.
#' @param parallel Logical. If TRUE, parallel processing is used. Only used if `mixture = TRUE`.
#' @param n_cores Numerical. number of cores to use when `parallel = TRUE`. If `NULL`, two cores are kept free for other tasks.
#' @param seeds A numerical vector containing the starting seeds for each random start. Optional. If `NULL`, random seeds are automatically generated.  Only used if `mixture = TRUE`.
#' @param newdata A data frame. Note that the output of `step2()` already contains the data, so using this argument is only required if the data have been manipulated after step 2 (e.g., removing outliers). Optional.
#' @param tryhard Should the algorithm run multiple times to obtain a solution (see [OpenMx::mxTryHard])? Does not work in conjunction with mixture modeling.
#' @param verbose Logical (FALSE by default). If TRUE, print progress messages.
#'
#' @returns A list containing the following elements:
#' @returns `estimates` A vector that includes the parameter estimates. In the case of mixture modeling, a list of vectors (one for each cluster).
#' @returns `model` The fitted model (an OpenMx object). In the case of mixture modeling, a list of models (one for each cluster).
#' @returns `data` The data set which served as input.
#' @returns If mixture modeling is used, the list also contains the following elements:
#' @returns `loglik` The observed data log likelihood of the best model.
#' @returns `posterior_probabilities` A data.frame with the posterior probabilities per person per cluster, and their modal cluster assignment (i.e., for which cluster they have the highest posterior probability).
#' @returns `class_proportions` The class proportions (prior probabilities).
#' @returns `n_nonconverged` The number of starts that did not converge in Phase 2.
#' @returns `seeds` A vector containing the seeds per random start.
#' @returns `best_seed` The seed of the best random start.
#'
#' @export
#' @importFrom OpenMx imxReportProgress

step3 <- function(step2output, A, Q,
                  B = NULL, D = NULL,
                  x0 = NULL, P0 = NULL,
                  u = NULL,
                  group_var = NULL,
                  mixture = FALSE,
                  n_clusters = NULL,
                  n_starts = 25,
                  n_best_starts = 5,
                  maxit_phase1 = 100,
                  maxit_phase2 = 100,
                  stablecluster_criterion = .1,
                  convergence_criterion = 1e-6,
                  GEM_iterations = 3,
                  parallel = FALSE,
                  n_cores = NULL,
                  seeds = NULL,
                  newdata = NULL,
                  tryhard = FALSE,
                  verbose = TRUE) {

  #### Errors and Warnings ####
  if (mixture) {
    if (!is.null(group_var)) {
      stop(glue::glue("You cannot combine multi-group modeling with mixture modeling. ",
                      "Set the 'group_var' argument to NULL or the 'mixture' argument to FALSE."))
    }

    if (is.null(n_clusters)) {
      stop("When using mixture modeling, 'n_clusters' cannot be empty.")
    }

    if (n_clusters < 2) {
      stop("The number of clusters (n_clusters) cannot be smaller than 2.")
    }

    if (!is.null(seeds) && seeds != n_starts) {
      stop("You need to provide as many seeds (length of 'seeds') as random starts (size of 'n_starts').")
    }

    if (parallel && !is.null(n_cores) && n_cores > parallel::detectCores()) {
      stop("The number of cores you specified in the 'n_cores' argument is greater than the number of cores on your machine.")
    }
  }

  #### 1) Preparations ####
  ## extract objects from step 1 output:
  if (is.null(newdata)) {
    data <- step2output$data
  } else {
    data <- newdata
  }

  id_var <- step2output$other$id_var
  unique_ids <- data[[id_var]] |> unique()
  n_persons <- length(unique_ids)
  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  factors <- step2output$other$factors
  factors_ind <- paste0(factors, "_ind")
  # names of factor score variables (single indicators)
  n_factors <- length(factors)

  # rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |>
    dplyr::rename_with(~ factors_ind, tidyselect::all_of(factors))

  #### 2) create OpenMx matrices ####
  # number of latent variables:
  xdim <- nrow(A$free)
  # exogenous covariates (ignored so far):
  udim <- 1
  # number of indicators (i.e., factor score variables)
  ydim <- n_factors


  # B matrix (= exogenous covariates on latent constructs)
  if (is.null(B)) {
    B <- OpenMx::mxMatrix("Zero", name = "B",
                          nrow = xdim, ncol = udim)
  }

  # D matrix (= exogenous covariates on observed variables)
  if (is.null(D)) {
    D <- OpenMx::mxMatrix("Zero", name = "D",
                          nrow = ydim, ncol = udim)
  }


  # x0 and P0 (= initial values and (co)variances of the latent constructs)
  if (is.null(x0)) {
    if (xdim > ydim) {
      labels <- c(paste0("ini_", factors),
                  paste0("m_intercept_", factors))
    } else {
      labels <- paste0("ini_", factors)
    }

    x0 <- OpenMx::mxMatrix("Full", name = "x0",
                           nrow = xdim, ncol = 1,
                           free = FALSE,
                           values = 0,
                           labels = labels)
  }

  if (is.null(P0)) {
    # create matrix that indiciates free parameters
    free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)

    # create label matrix
    labels <- matrix(NA, n_factors, n_factors)
    # Fill the matrix with symmetric labels
    for (i in 1:n_factors) {
      for (j in 1:n_factors) {
        labels[i, j] <- paste0("P0_", factors[min(i, j)], "_", factors[max(i, j)])
      }
    }

    # create startvalues
    startvalues <- matrix(10, nrow = n_factors, ncol = n_factors)
    diag(startvalues) <- 100

    # if random intercept is included, expand these matrices
    if (xdim > ydim) {
      others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
      free <- rbind(cbind(free, others),
                    cbind(others, free))

      others <- matrix(NA, nrow = n_factors, ncol = n_factors)
      intercepts <- matrix(NA, n_factors, n_factors)
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          intercepts[i, j] <- paste0("P0_icp_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
      labels <- rbind(cbind(labels, others),
                      cbind(others, intercepts))

      others <- matrix(0, nrow = n_factors, ncol = n_factors)
      startvalues <- rbind(cbind(startvalues, others),
                           cbind(others, startvalues))
    }

    P0 <- OpenMx::mxMatrix(type = "Full", name = "P0",
                           nrow = xdim, ncol = xdim,
                           free = free,
                           values = startvalues,
                           labels = labels,
                           byrow = TRUE)
  }

  # u (= covariates)
  if (is.null(u)) {
    u <- OpenMx::mxMatrix("Zero", name = "u", nrow = udim, ncol = 1)
  }

  #### 3) create person-models ####
  # create a list of models (one for each individual):
  personmodelnames <- paste0("id_", unique_ids)

  personmodel_list <- vector(mode = "list", length = n_persons)
  names(personmodel_list) <- unique_ids
  for (i in 1:n_persons) {
    # get person's id, data, lambda_star, theta_star
    id_i <- unique_ids[i]
    data_i <- data[data[[id_var]] == id_i, ] |> as.data.frame()
    lambda_star_i <- lambda_star[lambda_star[[id_var]] == id_i, factors] |>
      as.numeric()
    theta_star_i <- theta_star[theta_star[[id_var]] == id_i, factors] |>
      as.numeric()
    # fix the loadings (C matrix) to lambda_star
    values <- matrix(0, nrow = n_factors, ncol = n_factors)
    diag(values) <- lambda_star_i
    C_names <- list(factors_ind, factors)
    if (xdim > ydim) {
      values <- cbind(values, values)
      C_names <- list(factors_ind, c(paste0(factors),
                                     paste0("intercept_", factors)))
    }
    C <- OpenMx::mxMatrix("Full", name = "C",
                          nrow = ydim, ncol = xdim,
                          free = FALSE,
                          values = values,
                          labels = NA,
                          byrow = TRUE,
                          dimnames = C_names
    )

    # fix the residual variances (R matrix) to theta_star
    # R matrix (= measurement noise)
    R <- OpenMx::mxMatrix("Diag", name = "R",
                          nrow = ydim, ncol = ydim,
                          free = FALSE,
                          values = theta_star_i,
                          labels = NA
    )

    # create the model
    modelname <- personmodelnames[i]
    model_i <- OpenMx::mxModel(name = modelname,
                               A, B, C, D,
                               Q, R, x0, P0,
                               u,
                               OpenMx::mxExpectationStateSpace('A', 'B', 'C', 'D',
                                                               'Q', 'R', 'x0', 'P0',
                                                               'u'),
                               OpenMx::mxFitFunctionML(),
                               OpenMx::mxData(data_i[, factors_ind, drop = FALSE],
                                              'raw'))

    # if there is a grouping variable, relabel parameters
    # (individuals in the same group get the same labels)
    # (individuals in different groups get different labels)
    if (!is.null(group_var)) {
      # get group value of the individual:
      group_i <- stats::na.omit(data_i[[group_var]]) |> unique()
      # throw errors if the grouping variable for an individual is empty or has more than 1 value
      if (length(group_i) > 1) {
        stop(glue::glue("Individual {id_i} has more than 1 unique value on the grouping variable '{group_var}'."))
      }
      if (rlang::is_empty(group_i)) {
        stop(glue::glue("Individual {id_i} does not have a valid value on the grouping variable '{group_var}'."))
      }

      model_i <- OpenMx::omxSetParameters(model_i,
                                          labels = names(stats::coef(model_i)),
                                          newlabels = paste0(names(stats::coef(model_i)),
                                                             "_", group_i))
    }

    # add person-model to the list
    personmodel_list[[i]] <- model_i
  }

  #### 4) single- and multi-group modeling ####
  # combine person-models in a multi-group (i.e., multi-subject) model
  if (!mixture) {
    fullmodel <- OpenMx::mxModel("fullmodel",
                                 personmodel_list,
                                 OpenMx::mxFitFunctionMultigroup(personmodelnames))
    # fit the model
    if (tryhard) {
      fullmodelr <- OpenMx::mxTryHard(fullmodel, silent = !verbose)
    } else {
      fullmodelr <- OpenMx::mxRun(fullmodel, silent = !verbose)
    }

    estimates <- stats::coef(fullmodelr)

    ## build output
    output <- list("estimates" = estimates,
                   "model" = fullmodelr,
                   "data" = data)
  }

  #### 5) mixture modeling ####
  # The estimation happens in two parts:
  # First, the algorithm aims for a stable clustering with an imprecise M-Step
  # after a stable clustering is achieved, a final precise M-Step is performed
  # Then, the n_best_starts best starts are completed
  if (mixture) {
    # provide seeds for the multiple starts (for replicability)
    if (is.null(seeds)) {
      seeds <- sample.int(.Machine$integer.max, n_starts)
    }
    input_phase1 <- lapply(seeds, function(s) list(seed = s))

    ## run Phase 1 (Generalized EM for all starts):
    if (verbose) {
      message(glue::glue("==== Starting Phase 1 with {n_starts} starts. ===="))
    }

    if(parallel) {
      # add number of cores if not specified:
      if (is.null(n_cores)) {
        n_cores = parallel::detectCores() - 2
        if(n_cores < 2) {
          stop("Your machine does not have enough cores to keep 2 free. Please use sequential processing or specify a number of cores using the 'n_cores' argument.")
        }
      }

      # create cluster:
      backend <- parabar::start_backend(n_cores, cluster_type = "psock", backend_type = "async")

      # add variables to cluster:
      parabar::export(backend,
                      variables = c(
                        "personmodel_list", "n_clusters", "n_factors", "n_persons",
                        "stablecluster_criterion", "convergence_criterion",
                        "GEM_iterations", "maxit_phase1"
                      ),
                      envir = environment())

      # run random starts in parallel:
      results_phase1 <- parabar::par_lapply(backend,
                                            input_phase1,
                                            run_start,
                                            personmodel_list = personmodel_list,
                                            n_clusters = n_clusters,
                                            n_factors = n_factors,
                                            n_persons = n_persons,
                                            phase = 1,
                                            stablecluster_criterion = stablecluster_criterion,
                                            convergence_criterion = convergence_criterion,
                                            GEM_iterations = GEM_iterations,
                                            maxit = maxit_phase1,
                                            verbose = FALSE)

      # close cluster:
      parabar::stop_backend(backend)
    } else {
      results_phase1 <- lapply(input_phase1,
                               run_start,
                               personmodel_list = personmodel_list,
                               n_clusters = n_clusters,
                               n_factors = n_factors,
                               n_persons = n_persons,
                               phase = 1,
                               stablecluster_criterion = stablecluster_criterion,
                               convergence_criterion = convergence_criterion,
                               GEM_iterations = GEM_iterations,
                               maxit = maxit_phase1,
                               verbose = verbose)
    }

    if (verbose) {
      message(glue::glue("==== Phase 1 finished.",
                         "Proceeding to Phase 2 with {n_best_starts} starts. ===="))
    }

    # run Phase 2 (full EM for the best starts):
    top_starts <- results_phase1 |>
      purrr::map_dbl(~ .x$observed_data_LL) |>
      order(decreasing = TRUE) |>
      utils::head(n_best_starts)

    input_phase2 <- results_phase1[top_starts]

    if(parallel) {
      # create cluster:
      backend <- parabar::start_backend(n_cores, cluster_type = "psock", backend_type = "async")

      # add variables to cluster:
      parabar::export(backend,
                      variables = c(
                        "personmodel_list", "n_clusters", "n_factors", "n_persons",
                        "stablecluster_criterion", "convergence_criterion",
                        "GEM_iterations", "maxit_phase1"
                      ),
                      envir = environment())

      # run random starts in parallel:
      results_phase2 <- parabar::par_lapply(backend,
                                            input_phase2,
                                            run_start,
                                            personmodel_list = personmodel_list,
                                            n_clusters = n_clusters,
                                            n_factors = n_factors,
                                            n_persons = n_persons,
                                            phase = 2,
                                            stablecluster_criterion = stablecluster_criterion,
                                            convergence_criterion = convergence_criterion,
                                            GEM_iterations = GEM_iterations,
                                            maxit = maxit_phase2,
                                            verbose = FALSE)

      # close cluster:
      parabar::stop_backend(backend)
    } else {
      results_phase2 <- lapply(input_phase2,
                               run_start,
                               personmodel_list = personmodel_list,
                               n_clusters = n_clusters,
                               n_factors = n_factors,
                               n_persons = n_persons,
                               phase = 2,
                               stablecluster_criterion = stablecluster_criterion,
                               convergence_criterion = convergence_criterion,
                               GEM_iterations = GEM_iterations,
                               maxit = maxit_phase2,
                               verbose = verbose)
    }

    if (verbose) {
      message(glue::glue("==== Phase 2 finished. ===="))
    }

    # select best model:
    best_start <- results_phase2 |>
      purrr::map_dbl(~ .x$observed_data_LL) |>
      order(decreasing = TRUE) |>
      utils::head(1)

    selected_start <- results_phase2[[best_start]]


    ## build output
    cluster_names <- paste0("cluster", 1:n_clusters)
    estimates <- selected_start$clustermodels_run |>
      purrr::map(stats::coef)
    names(estimates) <- cluster_names
    best_model <- selected_start$clustermodels_run
    loglik <- selected_start$observed_data_LL
    post <- as.data.frame(selected_start$post) |> round(3)
    colnames(post) <- cluster_names
    post$modal <- apply(post, 1, function(x) names(x)[which.max(x)])
    post[[id_var]] <- unique_ids
    post <- post[, c(id_var, cluster_names, "modal")]
    class_proportions <- selected_start$class_proportions |> round(3)
    names(class_proportions) <- cluster_names
    convergence_counter <- results_phase2 |>
      purrr::map_lgl(~ .x$converged)
    n_nonconverged <- sum(!convergence_counter)
    best_seed <- selected_start$seed
    output <- list("estimates" = estimates,
                   "model" = best_model,
                   "data" = data,
                   "loglik" = loglik,
                   "posterior_probabilities" = post,
                   "class_proportions" = class_proportions,
                   "n_nonconverged" = n_nonconverged,
                   "seeds" = seeds,
                   "best_seed" = best_seed)
  }

  return(output)
}
