# Internal functions to run a single start of the EM algorithm in step3

#### run_start (internal helper) ####
# perform a single random start in either Phase 1 or Phase 2
run_start <- function(
  input_list,
  personmodel_list,
  n_clusters,
  n_factors,
  n_persons,
  phase,
  stablecluster_criterion,
  convergence_criterion,
  GEM_iterations,
  maxit,
  verbose
) {
  # Preparations ----
  # set seed:
  seed <- input_list$seed
  set.seed(seed)
  # in Phase 1, create random hard cluster assignment:
  if (phase == 1) {
    post <- matrix(0, nrow = n_persons, ncol = n_clusters)
    for (person in 1:n_persons) {
      if (person <= n_clusters) {
        # assign the first persons to certain clusters to avoid empty clusters
        # (first person in first cluster, second person in second cluster, ...)
        # until every cluster has one individual
        post[person, person] <- 1
      } else {
        post[person, sample(1:n_clusters, 1)] <- 1
      }
    }

    observed_data_LL0 <- -Inf
  }

  # in Phase 2, extract results from Phase 1
  if (phase == 2) {
    observed_data_LL0 <- input_list$observed_data_LL
    class_proportions <- input_list$class_proportions
    personLL <- input_list$personLL
    clustermodels_run <- input_list$clustermodels_run
  }

  # create list with names of objective functions (needed in M-Step)
  objectives <- sapply(personmodel_list, function(x) x@name) |>
    as.character() |>
    paste0(".objective")

  stable_clustering <- FALSE

  # EM Loop ----

  # loop over iterations until stable clustering is achieved or max iterations
  # are reached
  for (it in 1:maxit) {
    #### E-step: update class membership and class proportions ####
    if (phase == 1 && it > 1) {
      post0 <- post
      # update posteriors:
      post <- EStep(
        pi_ks = class_proportions,
        ngroup = n_persons,
        nclus = n_clusters,
        loglik = personLL
      )

      # compute (absolute) differences in posteriors
      delta <- abs(post - post0)
      person_delta <- rowSums(delta)

      # check if stable clustering has been achieved:
      if (mean(person_delta) < stablecluster_criterion) {
        stable_clustering <- TRUE
      }
    }

    if (phase == 2) {
      # update posteriors:
      post <- EStep(
        pi_ks = class_proportions,
        ngroup = n_persons,
        nclus = n_clusters,
        loglik = personLL
      )
    }

    # compute class proportions:
    class_proportions <- colMeans(post)

    # M-step: fitting SSM model and update parameter estimates ----
    # get start values from the coefficients of previous iteration:
    if (phase == 2 | it > 1) {
      startvalues <- clustermodels_run |>
        purrr::map(OpenMx::omxGetParameters)
    } else {
      # in the first iteration, generate random start values in each cluster
      startvalues <- generate_startvalues(
        n_clusters = n_clusters,
        n_factors = n_factors,
        personmodel_list = personmodel_list
      )
    }

    # create one multi-group (i.e., multi-subject) model per cluster
    # person-models are weighted by their posterior probabilities
    clustermodels <- vector(mode = "list", length = n_clusters)
    for (k in 1:n_clusters) {
      clustername <- paste0("model_k", k)
      weighted_objectives <- paste(post[, k], "*", objectives, collapse = " + ")
      model_k <- OpenMx::mxModel(
        clustername,
        personmodel_list,
        OpenMx::mxAlgebraFromString(weighted_objectives, name = "weightedfit"),
        OpenMx::mxFitFunctionAlgebra("weightedfit")
      )
      # add start values:
      model_k <- OpenMx::omxSetParameters(
        model_k,
        labels = names(startvalues[[k]]),
        values = startvalues[[k]]
      )

      # change options (skip computation of hessian and SEs)
      model_k <- OpenMx::mxOption(
        model_k,
        key = "Calculate Hessian",
        value = "No"
      )
      model_k <- OpenMx::mxOption(
        model_k,
        key = "Standard Errors",
        value = "No"
      )
      # if stable clustering has not yet been achieved in phase 1, reduce
      # iterations
      if (phase == 1 && !stable_clustering) {
        model_k <- OpenMx::mxOption(
          model_k,
          key = "Major iterations",
          value = GEM_iterations
        )
      }

      clustermodels[[k]] <- model_k
    }
    names(clustermodels) <- paste0("model_k", 1:n_clusters)

    # run the models
    clustermodels_run <- clustermodels |>
      purrr::map(OpenMx::mxRun, silent = !verbose, suppressWarnings = TRUE)

    # obtain person-wise LL in a n_persons x n_clusters matrix:
    personLL <- clustermodels_run |>
      # loop over all clusters
      purrr::map(function(clustermodel) {
        # loop over all persons per cluster and extract the person-wise fit
        # function value
        purrr::map_dbl(clustermodel$submodels, function(personmodel) {
          personmodel$fitfunction$result
        })
      }) |>
      simplify2array()
    personLL <- personLL / (-2)
    # openMx gives the minus2LL, so we divide by minus 2

    # compute observed-data log likelihood from person-wise LL:
    observed_data_LL <- compute_observed_data_LL(
      personLL = personLL,
      class_proportions = class_proportions
    )

    LL_change <- observed_data_LL - observed_data_LL0

    # check stable clustering and convergence and break loop if applicable:
    if (stable_clustering | LL_change < convergence_criterion) {
      converged <- TRUE
      break
    }

    # check if maximum number of iterations has been reached:
    if (it == maxit) {
      converged <- FALSE
    }

    observed_data_LL0 <- observed_data_LL
  }

  if (verbose) {
    message(glue::glue("Random start for seed {seed} finished."))
  }

  # save relevant objects in list ----
  output <- list(
    "observed_data_LL" = observed_data_LL,
    "post" = post,
    "class_proportions" = class_proportions,
    "personLL" = personLL,
    "clustermodels_run" = clustermodels_run,
    "seed" = seed,
    "converged" = converged
  )

  return(output)
}

#### EStep (internal helper) ####
# perform the E-Step of the EM algorithm
EStep <- function(pi_ks, ngroup, nclus, loglik) {
  max_g <- rep(0, ngroup)
  z_gks <- matrix(NA, nrow = ngroup, ncol = nclus)

  for (g in 1:ngroup) {
    for (k in 1:nclus) {
      z_gks[g, k] <- log(pi_ks[k]) + loglik[g, k]
    }
    max_g[g] <- max(z_gks[g, ]) # prevent arithmetic underflow
    z_gks[g, ] <- exp(z_gks[g, ] - rep(max_g[g], nclus))
  }

  # divide by the rowwise sum of the above calculated part
  z_gks <- diag(1 / apply(z_gks, 1, sum)) %*% z_gks

  return(z_gks)
}
# taken from https://github.com/AndresFPA/mmgsem/blob/main/R/E_Step.R

#### generate_startvalues (internal helper) ####
# function to check whether random values are in the bounds specified in the model matrices
within_bounds <- function(values, model_matrix) {
  free <- model_matrix$free
  x  <- values[free]
  lower <- model_matrix$lbound[free]
  upper <- model_matrix$ubound[free]

  all(is.na(lower) | x >= lower) &&
    all(is.na(upper) | x <= upper)
}

# generate random starting values for OpenMx
generate_startvalues <- function(n_clusters, n_factors, personmodel_list) {
  startvalues <- vector(mode = "list", length = n_clusters)
  A <- personmodel_list[[1]]$matrices[[1]]
  free_phi <- A$free
  n_free_phi <- sum(free_phi)
  labels_phi <- A$labels[free_phi]

  Q <- personmodel_list[[1]]$matrices[[5]]
  free_zeta <- Q$free
  n_free_zeta <- sum(free_zeta)
  # select only labels on the lower triangle of the matrix (including diagonal) to avoid duplicates:
  labels_zeta <- Q$labels[lower.tri(Q$labels, diag = TRUE) & free_zeta]

  P0 <- personmodel_list[[1]]$matrices[[8]]
  free_P0 <- P0$free
  n_free_P0 <- sum(free_P0)
  # select only labels on the lower triangle of the matrix (including diagonal) to avoid duplicates:
  labels_P0 <- P0$labels[lower.tri(P0$labels, diag = TRUE) & free_P0]

  for (k in 1:n_clusters) {
    # attempt to create random starting values that are within the bounds specified in the model matrices
    maxtries <- 1000
    for (attempt in 1:maxtries) {
      # generate a stationary matrix of regression coefficients
      # phi parameters ----
      # create multiplicative noise (.5 to 1.5) for free phi parameters:
      noise_multi <- matrix(1, nrow = nrow(A), ncol = ncol(A))
      noise_multi[free_phi] <- stats::runif(n_free_phi, 0.5, 1.5)
      # create additive noise (-0.5 to 0.5) for free phi parameters:
      noise_add <- matrix(0, nrow = nrow(A), ncol = ncol(A))
      noise_add[free_phi] <- stats::runif(n_free_phi, -0.5, 0.5)
      # add noise to starting values of phi:
      values_phi <- A$values * noise_multi + noise_add

      # check if within the bounds, then fix to bounds if not:
      below_bounds <- free_phi & !is.na(A$lbound) & values_phi < A$lbound
      values_phi[below_bounds] <- A$lbound[below_bounds]
      above_bounds <- free_phi & !is.na(A$ubound) & values_phi > A$ubound
      values_phi[above_bounds] <- A$ubound[above_bounds]

      # check if the largest eigenvalue is greater than .9 in modulus:
      # (select 1:n_factors to account for possible random intercept specification)
      ev <- eigen(values_phi[1:n_factors, 1:n_factors])$values
      max_modulus <- max(Mod(ev))
      if (max_modulus > .9) {
        # if yes, rescale the matrix:
        values_phi[1:n_factors, 1:n_factors] <- values_phi[1:n_factors, 1:n_factors] * (.9 / max_modulus)
      }

      # zeta parameters ----
      # generate a positive definitive matrix of innovation (co)variances
      # create multiplicative noise (.5 to 1.5) for free zeta parameters:
      noise_multi <- matrix(1, nrow = nrow(Q), ncol = ncol(Q))
      noise_multi[free_zeta] <- stats::runif(n_free_zeta, 0.5, 1.5)
      # create additive noise (-0.5 to 0.5) for free zeta parameters:
      noise_add <- matrix(0, nrow = nrow(Q), ncol = ncol(Q))
      noise_add[free_zeta] <- stats::runif(n_free_zeta, -0.5, 0.5)
      # add noise to starting values of zeta:
      values_zeta <- Q$values * noise_multi + noise_add

      # check if within the bounds, then fix to bounds if not:
      below_bounds <- free_zeta & !is.na(Q$lbound) & values_zeta < Q$lbound
      values_zeta[below_bounds] <- Q$lbound[below_bounds]
      above_bounds <- free_zeta & !is.na(Q$ubound) & values_zeta > Q$ubound
      values_zeta[above_bounds] <- Q$ubound[above_bounds]

      # turn the matrix into a positive definitive matrix:
      # (select 1:n_factors to account for possible random intercept specification)
      values_zeta[1:n_factors, 1:n_factors] <- Matrix::nearPD(values_zeta[1:n_factors, 1:n_factors])$mat |> as.matrix()

      # P0 parameters ----
      # generate a positive definitive matrix of initial (co)variances (P0)
      # create multiplicative noise (.5 to 1.5) for free P0 parameters:
      noise_multi <- matrix(1, nrow = nrow(P0), ncol = ncol(P0))
      noise_multi[free_P0] <- stats::runif(n_free_P0, 0.5, 1.5)
      # create additive noise (-0.5 to 0.5) for free P0 parameters:
      noise_add <- matrix(0, nrow = nrow(P0), ncol = ncol(P0))
      noise_add[free_P0] <- stats::runif(n_free_P0, -0.5, 0.5)
      # add noise to starting values of P0:
      values_P0 <- P0$values * noise_multi + noise_add

      # check if within the bounds, then fix to bounds if not:
      below_bounds <- free_P0 & !is.na(P0$lbound) & values_P0 < P0$lbound
      values_P0[below_bounds] <- P0$lbound[below_bounds]
      above_bounds <- free_P0 & !is.na(P0$ubound) & values_P0 > P0$ubound
      values_P0[above_bounds] <- P0$ubound[above_bounds]

      # turn the matrix into a positive definitive matrix:
      # (select 1:n_factors to account for possible random intercept specification)
      values_P0[1:n_factors, 1:n_factors] <- Matrix::nearPD(values_P0[1:n_factors, 1:n_factors])$mat |> as.matrix()
      # if there's a random intercept, also do this for the lower right sub-matrix
      if (nrow(values_P0) > n_factors && ncol(values_P0) > n_factors) {
        values_P0[(n_factors + 1):nrow(values_P0), (n_factors + 1):ncol(values_P0)] <- Matrix::nearPD(values_P0[(n_factors + 1):nrow(values_P0), (n_factors + 1):ncol(values_P0)])$mat |> as.matrix()
      }

      # check again if the generated starting values are within the bounds specified in the model matrices
      # (because the rescaling and nearPD function might have changed the values)
      if (within_bounds(values_phi, A) &&
          within_bounds(values_zeta, Q) &&
          within_bounds(values_P0, P0)) {

        ## create named vector of starting values:
        # extract the new starting values assign the corresponding labels as names:
        startvalues_phi <- values_phi[free_phi]
        names(startvalues_phi) <- labels_phi

        # only obtain lower triangle to avoid duplicates in the covariances:
        startvalues_zeta <- values_zeta[lower.tri(values_zeta, diag = TRUE) & free_zeta]
        names(startvalues_zeta) <- labels_zeta

        startvalues_P0 <- values_P0[lower.tri(values_P0, diag = TRUE) & free_P0]
        names(startvalues_P0) <- labels_P0

        break
      }

      # if the maximum number of attempts is reached and the generated starting values are still not within bounds, use the pre-specified starting values from the model matrices
      if (attempt == maxtries) {
        # create named vector of starting values, using the original starting values without random perturbation where needed:
        if (!within_bounds(values_phi, A)) {
          startvalues_phi <- A$values[free_phi]
        } else {
          startvalues_phi <- values_phi[free_phi]
        }
        names(startvalues_phi) <- labels_phi
        if (!within_bounds(values_zeta, Q)) {
          startvalues_zeta <- Q$values[lower.tri(Q$values, diag = TRUE) & free_zeta]
        } else {
          startvalues_zeta <- values_zeta[lower.tri(values_zeta, diag = TRUE) & free_zeta]
        }
        names(startvalues_zeta) <- labels_zeta

        if (!within_bounds(values_P0, P0)) {
          startvalues_P0 <- P0$values[lower.tri(P0$values, diag = TRUE) & free_P0]
        } else {
          startvalues_P0 <- values_P0[lower.tri(values_P0, diag = TRUE) & free_P0]
        }
        names(startvalues_P0) <- labels_P0
      }
    }

    startvalues[[k]] <- c(
      startvalues_phi,
      startvalues_zeta,
      startvalues_P0
    )

  }
  return(startvalues)
}

#### compute_observed_data_LL (internal helper) ####
# compute observed LL from personwise LL and class proportions
compute_observed_data_LL <- function(personLL, class_proportions) {
  # sum with the log of class proportions:
  personLL_weighted <- sweep(personLL, 2, log(class_proportions), "+")
  # get the max value per row:
  max_i <- apply(personLL_weighted, 1, max)
  # subtract max value from each row (prevent arithmetic underflow):
  minus_max <- sweep(personLL_weighted, 1, max_i, "-")
  # exp to get out of log space:
  personLL_exp <- exp(minus_max)
  # sum per row and then take the log again:
  personLL_summed <- log(rowSums(personLL_exp))
  # re-add the max value and then sum to obtain observed data LL
  observed_data_LL <- sum(personLL_summed + max_i)

  return(observed_data_LL)
}
