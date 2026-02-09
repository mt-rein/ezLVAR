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

  #### EM Loop ####

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

    #### M-step: fitting SSM model and update parameter estimates ####
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
      purrr::map(OpenMx::mxRun, silent = TRUE, suppressWarnings = TRUE)

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

  # save relevant objects in list:
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
# generate random starting values for OpenMx
generate_startvalues <- function(n_clusters, n_factors, personmodel_list) {
  startvalues <- vector(mode = "list", length = n_clusters)
  free_phi <- personmodel_list[[1]]$A$free
  labels_phi <- personmodel_list[[1]]$A$labels[free_phi]
  values_phi <- personmodel_list[[1]]$A$values[free_phi]

  free_zeta <- personmodel_list[[1]]$Q$free
  labels_zeta <- personmodel_list[[1]]$Q$labels[free_zeta]
  values_zeta <- personmodel_list[[1]]$Q$values[free_zeta]

  for (k in 1:n_clusters) {
    # generate a stationary matrix of regression coefficients
    # create additive noise (-.5 to .5) for phi:
    noise_multi <- stats::runif(length(values_phi), 0.5, 1.5)
    phistart <- matrix(values_phi * noise_multi, nrow = n_factors)
    # check if the largest eigenvalue is greater than .9 in modulus:
    ev <- eigen(phistart)$values
    max_modulus <- max(Mod(ev))
    if (max_modulus > .9) {
      # if yes, rescale the matrix:
      phistart_scaled <- phistart * (.9 / max_modulus)
    } else {
      phistart_scaled <- phistart
    }

    # generate a positive definitive matrix of innovation (co)variances
    noise_multi <- stats::runif(length(values_phi), 0.5, 1.5)
    zetastart <- matrix(values_zeta * noise_multi, nrow = n_factors)
    zetastartPD <- Matrix::nearPD(zetastart)$mat |> as.matrix()

    startvalues[[k]] <- c(
      c(phistart_scaled),
      zetastartPD[lower.tri(zetastartPD, diag = TRUE)]
    )
    names(startvalues[[k]]) <- c(labels_phi, unique(labels_zeta))
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
