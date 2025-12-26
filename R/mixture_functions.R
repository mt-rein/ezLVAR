mixture_phase1 <- function(personmodel_list,
                           n_clusters,
                           n_factors,
                           n_persons,
                           n_starts,
                           stablecluster_criterion,
                           GEM_iterations,
                           maxit,
                           seeds,
                           verbose) {

  # provide seeds for the multiple starts (for replicability)
  if (is.null(seeds)) {
    all_seeds <- sample(1:100000000, n_starts)
  }

  # create list to store the results per random start:
  all_starts <- vector(mode = "list", length = n_starts)
  # create list with names of objective functions (needed in M-Step)
  objectives <- sapply(personmodel_list, function(x) x$name) |>
    as.character() |>
    paste0(".objective")

  #### Random Start Loop ####
  for (random_start in 1:n_starts) {
    if (verbose) {
      message(glue::glue("==== Start {random_start} ===="))
    }

    # set seed:
    set.seed(all_seeds[random_start])
    # create random hard cluster assignment:
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

    # reset the observed data LL:
    observed_data_LL0 <- -Inf
    # reset the check for stable clustering:
    stable_clustering <- FALSE

    #### EM Loop ####
    # loop over iterations until stable clustering is achieved or max iterations are reached
    for (it in 1:maxit) {
      #### E-step: update class membership and class proportions ####
      if (it > 1) {
        post0 <- post
        # update posteriors:
        post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                      nclus = n_clusters, loglik = personLL)
        # compute (absolute) differences in posteriors
        delta <- abs(post - post0)
        person_delta <- rowSums(delta)
        if (verbose) {
          message(glue::glue(
            "Highest person-delta in posterior probabilities: {round(max(person_delta, na.rm = TRUE), 3)}. ",
            "Average person-delta: {round(mean(person_delta, na.rm = TRUE), 3)}."
          ))
        }

        # check if stable clustering has been achieved:
        if (mean(person_delta) < stablecluster_criterion) {
          stable_clustering <- TRUE
          if (verbose) {
            message(glue::glue(">> Start {random_start} has a stable clustering. Switched to full M-Step."))
          }
        }
      }

      # compute class proportions:
      class_proportions <- colMeans(post)

      #### M-step: fitting SSM model and update parameter estimates ####
      # get start values from the coefficients of previous iteration:
      if (it > 1) {
        startvalues <- clustermodels_run |>
          purrr::map(coef)
      } else {
        # in the first iteration, generate random start values in each cluster
        startvalues <- generate_startvalues(n_clusters = n_clusters,
                                            n_factors = n_factors,
                                            personmodel_list = personmodel_list)
      }


      # create one multi-group (i.e., multi-subject) model per cluster
      # person-models are weighted by their posterior probabilities
      clustermodels <- vector(mode = "list", length = n_clusters)
      for (k in 1:n_clusters) {
        clustername <- paste0("model_k", k)
        weighted_objectives <- paste(post[, k], "*", objectives,
                                     collapse = " + ")
        model_k <- OpenMx::mxModel(clustername, personmodel_list,
                                   OpenMx::mxAlgebraFromString(weighted_objectives,
                                                               name = "weightedfit"),
                                   OpenMx::mxFitFunctionAlgebra("weightedfit"))
        # add start values:
        model_k <- OpenMx::omxSetParameters(model_k,
                                            labels = names(startvalues[[k]]),
                                            values = startvalues[[k]])

        # change options (skip computation of hessian and SEs)
        model_k <- OpenMx::mxOption(model_k,
                                    key = "Calculate Hessian",
                                    value = "No")
        model_k <- OpenMx::mxOption(model_k,
                                    key = "Standard Errors",
                                    value = "No")
        # if stable clustering has not yet been achieved, reduce iterations
        if (!stable_clustering) {
          model_k <- OpenMx::mxOption(model_k,
                                      key = "Major iterations",
                                      value = GEM_iterations)
        }

        clustermodels[[k]] <- model_k
      }
      names(clustermodels) <- paste0("model_k", 1:n_clusters)

      # run the models
      clustermodels_run <- clustermodels |>
        purrr::map(OpenMx::mxRun,
                   silent = TRUE,
                   suppressWarnings = TRUE)

      # obtain person-wise LL in a n_persons x n_clusters matrix:
      personLL <- clustermodels_run |>
        purrr::map(~ purrr::map_dbl(.x$submodels, ~ .x$fitfunction$result)) |>
        simplify2array()
      personLL <- personLL / (-2)
      # openMx gives the minus2LL, so we divide by minus 2

      # compute observed-data log likelihood from person-wise LL:
      observed_data_LL <- compute_observed_data_LL(personLL = personLL,
                                                   class_proportions = class_proportions)

      LL_change <- observed_data_LL - observed_data_LL0

      if (verbose) {
        message(glue::glue(
          "Iteration {it}: ",
          "Log Likelihood: {round(observed_data_LL, 3)}; ",
          "Change: {round(LL_change, 6)}."
        ))
      }

      # check stable clustering and break loop if applicable:
      if (stable_clustering) {
        break
      }

      observed_data_LL0 <- observed_data_LL

      # check if maxit has been reached:
      if (verbose && it == maxit) {
        message(glue::glue(
          "Start {random_start} did not arrive at a stable clustering."
        ))
      }

    }

    # save relevant objects in list:
    all_starts[[random_start]] <- list("observed_data_LL" = observed_data_LL,
                                       "class_proportions" = class_proportions,
                                       "personLL" = personLL,
                                       "clustermodels_run" = clustermodels_run,
                                       "seed" = all_seeds[random_start])
  }

  return(all_starts)
}

mixture_phase2 <- function(output_phase1,
                           personmodel_list,
                           n_clusters,
                           n_factors,
                           n_persons,
                           n_best_starts,
                           convergence_criterion,
                           maxit,
                           verbose) {

  all_starts <- output_phase1

  best_starts <- all_starts |>
    purrr::map_dbl(~ .x$observed_data_LL) |>
    order(decreasing = TRUE) |>
    head(n_best_starts)


  # create list with name of objective functions (needed in M-Step)
  objectives <- sapply(personmodel_list, function(x) x$name) |>
    as.character() |>
    paste0(".objective")

  best_loglik <- -Inf
  nonconvergences <- 0

  #### Best start loops ####
  for (random_start in 1:n_best_starts) {
    if (verbose) {
      message(glue::glue("==== Start {random_start} ===="))
    }
    start_number <- best_starts[random_start]

    # extract outputs from corresponding start in Phase 1:
    observed_data_LL0 <- all_starts[[start_number]]$observed_data_LL
    class_proportions <- all_starts[[start_number]]$class_proportions
    personLL <- all_starts[[start_number]]$personLL
    clustermodels_run <- all_starts[[start_number]]$clustermodels_run
    start_seed <- all_starts[[start_number]]$seed
    # (i.e., pick up the algorithm from Phase 1)
    set.seed(start_seed)

    #### EM Loop ####
    # loop over iterations until convergence or max iterations are reached
    for (it in 1:maxit) {
      #### E-step: update class membership and class proportions ####
      post <- EStep(pi_ks = class_proportions, ngroup = n_persons,
                    nclus = n_clusters, loglik = personLL)

      # compute class proportions:
      class_proportions <- colMeans(post)

      #### M-step: fitting SSM model and update parameter estimates ####
      # get start values from the coefficients of previous iteration:
      startvalues <- clustermodels_run |>
        purrr::map(coef)

      # create one multi-group (i.e., multi-subject) model per cluster
      # person-models are weighted by their posterior probabilities
      clustermodels <- vector(mode = "list", length = n_clusters)
      for (k in 1:n_clusters) {
        clustername <- paste0("model_k", k)
        weighted_objectives <- paste(post[, k], "*", objectives,
                                     collapse = " + ")
        model_k <- OpenMx::mxModel(clustername, personmodel_list,
                                   OpenMx::mxAlgebraFromString(weighted_objectives,
                                                               name = "weightedfit"),
                                   OpenMx::mxFitFunctionAlgebra("weightedfit"))
        # add start values:
        model_k <- OpenMx::omxSetParameters(model_k,
                                            labels = names(startvalues[[k]]),
                                            values = startvalues[[k]])

        # change options (skip computation of hessian and SEs)
        model_k <- OpenMx::mxOption(model_k,
                                    key = "Calculate Hessian",
                                    value = "No")
        model_k <- OpenMx::mxOption(model_k,
                                    key = "Standard Errors",
                                    value = "No")

        clustermodels[[k]] <- model_k
      }
      names(clustermodels) <- paste0("model_k", 1:n_clusters)

      # run the models
      clustermodels_run <- clustermodels |>
        purrr::map(OpenMx::mxRun,
                   silent = TRUE,
                   suppressWarnings = TRUE)

      # obtain person-wise LL in a n_persons x n_clusters matrix:
      personLL <- clustermodels_run |>
        purrr::map(~ purrr::map_dbl(.x$submodels, ~ .x$fitfunction$result)) |>
        simplify2array()
      personLL <- personLL / (-2)

      # compute observed-data log likelihood from person-wise LL:
      observed_data_LL <- compute_observed_data_LL(personLL = personLL,
                                                   class_proportions = class_proportions)

      LL_change <- observed_data_LL - observed_data_LL0

      if (verbose) {
        message(glue::glue(
          "Iteration {it}: ",
          "Log Likelihood: {round(observed_data_LL, 3)}; ",
          "Change: {round(LL_change, 6)}."
        ))
      }

      # check convergence and break loop if applicable:
      if (LL_change < convergence_criterion) {
        if (verbose) {
          message("Convergence achieved.")
        }
        break
      }

      observed_data_LL0 <- observed_data_LL

      # check if maximum number of iterations has been reached:
      if (it == maxit) {
        if (verbose) {
          message(glue::glue(
            "Max iterations reached without convergence. Start: {random_start}"
          ))
        }
        nonconvergences <- nonconvergences + 1
      }
    }

    # check if the new fit is better than the previous best one
    if (observed_data_LL > best_loglik) {
      best_startnumber <- start_number
      best_seed <- start_seed
      best_loglik <- observed_data_LL
      best_post <- post
      best_model <- clustermodels_run
    }
  }

  output <- list("best_startnumber" = best_startnumber,
                 "best_model" = best_model,
                 "best_post" = best_post,
                 "best_loglik" = best_loglik,
                 "best_seed" = best_seed,
                 "nonconvergenced_starts" = nonconvergences)

  return(output)
}
