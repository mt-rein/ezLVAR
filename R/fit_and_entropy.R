#### functions for fit and entropy ####
#### r2_entropy ####

#' Compute the R-squared entropy of a mixture Step 3 model.
#'
#' @description Computes the R-squared entropy index for a mixture model, which
#'   quantifies the proportional reduction in classification uncertainty gained
#'   from using posterior class probabilities instead of prior class
#'   proportions. A value close to 1 indicates high classification certainty,
#'   whereas a value close to 0 indicates substantial uncertainty in class
#'   assignments.
#'
#' @references J.K. Vermunt and J. Magidson (2016). Technical Guide for Latent
#'   GOLD 5.1: Basic, Advanced, and Syntax. Belmont, MA: Statistical Innovations
#'   Inc.
#'
#' @param step3output An object obtained with the [step3] function that includes
#'   a mixture model (i.e., `mixture = TRUE`).
#'
#' @returns A numeric value.
#'
#' @export
r2_entropy <- function(step3output) {
  if (!inherits(step3output, "3slvar_step3")) {
    stop("The input must be the output of the step3() function.")
  }

  if (step3output[["type"]] != "mixture") {
    stop(
      "The input must be the output of the step3() function including a mixture model."
    )
  }

  prior <- step3output$class_proportions
  posterior <- step3output$posterior_probabilities[,
    2:(1 + step3output$n_groups)
  ] |>
    as.matrix()
  # avoid log(0): replace 0 with a tiny number
  posterior[posterior < 1e-12] <- 1e-12

  entropy_prior <- sum(-prior * log(prior))
  entropy_posterior <- (-posterior * log(posterior)) |>
    rowSums() |>
    mean()
  r2 <- (entropy_prior - entropy_posterior) / entropy_prior
  return(r2)
}

#### fit_indices ####

#' Compute fit measures of Step 3 model
#'
#' @description Computes likelihood based fit indices for a Step 3 model:
#'   log-likelihood, AIC, AIC3, BIC, and ICL (only for mixture models).
#'
#' @references J.K. Vermunt and J. Magidson (2016). Technical Guide for Latent
#'   GOLD 5.1: Basic, Advanced, and Syntax. Belmont, MA: Statistical Innovations
#'   Inc.
#'
#' @param step3output An object obtained with the [step3] function that includes
#'   a mixture model (i.e., `mixture = TRUE`).
#'
#' @returns A named vector of numeric values.
#'
#' @export
fit_indices <- function(step3output) {
  if (!inherits(step3output, "3slvar_step3")) {
    stop("The input must be the output of the step3() function.")
  }

  n_parameters <- step3output$n_parameters
  n_persons <- step3output$n_persons
  logLik <- step3output$logLik

  AIC <- 2 * n_parameters - 2 * logLik
  AIC3 <- 3 * n_parameters - 2 * logLik
  BIC <- n_parameters * log(n_persons) - 2 * logLik

  if (step3output[["type"]] == "mixture") {
    # posterior (needed for entropy)
    posterior <- step3output$posterior_probabilities[,
      2:(1 + step3output$n_groups)
    ] |>
      as.matrix()
    # avoid log(0): replace 0 with a tiny number
    posterior[posterior < 1e-12] <- 1e-12

    # compute total entropy
    entropy <- (-posterior * log(posterior)) |> sum()
    ICL <- BIC + 2 * entropy
  } else {
    # when there is no mixture modeling, the entropy is 0
    ICL <- NULL
  }

  return(c(
    "logLik" = logLik,
    "AIC" = AIC,
    "AIC3" = AIC3,
    "BIC" = BIC,
    "ICL" = ICL
  ))
}
