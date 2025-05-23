#' Step 3: Estimate the Structural Model
#'
#' @description
#' This function performs Step 3 of 3S-LVAR: Estimate the Structural Model using the factor scores as single indicators while accounting for their inherent uncertainty.
#' The estimation is implemented in the State Space Model framework in `OpenMx`.
#'
#' @param step2output The output obtained with the `step2()` function.
#' @param id A string containing the name of the ID variable.
#' @param A An `OpenMx` matrix object that describes the regression coefficients in the State Space model. Diagonal entries represent autoregressive effects, and off-diagonal entries represent cross-lagged effects. See [OpenMx::mxExpectationStateSpace]. The helper function [create_A()] is available to create this object.
#' @param Q An `OpenMx` matrix object that describes the innovation (co)variance matrix. See [OpenMx::mxExpectationStateSpace]. The helper function [create_Q()] is available to create this object.
#' @param B An `OpenMx` matrix object that describes the covariate effects on the latent variables. Optional. If NULL, a zero matrix (i.e., no covariate effects) is automatically created.
#' @param D An `OpenMx` matrix object that describes the covariate effects on the observed items. Optional. If NULL, a zero matrix (i.e., no covariate effects) is automatically created.
#' @param x0 An `OpenMx` matrix object that describes the initial latent variable vector that initiates the Kalman Filter. Optional. If NULL, an object is automatically created where all values are fixed to 0.
#' @param P0 An `OpenMx` matrix object that describes the initial error covariance matrix that initiates the Kalman Filter. Optional. If NULL, a matrix will be automatically created. The parameters are freely estimated, with large starting values (100 for the diagonal and 10 for the off-diagonal).
#' @param u An `OpenMx` matrix object that indicates the covariates. Required if `B` or `D` are included. If NULL, an object with no covariates is automatically created.
#' @param step3group A string containing the name of the grouping variable for the structural model.
#' @param newdata A data frame. Note that the output of `step2()` already contains the data, so using this argument is only required if the data have been manipulated after step 2 (e.g., removing outliers).
#' @param tryhard Should the algorithm run multiple times to obtain a solution? (See [OpenMx::mxTryHard])
#'
#' @return A list containing the following elements:
#'
#' `data` The data set which served as input.
#'
#' `estimates` A vector that includes the parameter estimates.
#'
#' `model` The fitted model (an OpenMx object).
#'
#' @export
#' @importFrom OpenMx imxReportProgress


step3 <- function(step2output, id, A, Q,
                  B = NULL, D = NULL,
                  x0 = NULL, P0 = NULL,
                  u = NULL,
                  step3group = NULL,
                  newdata = NULL,
                  tryhard = FALSE){

  #### 1) Preparations ####
  ## extract objects from step 1 output:
  if(is.null(newdata)){
    data <- step2output$data |> as.data.frame()
  } else {
    data <- newdata
  }

  data <- data |> as.data.frame()
  # there are known issues when the supplied data are a tibble, so we transform it into a plain data.frame

  if(!is.character(data[[id]])){
    warning("Your id variable has been transformed into a character.")
    data[id] <- data[[id]] |> as.character()
  }

  lambda_star <- step2output$lambda_star
  theta_star <- step2output$theta_star
  factors <- step2output$other$factors
  factors_ind <- paste0(factors, "_ind")                                        # names of factor score variables (single indicators)
  n_factors <- length(factors)
  step1group <- step2output$other$step1group                                    # name of the grouping variable in step1
  unique_ids <- data[[id]] |> unique()                                          # vector of unique ids
  N <- length(unique_ids)

  #### 2) data manipulation ####
  ## rename the factor score variables in the data
  # to use them as indicators of the latent variables
  data <- data |>
    dplyr::rename_with(~ factors_ind, all_of(factors))

  #### 3) create OpenMx matrices ####
  xdim <- nrow(A$free)                                                          # number of latent variables (depends whether there is a random intercept or not)
  udim <- 1                                                                     # exogenous covariates (ignored so far)
  ydim <- n_factors                                                             # number of indicators (i.e., factor score variables)

  # B matrix (= exogenous covariates on latent constructs)
  if(is.null(B)){
    B <- OpenMx::mxMatrix("Zero", name = "B",
                          nrow = xdim, ncol = udim)
  }

  # D matrix (= exogenous covariates on observed variables)
  if(is.null(D)){
    D <- OpenMx::mxMatrix("Zero", name = "D",
                          nrow = ydim, ncol = udim)
  }


  # x0 and P0 (= initial values and (co)variances of the latent constructs)
  if(is.null(x0)){
    if(xdim > ydim){
      labels = c(paste0("ini_", factors),
                 paste0("m_intercept_", factors))
    } else {
      labels = paste0("ini_", factors)
    }

    x0 <- OpenMx::mxMatrix("Full", name = "x0",
                           nrow = xdim, ncol = 1,
                           free = FALSE,
                           values = 0,
                           labels = labels)
  }

  if(is.null(P0)){
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
    if(xdim > ydim){
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
  if(is.null(u)){
    u <- OpenMx::mxMatrix("Zero", name = "u", nrow = udim, ncol = 1)
  }

  #### 4) create OpenMx models ####
  ## if there is no grouping variable in step3:
  if(purrr::is_empty(step3group)){
    # create a list of models (one for each individual) for each latent class:
    personmodelnames <- paste0("id_", unique_ids)
    names(personmodelnames) <- unique_ids

    personmodel_list <- vector(mode = "list", length = N)
    names(personmodel_list) <- unique_ids
    for(i in unique_ids){
      # which step1 group (if any) does the individual belong to?
      if(!purrr::is_empty(step1group)){
        s1g <- data[data[[id]] == i, step1group] |> unique()
      } else {
        s1g <- 1
      }

      # fix the loadings to lambda_star (potentially depending on step1 group)
      # C matrix
      values <- matrix(0, nrow = n_factors, ncol = n_factors)
      diag(values) <- lambda_star[s1g, ]
      C_names <- list(factors_ind, factors)
      if(xdim > ydim){
        values = cbind(values, values)
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

      # fix the residual variances to theta_star (potentially depending on step 1 group)
      # R matrix (= measurement noise)
      R <- OpenMx::mxMatrix("Diag", name = "R",
                            nrow = ydim, ncol = ydim,
                            free = FALSE,
                            values = theta_star[s1g, ],
                            labels = NA
      )

      # create the model
      modelname <- personmodelnames[i]  |> as.character()
      dat_person <- data[data[[id]] == i, factors_ind]
      personmodel_list[[i]] <- OpenMx::mxModel(name = modelname,
                                               A, B, C, D,
                                               Q, R, x0, P0,
                                               u,
                                               OpenMx::mxExpectationStateSpace('A', 'B', 'C', 'D',
                                                                               'Q', 'R', 'x0', 'P0',
                                                                               'u'),
                                               OpenMx::mxFitFunctionML(),
                                               OpenMx::mxData(dat_person,
                                                              'raw'))
    }

    # combine the person-models to a multi-subject model
    fullmodel <- OpenMx::mxModel("fullmodel",
                                 personmodel_list,
                                 OpenMx::mxFitFunctionMultigroup(personmodelnames))
    # fit the model
    if(tryhard){
      fullmodelr <- OpenMx::mxTryHard(fullmodel)
    } else {
      fullmodelr <- OpenMx::mxRun(fullmodel)
    }
  }




  ## if there is a grouping variable in step 3
  if(!purrr::is_empty(step3group)){
    if(!is.character(data[[step3group]])){
      warning("Your grouping variable has been transformed into a character.")
      data[step3group] <- data[[step3group]] |> as.character()
    }

    # create a list of models (one for each individual) for each latent class:
    personmodelnames <- paste0("id_", unique_ids)
    names(personmodelnames) <- unique_ids

    personmodel_list <- vector(mode = "list", length = N)
    names(personmodel_list) <- unique_ids
    # create the person-models:
    for(g in unique(data[[step3group]])){
      group_ids <- data[data[[step3group]] == g, id] |> unique()
      for(i in group_ids){
        if(!purrr::is_empty(step1group)){
          s1g <- data[data[[id]] == i, step1group] |> unique()
        } else {
          s1g <- 1
        }

        # fix the loadings to lambda_star (potentially depending on step1 group)
        # C matrix
        values <- matrix(0, nrow = n_factors, ncol = n_factors)
        diag(values) <- lambda_star[s1g, ]
        C_names <- list(factors_ind, factors)
        if(xdim > ydim){
          values = cbind(values, values)
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

        # fix the residual variances to theta_star (potentially depending on step 1 group)
        # R matrix (= measurement noise)
        R <- OpenMx::mxMatrix("Diag", name = "R",
                              nrow = ydim, ncol = ydim,
                              free = FALSE,
                              values = theta_star[s1g, ],
                              labels = NA
        )

        modelname <- personmodelnames[i]  |> as.character()
        dat_person <- data[data[[id]] == i, factors_ind]
        temp_model <- OpenMx::mxModel(name = modelname,
                                      A, B, C, D,
                                      Q, R, x0, P0,
                                      u,
                                      OpenMx::mxExpectationStateSpace('A', 'B', 'C', 'D',
                                                                      'Q', 'R', 'x0', 'P0',
                                                                      'u'),
                                      OpenMx::mxFitFunctionML(),
                                      OpenMx::mxData(dat_person,
                                                     'raw'))
        temp_model <- OpenMx::omxSetParameters(temp_model,
                                               labels = names(coef(temp_model)),
                                               newlabels = paste0(names(coef(temp_model)),
                                                                  "_", g))
        personmodel_list[[i]] <- temp_model
      }
    }

    fullmodel <- OpenMx::mxModel("fullmodel", personmodel_list,
                                 OpenMx::mxFitFunctionMultigroup(personmodelnames))
    if(tryhard){
      fullmodelr <- OpenMx::mxTryHard(fullmodel)
    } else {
      fullmodelr <- OpenMx::mxRun(fullmodel)
    }

  }




  #### 7) build the output ####
  estimates <- coef(fullmodelr)

  output <- list("data" = data,
                 "estimates" = estimates,
                 "model" = fullmodelr)


  return(output)
}
