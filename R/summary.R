#### custom summary functions ####
#### Step 1 ####
#' Summary of Step 1
#'
#' @param step1output The output of the `step1` function, an object of class `3slvar_step1`
#' @param print_MM_summary Should the summary of the estimated measurement models be printed (i.e, call `lavaan::summary` on all blocks?
#' @param ... Arguments that are forwarded to `lavaan::summary()`.
#'
#' @export
summary.3slvar_step1 <- function(step1output, print_MM_summary = TRUE, ...) {
  fit_step1 <- step1output$mm_output
  n_blocks <- length(fit_step1)
  group_var <- lavaan::lavInspect(fit_step1[[1]], "group") |> as.character()
  group_names <- lavaan::lavInspect(fit_step1[[1]], "group.label")
  if (rlang::is_empty(group_var)) {
    n_groups <- 1
  } else {
    n_groups <- length(group_names)
  }

  cat("Number of measurement blocks: ", n_blocks, "\n")
  if (!rlang::is_empty(group_var)) {
    cat("Grouping variable (number of groups): ", group_var, "(", n_groups, ")\n")
  }

  if(print_MM_summary) {
    for (block in 1:n_blocks) {
      cat("Summary of measurement block ", block, ":\n")
      print(lavaan::summary(fit_step1[[block]], ...))
      cat("\n\n")
    }
  }

  invisible(step1output)
}

#### Step 2 ####
#' Summary of Step 2
#'
#' @param step2output The output of the `step2` function, an object of class `3slvar_step2`
#'
#' @export
summary.3slvar_step2 <- function(step2output) {
  factors <- step2output$other$factors

  descriptives_lambdastar <- data.frame(factor = character(),
                                        average = numeric(),
                                        min = numeric(),
                                        max = numeric())
  for (i in 1:length(factors)) {
    fac <- factors[i]
    average <- mean(step2output$lambda_star[[fac]]) |> round(2)
    minimum <- min(step2output$lambda_star[[fac]]) |> round(2)
    maximum <- max(step2output$lambda_star[[fac]]) |> round(2)
    descriptives_lambdastar[i, ] <- c(fac, average, minimum, maximum)
  }

  descriptives_thetastar <- data.frame(factor = character(),
                                        average = numeric(),
                                        min = numeric(),
                                        max = numeric())
  for (i in 1:length(factors)) {
    fac <- factors[i]
    average <- mean(step2output$theta_star[[fac]]) |> round(2)
    minimum <- min(step2output$theta_star[[fac]]) |> round(2)
    maximum <- max(step2output$theta_star[[fac]]) |> round(2)
    descriptives_thetastar[i, ] <- c(fac, average, minimum, maximum)
  }

  cat("Descriptives for Lambda* (model-based reliability)\n")
  print(descriptives_lambdastar)
  cat("\n")
  cat("Descriptives for Theta* (residual variance in the factor scores):\n")
  print(descriptives_thetastar)

  invisible(step2output)
}
