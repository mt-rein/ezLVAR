#### model selection core function ####

#' Model selection for Step 3 mixture models
#'
#' @description Compares multiple Step 3 mixture models with different numbers
#'   of clusters using various information criteria (log-likelihood, AIC, AIC3,
#'   BIC, ICL) and the convex hull method.
#'
#' @details The [summary] function provides an overview of the model selection
#'   results, including which model is selected by each criterion. The [plot]
#'   function visualizes the information criteria across models, allowing for a
#'   visual comparison of model fit and complexity.
#'
#' @param model_list A list of objects obtained with the [step3] function. Each
#'   object must include a mixture model (i.e., `mixture = TRUE`) with a
#'   different number of clusters. The list can be in any order, but all models
#'   must have been fitted on the same data and with the same A and Q matrices.
#' @param use_CHull A logical value indicating whether to perform convex hull
#'   model selection. Default is TRUE. If FALSE, only the information criteria
#'   will be computed and compared.
#' @param CHull_criterion A character string specifying the criterion to use for
#'   the convex hull model selection. Must be one of: "logLik", "AIC", "AIC3",
#'   "BIC", or "ICL". Default is "logLik".
#'
#' @returns A list of class "model_selection" with the following components:
#'   \item{selection_df}{A data frame with the number of clusters, information
#'   criteria, and scree-test values for each model.}
#'   \item{CHull_criterion}{The criterion used for the convex hull model
#'   selection.}
#'
#' @export
model_selection <- function(
  model_list,
  use_CHull = TRUE,
  CHull_criterion = "logLik"
) {
  #### error when not a list of step3 models
  if (
    !is.list(model_list) ||
      !all(sapply(model_list, function(x) inherits(x, "3slvar_step3")))
  ) {
    stop(
      "model_list must be a list of objects obtained with the step3() function."
    )
  }

  # sort the models by n_groups (in case the user just creates the list in
  # random order)
  idx <- model_list |>
    purrr::map_int(~ .x[["n_groups"]]) |>
    order()
  model_list_sorted <- model_list[idx]

  # create df with n_clusters, LL, n_par, BIC, AIC, AIC3, ICL, R2entropy
  selection_df <- lapply(model_list_sorted, function(x) {
    fitind <- fit_indices(x)
    if (x[["n_groups"]] == 1) {
      fitind["ICL"] = fitind["BIC"]
      rsqu_ent <- NA
    } else {
      rsqu_ent <- r2_entropy(x)
    }
    out <- c(
      "n_clusters" = x[["n_groups"]],
      "logLik" = fitind[["logLik"]],
      "n_par" = x[["n_parameters"]],
      "scree_test" = NA,
      "AIC" = fitind[["AIC"]],
      "AIC3" = fitind[["AIC3"]],
      "BIC" = fitind[["BIC"]],
      "ICL" = fitind[["ICL"]],
      "r2_entropy" = rsqu_ent
    )
  })
  selection_df <- do.call(rbind, selection_df) |> as.data.frame()

  # perform CHull model selection:
  if (use_CHull) {
    output_chull <- CHull(model_list_sorted, criterion = CHull_criterion)
    selection_df$scree_test <- output_chull$df_all$scree_test
  }

  output <- list(
    "selection_df" = selection_df,
    "CHull_criterion" = CHull_criterion
  )
  class(output) <- "model_selection"

  return(output)
}

#### summary function ####

#' @param object An object obtained with the [model_selection] function.
#' @param round_digits Numerical. The number of digits to round numeric values
#'   in the output. Default is 3.

#' @export
#'
#' @rdname model_selection
summary.model_selection <- function(object, round_digits = 3) {
  if (!inherits(object, "model_selection")) {
    stop("The input must be the output of the model_selection() function.")
  }

  cat("Model selection results\n")
  cat(paste0(
    "Comparing models with ",
    paste0(object$selection_df$n_clusters, collapse = ", "),
    " clusters.\n"
  ))
  cat("------------------------------\n")
  # only print CHull if at least one scree-test value is not NA
  if (!all(is.na(object$selection_df$scree_test))) {
    cat(paste0(
      "Convex hull (based on ",
      object$CHull_criterion,
      ") selected the model with ",
      object$selection_df$n_clusters[which.max(object$selection_df$scree_test)],
      " clusters.\n"
    ))
  } else {
    cat(
      "Convex hull model selection could not be performed because no scree-test values were computed.\n"
    )
  }
  cat(paste0(
    "AIC selected the model with ",
    object$selection_df$n_clusters[which.min(object$selection_df$AIC)],
    " clusters.\n"
  ))
  cat(paste0(
    "AIC3 selected the model with ",
    object$selection_df$n_clusters[which.min(object$selection_df$AIC3)],
    " clusters.\n"
  ))
  cat(paste0(
    "BIC selected the model with ",
    object$selection_df$n_clusters[which.min(object$selection_df$BIC)],
    " clusters.\n"
  ))
  cat(paste0(
    "ICL selected the model with ",
    object$selection_df$n_clusters[which.min(object$selection_df$ICL)],
    " clusters.\n"
  ))
  cat("------------------------------\n")
  print_rounded(object$selection_df, digits = round_digits)
}

#### plot function ####

#' @param object An object obtained with the [model_selection] function.
#' @param criterion A character vector specifying which criteria to plot. Must
#'   be one or more of: "AIC", "AIC3", "BIC", "ICL". Default is all four
#'   criteria.
#'
#' @export
#'
#' @rdname model_selection
plot.model_selection <- function(
  object,
  criterion = c("AIC", "AIC3", "BIC", "ICL")
) {
  if (!inherits(object, "model_selection")) {
    stop("The input must be the output of the model_selection() function.")
  }

  df <- object$selection_df

  # validate requested criteria
  if (!all(criterion %in% c("AIC", "AIC3", "BIC", "ICL"))) {
    stop("`criterion` must include one or more of: AIC, AIC3, BIC, ICL.")
  }

  # reshape to long format
  to_keep <- c("n_clusters", criterion)
  long <- reshape(
    df[, to_keep, drop = FALSE],
    varying = criterion,
    v.names = "value",
    timevar = "criterion",
    times = criterion,
    direction = "long"
  )
  rownames(long) <- NULL
  long$criterion <- factor(long$criterion, levels = criterion)

  # prepare colors (one per criterion)
  cols <- seq_along(criterion)

  # empty plot with correct ranges
  plot(
    x = range(long$n_clusters),
    y = range(long$value),
    type = "n",
    xlab = "Number of clusters",
    ylab = "Information criterion",
    xaxt = "n"
  )
  axis(1, at = df$n_clusters, labels = df$n_clusters)

  # draw lines and points for each criterion
  for (i in seq_along(criterion)) {
    crit <- criterion[i]
    sub <- long[long$criterion == crit, ]

    lines(sub$n_clusters, sub$value, col = cols[i], lwd = 2)
    points(sub$n_clusters, sub$value, col = cols[i], pch = 16)
  }

  legend(
    "topright",
    legend = criterion,
    col = cols,
    lwd = 2,
    pch = 16,
    bty = "n"
  )
}
