#' Title
#'
#' @param n_factors
#' @param random_intercept
#' @param free
#' @param startvalues
#' @param labels
#' @param lbound
#' @param ubound
#'
#' @return
#' @export
#'
#' @examples
create_Q <- function(step2output, startvalues,
                     random_intercept = FALSE, free = NULL,
                     labels = NULL, lbound = NA, ubound = NA){

  factors <- step2output$other$factors
  n_factors <- length(factors)

  if(!random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }

    # create label matrix
    if(is.null(labels)){
      labels <- matrix(NA, n_factors, n_factors)

      # Fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels[i, j] <- paste0("zeta_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
    }

    # create OpenMx model object
    Q <- OpenMx::mxMatrix(type = "Full", name = "Q",
                          nrow = n_factors, ncol = n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

  if(random_intercept){
    # create matrix that indicates free parameters
    if(is.null(free)){
      free <- matrix(TRUE, nrow = n_factors, ncol = n_factors)
    }
    others <- matrix(FALSE, nrow = n_factors, ncol = n_factors)
    free <- rbind(cbind(free, others),
                  cbind(others, others))

    # create label matrix
    if(is.null(labels)){
      labels <- matrix(NA, n_factors, n_factors)

      # Fill the matrix with symmetric labels
      for (i in 1:n_factors) {
        for (j in 1:n_factors) {
          labels[i, j] <- paste0("zeta_", factors[min(i, j)], "_", factors[max(i, j)])
        }
      }
    }
    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, others),
                    cbind(others, others))

    # create (start)value matrix
    others <- matrix(0, nrow = n_factors, ncol = n_factors)
    startvalues <- rbind(cbind(startvalues, others),
                    cbind(others, others))


    # create OpenMx model object
    Q <- OpenMx::mxMatrix(type = "Full", name = "Q",
                          nrow = 2*n_factors, ncol = 2*n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

  return(Q)
}

# to do: add lbound and ubound to random intercept part
