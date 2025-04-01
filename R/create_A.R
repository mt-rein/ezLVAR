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
create_A <- function(step2output, startvalues,
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
      labels <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
    }

    # create OpenMx model object
    A <- OpenMx::mxMatrix(type = "Full", name = "A",
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
      labels <- outer(factors, factors,
                      FUN = function(i, j) paste0("phi_", i, "_", j))
    }
    others <- matrix(NA, nrow = n_factors, ncol = n_factors)
    labels <- rbind(cbind(labels, others),
                    cbind(others, others))

    # create (start)value matrix
    startvalues <- rbind(cbind(startvalues, matrix(0, nrow = n_factors, ncol = n_factors)),
                         cbind(matrix(0, nrow = n_factors, ncol = n_factors), diag(n_factors)))


    # create OpenMx model object
    A <- OpenMx::mxMatrix(type = "Full", name = "A",
                          nrow = 2*n_factors, ncol = 2*n_factors,
                          free = free,
                          values = startvalues,
                          labels = labels,
                          lbound = lbound,
                          ubound = ubound,
                          byrow = TRUE)
  }

  return(A)
}

# to do: add lbound and ubound to random intercept part
