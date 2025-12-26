#### This script defines auxiliary functions for the simulation ####

#### EStep ####
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

#### generate_startvalues ####
# generate random starting values for OpenMx
generate_startvalues <- function(n_clusters, n_factors, personmodel_list) {
  startvalues <- vector(mode = "list", length = n_clusters)
  labels_phi <- personmodel_list[[1]]$A$labels |> as.character()
  labels_phi <- labels_phi[!is.na(labels_phi)]
  labels_zeta <- personmodel_list[[1]]$Q$labels |> as.character() |> unique()
  labels_zeta <- labels_zeta[!is.na(labels_zeta)]
  for (k in 1:n_clusters) {
    # generate a stationary matrix of regression coefficients
    phistart <- matrix(runif(n_factors * n_factors, 0.05, .6), nrow = n_factors)
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
    zetastart <- matrix(runif(n_factors * n_factors, 0.3, 1.5),
                        nrow = n_factors)
    zetastartPD <- Matrix::nearPD(zetastart)$mat |> as.matrix()

    startvalues[[k]] <- c(
      c(phistart_scaled),
      zetastartPD[lower.tri(zetastartPD, diag = TRUE)]
    )
    names(startvalues[[k]]) <- c(labels_phi, unique(labels_zeta))
  }
  return(startvalues)
}

#### compute_observed_data_LL ####
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
