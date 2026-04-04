#' Negative Binomial Math Utilities for scComBat-Seq
#'
#' Core statistical functions for NB distribution parameter mapping,
#' used internally by the scDesign3 simulation pipeline.

#' Map count values between two NB distributions using quantile matching
#'
#' @param counts_sub Matrix of observed counts (genes x cells)
#' @param old_mu Matrix of old NB means
#' @param old_phi Vector of old NB dispersions (length = nrow(counts_sub))
#' @param new_mu Matrix of new NB means
#' @param new_phi Vector of new NB dispersions
#' @return Matrix of remapped counts preserving distributional properties
match_quantiles <- function(counts_sub, old_mu, old_phi, new_mu, new_phi) {
  new_counts_sub <- matrix(NA, nrow = nrow(counts_sub), ncol = ncol(counts_sub))
  for (a in 1:nrow(counts_sub)) {
    for (b in 1:ncol(counts_sub)) {
      if (counts_sub[a, b] <= 1) {
        new_counts_sub[a, b] <- counts_sub[a, b]
      } else {
        tmp_p <- pnbinom(counts_sub[a, b] - 1, mu = old_mu[a, b], size = 1 / old_phi[a])
        if (abs(tmp_p - 1) < 1e-4) {
          # Outlier count: p≈1 would give Inf via qnbinom; preserve original
          new_counts_sub[a, b] <- counts_sub[a, b]
        } else {
          new_counts_sub[a, b] <- 1 + qnbinom(tmp_p, mu = new_mu[a, b], size = 1 / new_phi[a])
        }
      }
    }
  }
  return(new_counts_sub)
}

#' Monte Carlo integration over NB distributions for empirical Bayes shrinkage
#'
#' Estimates posterior gamma (additive batch effect) and phi (dispersion) for
#' each gene using a gene-subset weighted average over the empirical likelihood.
#'
#' @param dat Data matrix (genes x cells)
#' @param mu Mean matrix
#' @param gamma Vector of additive batch effects
#' @param phi Vector of dispersions
#' @param gene.subset.n Number of genes to sample for MC integration (NULL = use all, slow)
#' @return List with gamma_star, phi_star, and weights
monte_carlo_int_NB <- function(dat, mu, gamma, phi, gene.subset.n) {
  weights <- pos_res <- list()
  for (i in 1:nrow(dat)) {
    m <- mu[-i, !is.na(dat[i, ])]
    x <- dat[i, !is.na(dat[i, ])]
    gamma_sub <- gamma[-i]
    phi_sub <- phi[-i]

    if (!is.null(gene.subset.n) && is.numeric(gene.subset.n) && length(gene.subset.n) == 1) {
      if (i == 1) cat(sprintf("Using %s random genes for Monte Carlo integration\n", gene.subset.n))
      mcint_ind <- sample(1:(nrow(dat) - 1), gene.subset.n, replace = FALSE)
      m <- m[mcint_ind, ]
      gamma_sub <- gamma_sub[mcint_ind]
      phi_sub <- phi_sub[mcint_ind]
      G_sub <- gene.subset.n
    } else {
      if (i == 1) cat("Using all genes for Monte Carlo integration; slow for large gene sets\n")
      G_sub <- nrow(dat) - 1
    }

    LH <- sapply(1:G_sub, function(j) {
      prod(dnbinom(x, mu = m[j, ], size = 1 / phi_sub[j]))
    })
    LH[is.nan(LH)] <- 0

    if (sum(LH) == 0 || is.na(sum(LH))) {
      pos_res[[i]] <- c(gamma.star = as.numeric(gamma[i]), phi.star = as.numeric(phi[i]))
    } else {
      pos_res[[i]] <- c(
        gamma.star = sum(gamma_sub * LH) / sum(LH),
        phi.star   = sum(phi_sub * LH) / sum(LH)
      )
    }
    weights[[i]] <- as.matrix(LH / sum(LH))
  }
  pos_res <- do.call(rbind, pos_res)
  weights  <- do.call(cbind, weights)
  list(
    gamma_star = pos_res[, "gamma.star"],
    phi_star   = pos_res[, "phi.star"],
    weights    = weights
  )
}

#' Map dispersion parameters between distributions given a mean change
#'
#' Preserves variance in the original scale when mean shifts occur.
#'
#' @param old_mu Matrix of old means
#' @param new_mu Matrix of new means
#' @param old_phi Matrix of old dispersions
#' @param divider Matrix of divisors (e.g. fold-changes applied to mean)
#' @return Matrix of new dispersions
mapDisp <- function(old_mu, new_mu, old_phi, divider) {
  new_phi <- matrix(NA, nrow = nrow(old_mu), ncol = ncol(old_mu))
  for (a in 1:nrow(old_mu)) {
    for (b in 1:ncol(old_mu)) {
      old_var        <- old_mu[a, b] + old_mu[a, b]^2 * old_phi[a, b]
      new_var        <- old_var / (divider[a, b]^2)
      new_phi[a, b]  <- (new_var - new_mu[a, b]) / (new_mu[a, b]^2)
    }
  }
  return(new_phi)
}
