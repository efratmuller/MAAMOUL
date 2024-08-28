#' Use functions from the BioNet package to fit a beta-uniform mixture model
#' The df table is expected to have a "pval" column and a "feature" column
#'
#' @NoRd
fit_bum_model <- function(df, node_type, plot_dir, EPS = 0.000000001) {
  require(BioNet)
  pvals <- df$pval
  pvals <- sapply(pvals, function(p, EPS) max(p, EPS), EPS) # P-values cannot be zeros, so add epsilon
  names(pvals) <- df$feature

  # Save a plot of BUM model fit
  png(filename = file.path(plot_dir, paste(node_type, 'bum_fit.png', sep='_')), width = 700, height = 350)
  bum <- fitBumModel(pvals, starts=100)
  dev.off()

  # Save a plot of log likelihood surface
  # png(filename = file.path(plot_dir, paste(node_type, 'll_surface.png', sep='_')), width = 450, height = 350)
  # plotLLSurface(pvals, bum)
  # dev.off()

  return(bum)
}

#' Randomly sample a value given BUM parameters
#'
#' @NoRd
sample_from_bum <- function(a, lambda, N = 1) {
  beta_shape1 <- c(a, 1)
  # beta_shape1[1] is the beta distribution,
  #  beta_shape1[2] is the uniform distribution

  # Sample which mixture model component (noise/signal) each variable comes from
  mixt_comp <- sample(1:2, prob=c(1-lambda, lambda), size = N, replace = TRUE)

  # Then, sample the actual p-value
  sampled_pval <- rbeta(N, shape1 = beta_shape1[mixt_comp], shape2 = 1)
  return(sampled_pval)
}
