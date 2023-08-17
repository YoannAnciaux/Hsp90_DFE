###############################################################
# Functions for FGM DFE and mixture of FGM DFE and a gaussian #
# Written for documentation with roxygen2                     #
###############################################################

#' Distribution of fitness effect from Fisher's geometric Model
#'
#' @description Density, distribution function, quantile function and random
#' generation for the distribution of fitness effects of random mutations from \href{https://doi.org/10.1111/evo.12671}{Martin
#' & Lenormand (2015)}.
#'
#' @details The FGM DFE has a density of :
#' \deqn{f_x(x, n_dim, \lambda, s_o) = \frac{2}{\lambda}%
#' f_{\chi_{n}^2} (\frac{2 (s_o-x)}{\lambda}, \frac{2 s_o}{\lambda})}{%
#' fx(x, n, \lambda, so) = 2/\lambda d\chi[2,n](2 (so - x) / \lambda, 2 so / \lambda)}
#' The density is encoded using the functions from the package \code{\link[distr]{distr}}
#'
#' @param x,q vector of quantiles. The density is 0 if \eqn{s \ge so}.
#' @param p vector of probabilities
#' @param n number of observations
#' @param n_dim natural number. The dimensionality, i.e. the number of phenotypic
#' dimensions (traits) under selection.
#' @param lambda positive real number. Mutational variance per trait among random
#' mutations scaled by the strength of selection. It scales the mean fitness effect
#' of mutations.
#' @param so positive real number. Fitness distance between the wild-type (phenotype
#' without mutations) and the optimum.
#'
#' @return dfgm gives the density, pfgm gives the distribution function, qfgm
#' gives the quantile function, and rfgm generates random deviates.\n
#' The length of the result is determined by n for rfgm, and is the maximum of
#' the lengths of x,q or p for the other functions.
#'
#' @source Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#' across environments: Fisher's geometrical model with multiple optima.
#' Evolution, 69(6), 1433-1447.
#'
#' @name fgm
NULL
#> NULL

#' @rdname fgm
dfgm <- function(x, n_dim, lambda, so) {
  C <- distr::Chisq(df = n_dim, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::d(S)(x))
}

#' @rdname fgm
pfgm <- function(q, n_dim, lambda, so) {
  C <- distr::Chisq(df = n_dim, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::p(S)(q))
}

#' @rdname fgm
qfgm <- function(p, n_dim, lambda, so) {
  C <- distr::Chisq(df = n_dim, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::q(S)(p))
}

#' @rdname fgm
rfgm <- function(n, n_dim, lambda, so) {
  C <- distr::Chisq(df = n_dim, ncp = 2 * so / lambda)
  S <- so - lambda / 2 * C
  return(distr::r(S)(n))
}


#' Mixture of the distribution of fitness effect from Fisher's geometric Model
#' and a normal distribution
#'
#' @description Density, distribution function, quantile function and random
#' generation for the distribution of fitness effects of random mutations from \href{https://doi.org/10.1111/evo.12671}{Martin
#' & Lenormand (2015)}.
#'
#' @details The FGM DFE has a density of :
#' \deqn{f_x(x, n_dim, \lambda, s_o) = \frac{2}{\lambda}%
#' f_{\chi_{n}^2} (\frac{2 (s_o-x)}{\lambda}, \frac{2 s_o}{\lambda})}{%
#' fx(x, n, \lambda, so) = 2/\lambda d\chi[2,n](2 (so - x) / \lambda, 2 so / \lambda)}
#' The density is encoded using the functions from the package \code{\link[distr]{distr}}
#'
#' @param x,q vector of quantiles. The density is 0 if \eqn{s \ge so}.
#' @param p vector of probabilities
#' @param n number of observations
#' @param mean mean of the normal distribution.
#' @param sd standard deviation of the normal distribution.
#' @param n_dim natural number. The dimensionality, i.e. the number of phenotypic
#' dimensions (traits) under selection, in FGM.
#' @param lambda positive real number. Mutational variance per trait among random
#' mutations scaled by the strength of selection. It scales the mean fitness effect
#' of mutations, in FGM.
#' @param so positive real number. Fitness distance between the wild-type (phenotype
#' without mutations) and the optimum in FGM.
#' @param prob positive real numberbetween 0 and 1. Proportion of the distribution
#' belonging to the FGM. 1-prob is the porportion of the distribution belonging
#' to the gaussian.
#'
#' @return dnormfgmmix gives the density, pnormfgmmix gives the distribution function, qnormfgmmix
#' gives the quantile function, and rnormfgmmix generates random deviates.\n
#' The length of the result is determined by n for rnormfgmmix, and is the maximum of
#' the lengths of x,q or p for the other functions.
#'
#' @source Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#' across environments: Fisher's geometrical model with multiple optima.
#' Evolution, 69(6), 1433-1447.
#'
#' @name normfgmmix
NULL
#> NULL

#' @rdname normfgmmix
dnormfgmmix <- function(x, mean, sd, n_dim, lambda, so, prob) {
  if(n_dim > 0 & lambda > 0 & so > 0) {
    prob * dfgm(x, n_dim, lambda, so) + (1 - prob) * dnorm(x, mean, sd)
  } else {
    rep(NA, length(x))
  }
}

#' @rdname normfgmmix
pnormfgmmix <- function(q, mean, sd, n_dim, lambda, so, prob) {
  if(n_dim > 0 & lambda > 0 & so > 0) {
    prob * pfgm(q, n_dim, lambda, so) + (1 - prob) * pnorm(q, mean, sd)
  } else {
    rep(NA, length(q))
  }
}

#' @rdname normfgmmix
qnormfgmmix <- function(p, mean, sd, n_dim, lambda, so, prob){
  # https://stats.stackexchange.com/questions/390931/compute-quantile-function-from-a-mixture-of-normal-distribution
  # the min minmax below is computed to supply a range to the solver
  # the solution must be between the min and max
  # quantile of the mixed distributions
  if(n_dim > 0 & lambda > 0 & so > 0) {
    unlist(parallel::mclapply(p, function(pi) {
      minmax <- c(-10,#-.Machine$double.xmax,
                  min(c(qfgm(pi, n_dim, lambda, so),.Machine$double.xmax)));
      uniroot(function(q) pnormfgmmix(q, mean, sd, n_dim, lambda, so, prob) - pi,
              interval = minmax,
              tol = 10^{-16})$root},
      mc.cores = parallel::detectCores()-1))
  } else {
    rep(NA, length(q))
  }
}

#' @rdname normfgmmix
rnormfgmmix <- function(n, mean, sd, n_dim, lambda, so, prob) {
  ifelse(runif(n) < prob,
         rfgm(n, n_dim, lambda, so),
         rnorm(n, mean, sd))
}


#' Mixture of the distribution of fitness effect from Fisher's geometric Model
#' and a normal distribution with lambda proportional to the mean of the selection
#' coefficients of the mutants in the FGM DFE and not in the gaussian.
#'
#' @description Density, distribution function, quantile function and random
#' generation for the distribution of fitness effects of random mutations from \href{https://doi.org/10.1111/evo.12671}{Martin
#' & Lenormand (2015)}.
#'
#' @details The FGM DFE has a density of :
#' \deqn{f_x(x, n_dim, \lambda, s_o) = \frac{2}{\lambda}%
#' f_{\chi_{n}^2} (\frac{2 (s_o-x)}{\lambda}, \frac{2 s_o}{\lambda})}{%
#' fx(x, n, \lambda, so) = 2/\lambda d\chi[2,n](2 (so - x) / \lambda, 2 so / \lambda)}
#' The density is encoded using the functions from the package \code{\link[distr]{distr}}
#'
#' @param x,q vector of quantiles. The density is 0 if \eqn{s \ge so}.
#' @param p vector of probabilities
#' @param n number of observations
#' @param mean mean of the normal distribution.
#' @param sd standard deviation of the normal distribution.
#' @param n_dim natural number. The dimensionality, i.e. the number of phenotypic
#' dimensions (traits) under selection, in FGM.
#' @param mean_data mean of the selection coefficients of the data.
#' @param so positive real number. Fitness distance between the wild-type (phenotype
#' without mutations) and the optimum in FGM.
#' @param prob positive real numberbetween 0 and 1. Proportion of the distribution
#' belonging to the FGM. 1-prob is the porportion of the distribution belonging
#' to the gaussian.
#'
#' @return dnormfgmsmix gives the density, pnormfgmsmix gives the distribution function, qnormfgmsmix
#' gives the quantile function, and rnormfgmsmix generates random deviates.\n
#' The length of the result is determined by n for rnormfgmsmix, and is the maximum of
#' the lengths of x,q or p for the other functions.
#'
#' @source Martin, G., & Lenormand, T. (2015). The fitness effect of mutations
#' across environments: Fisher's geometrical model with multiple optima.
#' Evolution, 69(6), 1433-1447.
#'
#' @name normfgmsmix
NULL
#> NULL

#' @rdname normfgmsmix
dnormfgmsmix <- function(x, mean, sd, n_dim, so, prob, mean_data) {
  if(n_dim > 0 & so > 0) {
    lambda <- - 2 / n_dim * (mean_data - mean * (1-prob)) / prob;
    prob * dfgm(x, n_dim, lambda, so) + (1 - prob) * dnorm(x, mean, sd)
  } else {
    rep(NA, length(x))
  }
}

#' @rdname normfgmsmix
pnormfgmsmix <- function(q, mean, sd, n_dim, so, prob, mean_data) {
  if(n_dim > 0 & so > 0) {
    lambda <- - 2 / n_dim * (mean_data - mean * (1-prob)) / prob;
    prob * pfgm(q, n_dim, lambda, so) + (1 - prob) * pnorm(q, mean, sd)
  } else {
    rep(NA, length(q))
  }
}

#' @rdname normfgmsmix
qnormfgmsmix <- function(p, mean, sd, n_dim, so, prob, mean_data){
  # https://stats.stackexchange.com/questions/390931/compute-quantile-function-from-a-mixture-of-normal-distribution
  # the min minmax below is computed to supply a range to the solver
  # the solution must be between the min and max
  # quantile of the mixed distributions
  if(n_dim > 0 & so > 0) {
    lambda <- - 2 / n_dim * (mean_data - mean * (1-prob)) / prob;
    unlist(parallel::mclapply(p, function(pi) {
      minmax <- c(-10,#-.Machine$double.xmax,
                  min(c(qfgm(pi, n_dim, lambda, so),.Machine$double.xmax)));
      uniroot(function(q) pnormfgmmix(q, mean, sd, n_dim, lambda, so, prob) - pi,
              interval = minmax,
              tol = 10^{-16})$root},
      mc.cores = parallel::detectCores()-1))
  } else {
    rep(NA, length(q))
  }
}

#' @rdname normfgmsmix
rnormfgmsmix <- function(n, mean, sd, n_dim, so, prob, mean_data) {
  lambda <- - 2 / n_dim * (mean_data - mean * (1-prob)) / prob;
  ifelse(runif(n) < prob,
         rfgm(n, n_dim, lambda, so),
         rnorm(n, mean, sd))
}