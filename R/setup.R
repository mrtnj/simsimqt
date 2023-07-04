

#' Sample founder genotypes with frequencies drawn from a random
#' distribution.
#'
#' @param n_ind Number of individuals.
#' @param n_loci Number of loci (note that they are not guaranteed to be polymorphic,
#' so this number should be bigger than the number of causative variants that
#' are planned)
#' @param distribution Distribution for allele frequencies; only beta is implemented.
#' @param parameters List of parameters for the random distribution (e.g., for beta
#' it would be list(shape1, shape2).
#'
#' @return Matrix of individual x locus genotypes coded 0, 1, 2
#' @export
#'
draw_founder_genotypes <- function(n_ind,
                                   n_loci,
                                   distribution,
                                   parameters) {

  if (distribution == "beta") {

    freq <- stats::rbeta(n = n_loci,
                         shape1 = parameters$shape1,
                         shape2 = parameters$shape2)

    geno <- sapply(freq,
                   function(p) stats::rbinom(n = n_ind, prob = p, size = 2))

  } else {
     stop("Distribution not implemented.")
  }


  geno
}



#' Create a quantitative trait based on founder genotypes.
#'
#' @param founder_geno Matrix of individual x locus founder genotypes coded 0, 1, 2.
#' @param n_qtl Number of causative variants desired.
#' @param distribution Distribution for additive and dominance coefficients; only
#' "gaussian" is implemented.
#' @param parameters List of parameters with named components: mean_a and
#' var_a are the mean and variance for he additive genetic coefficients;
#' mean_dd and var_dd are the, optional, mean and variance for dominance degrees.
#' @param Vg Desired genetic variance.
#'
#' @return A list describing the trait, with indices to the loci, and additive
#' genetic coefficients and dominance coefficients for each locus.
#' @export
#'
make_trait <- function(founder_geno,
                       n_qtl,
                       distribution,
                       parameters,
                       Vg = 1,
                       trait_mean = 0) {

  loci_ix <- pick_loci(founder_geno, n_qtl)

  if (distribution == "gaussian") {

    a <- stats::rnorm(n_qtl,
                      mean = parameters$mean_a,
                      sd = sqrt(parameters$var_a))

    if ("mean_dd" %in% names(parameters)) {

      dd <- stats::rnorm(n_qtl,
                         mean = parameters$mean_dd,
                         sd = sqrt(parameters$var_dd))

      d <- abs(a) * dd

    } else {
      d <- rep(0, n_qtl)
    }

  } else {
    stop("Distribution not implemented.")
  }

  genetic_values <- calculate_genetic_values(founder_geno[, loci_ix],
                                             a,
                                             d)
  variance_initial <- stats::var(genetic_values)
  scaling_constant <- sqrt(Vg) / sqrt(variance_initial)

  intercept <- trait_mean - mean(genetic_values) * scaling_constant

  list(loci_ix = loci_ix,
       a = a * scaling_constant,
       d = d * scaling_constant,
       intercept = intercept)
}



#' Pick candidate loci to be causative variants from a genotype matrix.
#'
#' @param geno Matrix of individual x locus genotypes coded 0, 1, 2.
#' @param n_qtl Number of causative variants desired.
#'
#' @return Indices to selected loci.
#' @export
pick_loci <- function(geno, n_qtl) {

  f <- colSums(geno)/nrow(geno)/2
  candidate_ix <- which(f > 0 & f < 1)

  stopifnot("Too few polymorphic loci" =
              length(candidate_ix) >= n_qtl)

  chosen <- sample(1:length(candidate_ix),
                   n_qtl)

  candidate_ix[chosen]
}



#' Get genetic values from genotypes and locus-level coefficients.
#'
#' @param geno Matrix of individual x locus genotypes coded 0, 1, 2.
#' @param a Vector of additive genetic coefficients.
#' @param d Vector of dominance coefficients.
#'
#' @return Vector of genetic values.
#' @export
calculate_genetic_values <- function(geno, a, d) {

  stopifnot("Number of loci differs from number of additive coefficients" =
              ncol(geno) == length(a))
  stopifnot("Number of loci differs from number of dominance coefficients" =
              ncol(geno) == length(d))

  n_ind <- nrow(geno)
  n_loc <- ncol(geno)

  geno_additive <- geno - 1

  genetic_value <- numeric(n_ind)

  for (locus_ix in 1:n_loc) {

    locus_genetic_value <- geno_additive[, locus_ix] * a[locus_ix] +
      as.numeric(geno_additive[, locus_ix] == 0) * d[locus_ix]

    genetic_value <- genetic_value + locus_genetic_value

  }

  genetic_value
}


#' Add environmental variation to get trait values.
#'
#' @param genetic_value Genetic values for individuals.
#' @param Ve Desired environmental variance.
#'
#' @return Vector of trait values.
#' @export
add_environmental_noise <- function(genetic_value,
                                    Ve) {

  n_ind <- length(genetic_value)

  pheno <- genetic_value + stats::rnorm(n_ind, 0, sqrt(Ve))

  pheno

}



#' Create the simulation parameters object. Currently very uninteresting.
#'
#' @param traits List of traits.
#' @param Ve Vector of environmental variances if intented to be constant.
#'
#' @return Simulation parameters object.
#' @export
# make_simparam <- function(traits,
#                           Ve,
#                           use_sexes = TRUE) {
#
#   list(traits = traits,
#        Ve = Ve,
#        use_sexes = use_sexes,
#        last_id = 0)
#
# }

SimParam <- R6::R6Class("SimParam",
                        public = list(
                          traits = NULL,
                          Ve = NULL,
                          use_sexes = NULL,
                          last_id = NULL,
                          initialize = function(traits = traits,
                                                Ve = Ve,
                                                use_sexes = TRUE) {
                            self$traits <- traits
                            self$Ve <- Ve
                            self$use_sexes <- use_sexes
                            self$last_id <- 0
                          }
                        ))
