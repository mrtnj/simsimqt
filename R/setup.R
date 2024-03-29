

#' Sample founder genotypes with frequencies drawn from a random
#' distribution
#'
#' Create a genotype matrix from a founder population that can be used to start
#' the simulation.
#'
#' Loci are assumed to be unlinked (i.e., genotypes are drawn independently)
#' and biallelic (i.e., there are two alleles), and the organism is assumed to
#' be diploid (i.e, there are three genotypes). The genotypes are coded as
#' 0 (homozygote), 1 (heterozygote), and 2 (homzoygote).
#'
#' First, allele frequencies are drawn from a specified distribution for a
#' given number of loci, and then the genotypes are drawn from a binomial
#' distribution with size 2 and probability of success equal to the allele
#' frequency.
#'
#' @param n_ind Number of individuals.
#' @param n_loci Number of loci (note that they are not guaranteed to be polymorphic,
#' so this number should be bigger than the number of causative variants that
#' are planned).
#' @param distribution Distribution for allele frequencies; only beta is implemented.
#' @param parameters List of parameters for the random distribution for allele frequencies
#'  (e.g., for beta it would be list(shape1, shape2).
#'
#' @return Matrix of individual x locus genotypes coded 0, 1, 2
#' @examples
#' # Genotypes for 50 individuals, 1000 loci, allele frequencies from Beta(1/2, 1/2).
#' founder_genotypes <- draw_founder_genotypes(n_ind = 50,
#'                                             n_loci = 1000,
#'                                             distribution = "beta",
#'                                             parameters = list(shape1 = 1/2,
#'                                                               shape2 = 1/2))
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

    geno <- draw_founder_genotypes_freq(n_ind, freq)

  } else {
     stop("Distribution not implemented.")
  }


  geno
}


#' Sample founder genotypes with given frequencies
#'
#' This function is used internally by draw_founder_genotypes, which also draws
#' allele frequencies from a distribution. It may be used if you want more control
#' over starting allele frequencies.
#'
#' Loci are assumed to be unlinked (i.e., genotypes are drawn independently)
#' and biallelic (i.e., there are two alleles), and the organism is assumed to
#' be diploid (i.e, there are three genotypes). The genotypes are coded as
#' 0 (homozygote), 1 (heterozygote), and 2 (homzoygote).
#'
#' For each locus, genotypes are drawn from a binomial
#' distribution with size 2 and probability of success equal to the allele
#' frequency provided for that locus.
#'
#' @param n_ind Number of individuals.
#' @param freq Vector of base population frequencies.
#'
#' @return Matrix of individual x locus genotypes coded 0, 1, 2
#' @export
draw_founder_genotypes_freq <- function(n_ind,
                                        freq) {

  sapply(freq,
         function(p) stats::rbinom(n = n_ind, prob = p, size = 2))

}



#' Create a quantitative trait based on founder genotypes
#'
#' Function to set up causative loci and genetic architecture for a quantitative
#' trait. The results are returned as a trait list which can be used to make a
#' simulation parameters object.
#'
#' Polymorphic loci from the founder genotypes are randomly chosen to be
#' causative and assigned additive genetic coefficients and optionally
#' dominance coefficients based on a probability distribution.
#'
#' The resulting genetic variance and trait mean are adjusted to
#' desired numbers.
#'
#' @param founder_geno Matrix of individual x locus founder genotypes coded 0, 1, 2.
#' @param n_qtl Number of causative variants desired.
#' @param distribution Distribution for additive and dominance coefficients; only
#' "gaussian" is implemented.
#' @param parameters List of parameters with named components: mean_a and
#' var_a are the mean and variance for he additive genetic coefficients;
#' mean_dd and var_dd are the, optional, mean and variance for dominance degrees.
#' @param Vg Desired genetic variance. The default value is a trait with variance 1.
#' @param trait_mean Desired trait mean. The default value is trait standardised
#' to trait mean 0.
#'
#' @return A list describing the trait, with indices to the loci, and additive
#' genetic coefficients and dominance coefficients for each locus.
#' @examples
#' # Founder genotypes
#' founder_genotypes <- draw_founder_genotypes(n_ind = 50,
#'                                             n_loci = 1000,
#'                                             distribution = "beta",
#'                                             parameters = list(shape1 = 1/2,
#'                                                               shape2 = 1/2))
#' # A purely additive trait
#' trait <- make_trait(founder_genotypes,
#'                     n_qtl = 50,
#'                     distribution = "gaussian",
#'                     parameters = list(mean_a = 0, var_a = 1))
#' # A trait with dominance
#' trait <- make_trait(founder_genotypes,
#'                     n_qtl = 50,
#'                     distribution = "gaussian",
#'                     parameters = list(mean_a = 0, var_a = 1,
#'                                       mean_dd = 0.1, var_dd = 0.2))
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

  make_trait_manual(founder_geno,
                    loci_ix,
                    a,
                    d,
                    Vg,
                    trait_mean)
}


#' Create a trait based on given effect sizes
#'
#' This function is used internally by make_trait, and can be used if you
#' want more control of genetic architecture.
#'
#' Polymorphic loci from the founder genotypes are randomly chosen to be
#' causative and assigned additive genetic coefficients and optionally
#' dominance coefficients based on a probability distribution.
#'
#' The resulting genetic variance and trait mean are adjusted to
#' desired numbers.
#'
#' @param founder_geno Matrix of individual x locus founder genotypes coded 0, 1, 2.
#' @param loci_ix Indices of causative variants in genotype matrix.
#' @param a Vector of additive genetic coefficients.
#' @param d Vector of dominance coefficients.
#' @param Vg Desired genetic variance; NULL for no adjustment.
#' @param trait_mean Desired trait mean; NULL for no adjustment.
#'
#' @return A list describing the trait, with indices to the loci, and additive
#' genetic coefficients and dominance coefficients for each locus.
#' @export
#'
make_trait_manual <- function(founder_geno,
                              loci_ix,
                              a,
                              d,
                              Vg = NULL,
                              trait_mean = NULL) {

  if (is.null(Vg)) {
    scaling_constant <- 1
  } else {
    genetic_values <- calculate_genetic_values(founder_geno[, loci_ix],
                                               a,
                                               d)
    variance_initial <- stats::var(genetic_values)
    scaling_constant <- sqrt(Vg) / sqrt(variance_initial)
  }

  if (is.null(trait_mean)) {
    intercept <- 0
  } else {
    intercept <- trait_mean - mean(genetic_values) * scaling_constant
  }

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



#' @title Simulation parameters object
#'
#' @description
#' Contains data on traits, environmental variance and settings for a
#' simsimqt simulation.
#' @export
SimParam <- R6::R6Class("SimParam",
                        public = list(
                          traits = NULL,
                          Ve = NULL,
                          use_sexes = NULL,
                          last_id = NULL,
                          initialize = function(traits = traits,
                                                Ve = Ve,
                                                use_sexes = TRUE) {
                            stopifnot("Traits not in list" =
                                        is.list(traits))
                            self$traits <- traits
                            self$Ve <- Ve
                            self$use_sexes <- use_sexes
                            self$last_id <- 0
                          }
                        ))
