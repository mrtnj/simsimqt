


#' Mutate the genotype matrix of a population by a per-locus mutation rate.
#'
#' This function introduces mutations to a population object. Mutations are
#' assumed to happen between the two possible alleles with the same probabilty
#' for forward and backward mutation. The mutation rate is per locus.
#'
#' @param pop Population object.
#' @param mutation_rate Mutation rate per locus and genome copy.
#'
#' @return Population object with mutated alleles changed.
#' @export
mutate_per_locus <- function(pop,
                             mutation_rate) {

  pop$geno <- mutate_genotypes_per_locus(pop$geno,
                                         mutation_rate)

  pop

}



#' Mutate the genotype matrix of a population by introducing new mutations
#'
#' This function introduces mutations to a population object. Mutations are assumed
#' to always introduce a new locus, with a causative allele in one copy. There
#' is no backward mutation. The mutation's additive effect is drawn from a normal
#' distribution with given variance. The dominance coefficient is always zero.
#'
#' @param pop Population object.
#' @param simparam Simulation parameters object.
#' @param mutation_rate Mutation rate per individual and generation.
#' @param var_a Variance of the additive genetic coefficient.
#' @param trait Trait number to be mutated. For now, always 1.
#'
#' @return List containing population object (first component) and simparam
#' object (second component) with mutated alleles changed and trait information
#' updated.
#' @export
mutate_new_locus <- function(pop,
                             simparam,
                             mutation_rate,
                             var_a,
                             trait = 1) {

  new_pop <- pop
  new_simparam <- simparam

  n_ind <- length(pop$id)
  mutation_draw <- runif(1, 0, 1)

  if (mutation_draw < mutation_rate * n_ind) {
    new_locus <- sample(c(1, rep(0, n_ind - 1)))
    new_pop$geno <- cbind(pop$geno, as.matrix(new_locus))

    new_a <- rnorm(1, 0, var_a)

    new_simparam$traits[[trait]]$a <- c(simparam$traits[[trait]]$a, new_a)
    new_simparam$traits[[trait]]$d <- c(simparam$traits[[trait]]$d, 0)
    new_simparam$traits[[trait]]$loci_ix <- c(simparam$traits[[trait]]$loci_ix,
                                              ncol(pop$geno))

  }

  list(pop = new_pop, simparam = new_simparam)
}

