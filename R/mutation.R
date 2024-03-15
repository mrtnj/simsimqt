


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

