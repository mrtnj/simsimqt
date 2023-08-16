


#' Mutate the genotype matrix of a population by a per-locus mutation rate.
#'
#' @param pop Population object.
#' @param mutation_rate Mutation rate per locus and genome copy.
#'
#' @return Population object.
#' @export
mutate_per_locus <- function(pop,
                             mutation_rate) {

  pop$geno <- mutate_genotypes_per_locus(pop$geno,
                                         mutation_rate)

  pop

}

