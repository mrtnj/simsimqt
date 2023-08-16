

#' Create a population object
#'
#' @param geno Matrix of genotypes.
#' @param simparam Simulation parameters object.
#'
#' @return
#' @export
new_population <- function(geno,
                           simparam) {

  stopifnot("Only one trait implemented yet. Go yell at Martin!" =
              length(simparam$trait) < 2)

  gv_pheno <- get_trait_values(geno, simparam)

  sex <- rep(NA_character_, nrow(geno))
  if (simparam$use_sexes) {
    sex <- rep(c("F", "M"), ceiling(nrow(geno)/2))
    sex <- sex[1:nrow(geno)]
  }

  simparam$last_id <- nrow(geno)

  structure(list(id = as.character(1:nrow(geno)),
                 sex = sex,
                 geno = geno,
                 gv = gv_pheno$gv,
                 pheno = gv_pheno$pheno),
            class = "population")

}


get_trait_values <- function(geno, simparam) {

  if (length(simparam$traits) > 0) {
    gv <- calculate_genetic_values(geno = geno[, simparam$traits[[1]]$loci_ix],
                                   a = simparam$traits[[1]]$a,
                                   d = simparam$traits[[1]]$d) +
      simparam$traits[[1]]$intercept

    gv <- matrix(gv, ncol = 1)

    pheno <- add_environmental_noise(gv, simparam$Ve[1])

    pheno <- matrix(pheno, ncol = 1)

  } else {
    gv <- matrix(ncol = 0, nrow = nrow(geno))
    pheno <- matrix(ncol = 0, nrow = nrow(geno))
  }

  list(gv = gv, pheno = pheno)
}

#' Check the validity of a population object.
#'
#' @param pop A population object.
#'
#' @return The object, unchanged.
#' @export
validate_population <- function(pop) {

  stopifnot(is.character(pop$id))
  stopifnot(!any(duplicated(pop$id)))
  stopifnot(is.character(pop$sex))
  stopifnot(is.matrix(pop$geno))
  stopifnot(is.matrix(pop$gv))
  stopifnot(is.matrix(pop$pheno))
  stopifnot(is.numeric(pop$gv))
  stopifnot(is.numeric(pop$pheno))

  pop
}


#' Hard bracket indexing on a population object takes a subset of individuals.
#'
#' @param pop A population object
#' @param i Numeric indices of individuals to select.
#'
#' @return A population object containing the subset.
#' @export
`[.population` <- function(pop, ids) {

  stopifnot("IDs not found in population" =
              all(ids %in% pop$id))

  i <- match(ids, pop$id)

  structure(list(id = pop$id[i],
                 sex = pop$sex[i],
                 geno = pop$geno[i, , drop = FALSE],
                 gv = pop$gv[i, , drop = FALSE],
                 pheno = pop$pheno[i, , drop = FALSE]),
            class = "population")

}

#' Concatenate two populations
#'
#' @param pop1 First population
#' @param pop2 Second population
#'
#' @return Combined population
#' @export
c.population <- function(pop1, pop2) {

  stopifnot("Repeated IDs in combined population." =
              length(intersect(pop1$id, pop2$id)) == 0)

  structure(list(id = c(pop1$id, pop2$id),
                 sex = c(pop1$sex, pop2$sex),
                 geno = rbind(pop1$geno, pop2$geno),
                 gv = rbind(pop1$gv, pop2$gv),
                 pheno = rbind(pop1$pheno, pop2$pheno)),
            class = "population")

}


#' Print a nice summary of a population object.
#'
#' @param pop
#'
#' @return Nothing.
#' @export
print.population <- function(pop) {

  cat("Population with", nrow(pop$geno), "individuals and", ncol(pop$geno), "loci.\n")

  cat("First IDs and sexes:\n")
  print(head(pop$id))
  print(head(pop$sex))

  cat("First genetic values:\n")
  print(head(pop$gv))

  cat("First trait values:\n")
  print(head(pop$pheno))

}


#' Extract females from a population.
#'
#' @param pop A population object.
#'
#' @return Population of females.
#' @export
get_females <- function(pop) {

  pop[pop$id[pop$sex == "F"]]

}

#' Extract males from a population.
#'
#' @param pop A population object.
#'
#' @return Population of males.
#' @export
get_males <- function(pop) {

  pop[pop$id[pop$sex == "M"]]

}


#' Get mean genetic value from population object.
#'
#' @param pop Population object.
#'
#' @return Mean.
#' @export
meanG <- function(pop) {

  mean(pop$gv)

}

#' Get genetic variance from population object
#'
#' @param pop Population object.
#'
#' @return Variance.
#' @export
varG <- function(pop) {

  var(pop$gv)

}


#' Get genotypes at causative variants from population object
#'
#' @param pop Population object.
#' @param trait Trait number.
#' @param simparam Simulation parameters object.
#'
#' @return Genotype matrix coded 0, 1, 2.
#' @export
pull_qtl_geno <- function(pop,
                          trait = 1,
                          simparam) {

  pop$geno[, simparam$traits[[trait]]$loci_ix]

}


#' Get frequencies of causative variants from population object
#'
#' @param pop Population object.
#' @param trait Trait number.
#' @param simparam Simulation parameters object.
#'
#' @return Allele frequencies
#' @export
pull_qtl_freq <- function(pop,
                          trait = 1,
                          simparam) {

  geno <- pop$geno[, simparam$traits[[trait]]$loci_ix]
  colSums(geno)/2/nrow(geno)

}
