

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

  gv <- calculate_genetic_values(geno = geno[, simparam$traits[[1]]$loci_ix],
                                 a = simparam$traits[[1]]$a,
                                 d = simparam$traits[[1]]$d) +
    simparam$traits[[1]]$intercept

  pheno <- add_environmental_noise(gv, simparam$Ve[1])

  sex <- rep(NA_character_, nrow(geno))
  if (simparam$use_sexes) {
    sex <- rep(c("F", "M"), ceiling(nrow(geno)/2))
    sex <- sex[1:nrow(geno)]
  }

  simparam$last_id <- nrow(geno)

  structure(list(id = as.character(1:nrow(geno)),
                 sex = sex,
                 geno = geno,
                 gv = matrix(gv, ncol = 1),
                 pheno = matrix(pheno, ncol = 1)),
            class = "population")

}

#' Check the validity of a population object.
#'
#' @param pop A population object.
#'
#' @return The object, unchanged.
#' @export
validate_population <- function(pop) {

  stopifnot(is.character(pop$id))
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
