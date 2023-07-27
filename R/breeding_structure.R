
#' Select individuals from a population based on trait
#'
#' @param pop Population to select from.
#' @param n_ind Number of individuals to select.
#' @param trait Trait number of the trait to select on.
#'
#' @return A population object of selected individuals.
#' @export
select_ind <- function(pop,
                       n_ind,
                       trait = 1) {

  ranking <- order(pop$pheno[, trait],
                   decreasing = TRUE)
  pop[pop$id[ranking[1:n_ind]]]

}


#' Select individuals from a population based on sampling weighted by
#' a fitness function.
#'
#' @param pop Population to select from.
#' @param n_ind Number of individuals to select.
#' @param trait Trait number of the trait to select on.
#' @param fitness_function A function that takes trait values and return
#' fitness values.
#'
#' @return A population object of selected individuals.
#' @export
select_ind_fitness <- function(pop,
                               n_ind,
                               fitness_function,
                               trait = 1) {

  fitness <- fitness_function(pop$pheno[, trait])
  average_fitness <- mean(fitness)

  stopifnot("Fitness function can't return NA" =
              !is.na(average_fitness))

  selected_ix <- sample(1:length(fitness),
                        n_ind,
                        prob = fitness/average_fitness,
                        replace = FALSE)

  pop[pop$id[selected_ix]]

}


get_gametes <- function(g) {
  stopifnot(length(g) == 1)
  if (g == 0) {
    gametes <- 0
  } else if (g == 1) {
    gametes <- c(0, 1)
  } else if (g == 2) {
    gametes <- 1
  }
  gametes
}


get_inherited_genotype <- function(mother_geno, father_geno) {
  offspring_geno <- NA_real_

  if (!is.na(mother_geno) & !is.na(father_geno)) {
    mother_gametes <- get_gametes_cpp(mother_geno)
    father_gametes <- get_gametes_cpp(father_geno)

    from_mother <- sample(mother_gametes, 1)
    from_father <- sample(father_gametes, 1)

    offspring_geno <- from_mother + from_father
  }
  offspring_geno
}


#' Make crosses according to a cross plan
#'
#' @param pop Population to take individuals from.
#' @param cross_plan Matrix of individual IDs, where the first column is
#' mothers and the second fathers.
#' @param n_progeny Number of progenies per cross.
#' @param simparam Simulation parameters.
#'
#' @return Population object of offspring.
#' @export
make_cross <- function(pop,
                       cross_plan,
                       n_progeny,
                       simparam) {

  cross_plan <- cbind(rep(cross_plan[,1],
                          each = n_progeny),
                      rep(cross_plan[,2],
                          each = n_progeny))

  offspring_geno <- matrix(0,
                           ncol = ncol(pop$geno),
                           nrow = nrow(cross_plan))

  for (cross_ix in 1:nrow(cross_plan)) {
    mother_geno <- pop$geno[pop$id == cross_plan[cross_ix, 1], ]
    father_geno <- pop$geno[pop$id == cross_plan[cross_ix, 2], ]

    offspring_geno[cross_ix, ] <- get_inherited_genotype_vector(mother_geno,
                                                                father_geno)
  }

  next_id <- simparam$last_id + 1
  simparam$last_id <- next_id + nrow(offspring_geno) - 1

  sex <- rep(c("F", "M"), ceiling(nrow(offspring_geno)/2))
  sex <- sex[1:nrow(offspring_geno)]

  gv <- calculate_genetic_values(geno = offspring_geno[, simparam$traits[[1]]$loci_ix],
                                 a = simparam$traits[[1]]$a,
                                 d = simparam$traits[[1]]$d) +
    simparam$traits[[1]]$intercept

  pheno <- add_environmental_noise(gv, simparam$Ve[1])

  offspring <- structure(list(id = as.character(next_id:(next_id + nrow(offspring_geno) - 1)),
                              sex = sex,
                              geno = offspring_geno,
                              gv = matrix(gv, ncol = 1),
                              pheno = matrix(pheno, ncol = 1)),
                         class = "population")

  offspring
}


#' Randomly mate individuals from a population while assigning the number
#' of offspring in a balanced way across parents and sexes.
#'
#' @param females Population object of mothers.
#' @param males Population object of fathers.
#' @param n_crosses Number of crosses.
#' @param n_progeny Number of offspring per cross.
#' @param simparam Simulation parameters object.
#'
#' @return Population object of offspring.
#' @export
balanced_cross <- function(females,
                           males,
                           n_crosses,
                           n_progeny,
                           simparam) {

  female_order <- sample(females$id,
                         length(females$id))
  female_order <- rep(female_order, length.out = n_crosses)

  male_order <- sample(males$id,
                       length(males$id))
  male_order <- rep(male_order, length.out = n_crosses)

  cross_plan <- cbind(female_order, male_order)

  make_cross(c(males, females),
             cross_plan,
             n_progeny,
             simparam)
}


#' Assign random crosses where each parent, within sex, is equally likely to be
#' chosen each time.
#'
#' @param females Population of mothers.
#' @param males Populations of fathers.
#' @param n_crosses Number of crosses.
#' @param n_progeny Number of offspring per cross.
#' @param simparam Simulation parameters object.
#'
#' @return Population object of offspring.
#' @export
rand_cross2 <- function(females,
                        males,
                        n_crosses,
                        n_progeny,
                        simparam) {

  stopifnot("Wrong sex in females." =
              all(females$sex == "F"))
  stopifnot("Wrong sex in males" =
              all(males$sex == "M"))

  female_order <- sample(females$id, n_crosses, replace = TRUE)
  male_order <- sample(males$id, n_crosses, replace = TRUE)

  cross_plan <- cbind(female_order, male_order)

  make_cross(c(males, females),
             cross_plan,
             n_progeny,
             simparam)
}


#' Assign random crosses where each parent, regardless of sex, is equally likely
#' to be chosen each time.
#'
#' @param pop Population of parents
#' @param n_crosses Number of crosses.
#' @param n_progeny Number of offspring per cross.
#' @param simparam Simulation parameters object.
#'
#' @return Population object of offspring.
#' @export
rand_cross <- function(pop,
                       n_crosses,
                       n_progeny,
                       simparam) {

  parent1 <- sample(pop$id, n_crosses, replace = TRUE)
  parent2 <- sample(pop$id, n_crosses, replace = TRUE)

  cross_plan <- cbind(parent1, parent2)

  make_cross(pop,
             cross_plan,
             n_progeny,
             simparam)
}



#' Cross individuals in a population based on sampling weighted by fitness
#'
#' @param pop Population of parents
#' @param n_crosses Number of crosses.
#' @param n_progeny Number of offspring per cross.
#' @param trait Trait number of the trait to select on.
#' @param fitness_function A function that takes trait values and return
#' fitness values.
#' @param simparam Simulation parameters object.
#'
#' @return A population object of selected individuals.
#' @export
select_cross_fitness <- function(pop,
                                 n_crosses,
                                 n_progeny,
                                 fitness_function,
                                 simparam,
                                 trait = 1,
                                 ...) {

  fitness <- fitness_function(pop$pheno[, trait], ...)
  average_fitness <- mean(fitness)

  stopifnot("Fitness function can't return NA" =
              !any(is.na(fitness)))

  stopifnot("Fitness must be a number" =
              all(is.numeric(fitness)))

  parent1_ix <- sample(1:length(fitness),
                       n_crosses,
                       prob = fitness/average_fitness,
                       replace = TRUE)
  parent2_ix <- sample(1:length(fitness),
                       n_crosses,
                       prob = fitness/average_fitness,
                       replace = TRUE)

  cross_plan <- cbind(pop$id[parent1_ix],
                      pop$id[parent2_ix])

  make_cross(pop,
             cross_plan,
             n_progeny,
             simparam)

}
