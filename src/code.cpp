#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_gametes_cpp(double g) {
  NumericVector v = {0, 0};
  if (g == 1) {
    v[1] = 1;
  } else if (g == 2) {
    v[1] = 1;
    v[0] = 1;
  }
  return v;
}


// [[Rcpp::export]]
NumericVector get_inherited_genotype_vector (NumericVector mother_geno,
                                             NumericVector father_geno) {
  int n_loci = mother_geno.length();
  NumericVector offspring_geno(n_loci);

  for (int locus_ix = 0; locus_ix < n_loci; locus_ix++) {
    NumericVector mother_gametes = get_gametes_cpp(mother_geno[locus_ix]);
    NumericVector father_gametes = get_gametes_cpp(father_geno[locus_ix]);

    double from_mother = sample(mother_gametes, 1)[0];
    double from_father = sample(father_gametes, 1)[0];

    offspring_geno[locus_ix] = from_mother + from_father;
  }

  return offspring_geno;
}



// get_inherited_genotype <- function(mother_geno, father_geno) {
//   offspring_geno <- NA_real_
//
//   if (!is.na(mother_geno) & !is.na(father_geno)) {
//     mother_gametes <- get_gametes_cpp(mother_geno)
//     father_gametes <- get_gametes_cpp(father_geno)
//
//     from_mother <- sample(mother_gametes, 1)
//     from_father <- sample(father_gametes, 1)
//
//     offspring_geno <- from_mother + from_father
//   }
//   offspring_geno
// }
