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



// [[Rcpp::export]]
int mutate_geno (int geno,
                 bool allele1_mutates,
                 bool allele2_mutates) {
  int geno_after_mutation = geno;

  if (geno == 0) {
    // Wildtype homozygote
    if (allele1_mutates && allele2_mutates) {
      geno_after_mutation = 2;
    } else if (allele1_mutates || allele2_mutates) {
      geno_after_mutation = 1;
    }
  } else if (geno == 2) {
    // Mutant homozygote
    if (allele1_mutates && allele2_mutates) {
      geno_after_mutation = 0;
    } else if (allele1_mutates || allele2_mutates) {
      geno_after_mutation = 1;
    }
  } else if (geno == 1) {
    // Heterozygote
    if (allele1_mutates && allele2_mutates) {
      geno_after_mutation = 1;
    } else if (allele1_mutates || allele2_mutates) {
      // Choose which allele mutates
      double allele_draw = R::runif(0, 1);
      if (allele_draw < 0.5) {
        // 0 -> 1 means genotype 1 -> 2
        geno_after_mutation = 2;
      } else {
        // 1 -> 0 means genotype 1 -> 0
        geno_after_mutation = 0;
      }
    }

  }

  return geno_after_mutation;

}



// [[Rcpp::export]]
NumericMatrix mutate_genotypes_per_locus (NumericMatrix geno,
                                          double mutation_rate) {

  int n_loci = geno.ncol();
  int n_ind = geno.nrow();

  for (int locus_ix = 0; locus_ix < n_loci; locus_ix++) {
    for (int ind_ix = 0; ind_ix < n_ind; ind_ix++) {

      NumericVector mutation_draws = Rcpp::runif(2, 0, 1);
      bool allele1_mutates = mutation_draws[0] < mutation_rate;
      bool allele2_mutates = mutation_draws[1] < mutation_rate;

      int element = geno(ind_ix, locus_ix);
      geno(ind_ix, locus_ix) = mutate_geno(element,
           allele1_mutates,
           allele2_mutates);

    }
  }

  return geno;
}

