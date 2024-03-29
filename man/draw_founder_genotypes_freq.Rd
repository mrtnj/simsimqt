% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup.R
\name{draw_founder_genotypes_freq}
\alias{draw_founder_genotypes_freq}
\title{Sample founder genotypes with given frequencies}
\usage{
draw_founder_genotypes_freq(n_ind, freq)
}
\arguments{
\item{n_ind}{Number of individuals.}

\item{freq}{Vector of base population frequencies.}
}
\value{
Matrix of individual x locus genotypes coded 0, 1, 2
}
\description{
This function is used internally by draw_founder_genotypes, which also draws
allele frequencies from a distribution. It may be used if you want more control
over starting allele frequencies.
}
\details{
Loci are assumed to be unlinked (i.e., genotypes are drawn independently)
and biallelic (i.e., there are two alleles), and the organism is assumed to
be diploid (i.e, there are three genotypes). The genotypes are coded as
0 (homozygote), 1 (heterozygote), and 2 (homzoygote).

For each locus, genotypes are drawn from a binomial
distribution with size 2 and probability of success equal to the allele
frequency provided for that locus.
}
