% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxil.R
\name{evo_permutations}
\alias{evo_permutations}
\title{Find pairwise permutations of attribute names}
\usage{
evo_permutations(int_waves)
}
\arguments{
\item{int_waves}{An integer vector of wave indices.}
}
\value{
A list of pairs of wave indices
}
\description{
Find all possible pairwise wave permuations of attribute names, e.g.
{(som_bmi_s0, som_bmi_s1), (som_bmi_s0, som_bmi_s2),
(som_bmi_s1, som_bmi_s2)}. Autmatically extracts the set of wave suffixes.
}
