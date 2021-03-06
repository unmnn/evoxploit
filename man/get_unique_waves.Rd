% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxil.R
\name{get_unique_waves}
\alias{get_unique_waves}
\title{Extract unique wave numbers}
\usage{
get_unique_waves(att_names, suffix = "_s")
}
\arguments{
\item{att_names}{A character vector of attribute names.}

\item{suffix}{A string indicating the start of the wave index suffix.}
}
\value{
An integer vector of (sorted) wave numbers.
}
\description{
Extracts unique wave numbers (e.g. 0 for "som_bmi_s0") from a vector of
attribute names.
}
