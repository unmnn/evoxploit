% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auxil.R
\name{get_global_attname_stem}
\alias{get_global_attname_stem}
\title{Extract (unique) stems from attribute names}
\usage{
get_global_attname_stem(
  att_names,
  min_wave_count = length(get_unique_waves(att_names, suffix)),
  suffix = "_s"
)
}
\arguments{
\item{att_names}{A character vector of attribute names.}

\item{min_wave_count}{The minimum number of waves an attribute must be present
in for its name stem to be returned. If not set, the parameter defaults to
the total number of waves.}

\item{suffix}{A string indicating the start of the wave index suffix.}
}
\value{
A character vector of attribute name stems.
}
\description{
Extracts unique stems from attribute names. Returns only the stems for
attributes which occur in at least \code{min_wave_count} waves.
}
