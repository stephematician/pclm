% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rpclmfit.R
\name{rpclmfit}
\alias{rpclmfit}
\title{Randomly draw samples from fitted PCLM}
\usage{
rpclmfit(D, pclm_fit, n_sample = 1L)
}
\arguments{
\item{D}{a two-column data.frame (or matrix) of integers, each row represents
an observation of an interval of bins, columns correspond to lower and
upper bins of the interval (inclusive)}

\item{pclm_fit}{a list as returned by pclmfit}

\item{n_sample}{number of samples to draw for each observation}
}
\value{
matrix of drawn bins with one row per observation, and columns
correspond to samples
}
\description{
Generate a sample for each observation (or otherwise) of a variable taken
from an interval of a PCLM.
}
\details{
Define any interval of bins (not necessarily the same as the original
groupings) for each observation and then draw a bin from the PCLM using the
probabilities obtained from \code{pclmfit}.

The flexibility to define any interval does come at a small cost in terms of
performance. This is most useful when incorporating some kind of censoring
into a prediction that provides additional information to narrow down the
range of possible bins. A simple example might be when 'age at survey' is
grouped, but 'age at event' is not, and therefore the participant must be at
least the latter age (in bins).
}
