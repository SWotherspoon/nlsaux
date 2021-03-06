% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlsaux.r
\name{nlsInfluence}
\alias{nlsInfluence}
\title{NLS Influence Measures}
\usage{
nlsInfluence(object)
}
\arguments{
\item{object}{a fitted nls object.}
}
\value{
Returns a matrix with DFBETAS as columns labelled by the
  corresponding parameter names and DFFITS labelled by the name of
  the response, where each row corresponds to an observation.
}
\description{
Calculate DFBETAS and DFFITS for a fitted nls model.
}
