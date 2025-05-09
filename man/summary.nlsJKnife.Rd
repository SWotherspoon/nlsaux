% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlsaux.r
\name{summary.nlsJKnife}
\alias{summary.nlsJKnife}
\alias{print.summary.nlsJKnife}
\title{Summary for nlsJKnife}
\usage{
\method{summary}{nlsJKnife}(object, level = 0.95, ...)

\method{print}{summary.nlsJKnife}(x, digits = max(3L, getOption("digits") - 3L), ...)
}
\arguments{
\item{object}{An object of class `nlsJKnife`.}

\item{level}{The confidence level required.}

\item{...}{Currently ignored.}

\item{x}{An object of class `summary.nlsJKnife`.}

\item{digits}{The number of significant digits to use when printing.}
}
\value{
Returns an object of class `nlsJKnife` containing

- `object`: the `nlsJKnife` object
- `theta`: estimate of the test statistic
- `corrected`: jackknife bias corrected estimate
- `bias`: jackknife estimate of bias
- `cov`: jackknife estimate of covariance
- `level`: requested confidence level
- `lwr`: jackknife lower confidence limit
- `upr`: jackknife upper confidence limit
}
\description{
Constructs a table of estimates and confidence intervals for
the test statistic.
}
