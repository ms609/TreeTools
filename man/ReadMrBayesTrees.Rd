% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ReadMrBayes.R
\name{ReadMrBayesTrees}
\alias{ReadMrBayesTrees}
\alias{ReadMrBayes}
\alias{MrBayesTrees}
\title{Read posterior tree sample produced by MrBayes}
\usage{
ReadMrBayesTrees(filepath, n = NULL, burninFrac = NULL)

ReadMrBayes(filepath, n = NULL, burninFrac = NULL)

MrBayesTrees(filepath, n = NULL, burninFrac = NULL)
}
\arguments{
\item{filepath}{character string specifying path to \code{.nex} input file used
to initialize the MrBayes analysis,
relative to the R working directory (visible with \code{getwd()}).}

\item{n}{Integer specifying number of trees to sample from posterior.}

\item{burninFrac}{Fraction of trees to discard from each run as burn-in.
If \code{NULL} (the default), this will be read from the last \code{mcmc} or \code{mcmcp}
command in \code{filepath}.}
}
\value{
\code{ReadMrBayesTrees()} returns a 'multiPhylo' object containing
\code{n} trees sampled evenly from all runs generated by analysis of \code{filepath},
or \code{NULL} if no trees are found.
}
\description{
Read posterior trees from \href{https://nbisweden.github.io/MrBayes/}{'MrBayes'} output files, discarding burn-in
generations.
}
\details{
\code{ReadMrBayesTrees()} samples trees from the posterior distributions
computed using the Bayesian inference software \href{https://nbisweden.github.io/MrBayes/}{'MrBayes'}
}
\examples{
\dontrun{ # Download will take a few seconds
  url <- 
  "https://raw.githubusercontent.com/ms609/hyoliths/master/MrBayes/hyo.nex"
  trees <- ReadMrBayesTrees(url, n = 40)
  plot(Consensus(trees, p = 0.5))
}
}
\seealso{
Other tree import functions: 
\code{\link{ReadTntTree}()}
}
\author{
\href{https://orcid.org/0000-0001-5660-1727}{Martin R. Smith}
(\href{mailto:martin.smith@durham.ac.uk}{martin.smith@durham.ac.uk})
}
\concept{tree import functions}
