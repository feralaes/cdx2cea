% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05b_deterministic_analysis_functions.R
\name{generate_params_bounds}
\alias{generate_params_bounds}
\title{Generate lower and upper bound of the CEA parameters}
\usage{
generate_params_bounds(l_params_all)
}
\arguments{
\item{l_params_all}{List with all parameters of cost-effectiveness model.}

\item{seed}{Seed for reproducibility of Monte Carlo sampling.}
}
\value{
A list with the lower and upper bounds of the model parameters.
}
\description{
\code{generate_params_bounds} generates the lower and upper bounds of the
model parameters.
}
\examples{
generate_params_bounds()
}
