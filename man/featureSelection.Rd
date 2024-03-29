% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/featureSelection.R
\name{featureSelection}
\alias{featureSelection}
\title{Feature Selection}
\usage{
featureSelection(features, biosigndata, model = 1, scoremin = "A")
}
\arguments{
\item{features}{An ExpressionSet or SummarizedExperiment object}

\item{biosigndata}{A biosign object}

\item{model}{A numeric/string vector indicating the models from which the scores will be used.}

\item{scoremin}{The minimum score to select the features "A", "S", "B".}
}
\value{
Returns a filtered ExpressionSet or SummarizedExperiment object.
}
\description{
Filters features depending on it's score in a biosign model
}
\author{
Jordi Rofes Herrera
}
