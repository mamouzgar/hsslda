% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/separationFunctions.R
\name{getScore_custom}
\alias{getScore_custom}
\title{getScore_custom}
\usage{
getScore_custom(x, y, cols, custom.score.method, ...)
}
\arguments{
\item{x}{dataframe of training data}

\item{y}{vector of class labels matching training data rows}

\item{cols}{vector of column names}

\item{custom.score.method}{a custom-function that takes in x, y, and cols, and outputs a score where the larger the value, the more optimal your separation criteria is.}
}
\description{
getScore_custom: placeholder for a custom metric
}
