\name{savFormat}
\alias{savFormat}
\title{Base class for formatters}
\description{
  Defines the necessary slots to create parse different
  binary files using the same generic parser.
}
\section{Slots}{
  \describe{ \item{\code{name}:}{vector of column names}
  \item{\code{type}:}{vector of data types of elements}
  \item{\code{lengths}:}{vector of byte lengths for each
  element} \item{\code{order}:}{vector of column names for
  sorting} }
}

