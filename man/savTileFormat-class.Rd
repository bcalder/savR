\docType{class}
\name{savTileFormat-class}
\alias{savTileFormat-class}
\title{Tile Metrics formatter}
\description{
Lane, tile, code, value.  Codes are:
}
\details{
\tabular{ll}{ 100 \tab Cluster Density \cr 101 \tab PF
Cluster Density \cr 102 \tab Number of clusters \cr 103
\tab Number of PF clusters \cr 400 \tab Control lane \cr }
}
\section{Slots}{
  \describe{ \item{\code{name}:}{vector of column names}
  \item{\code{type}:}{vector of data types of elements}
  \item{\code{lengths}:}{vector of byte lengths for each
  element} \item{\code{order}:}{vector of column names for
  sorting} \item{\code{version}:}{integer version number} }
}

