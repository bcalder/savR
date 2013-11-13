\name{savR}
\alias{savR}
\title{Build a SAV project}
\arguments{
  \item{object}{Path to Flowcell data}
}
\description{
  Method to build a \link{savProject} object and populate
  it
}
\examples{
fc <- savR(system.file("extdata", "HiSeq", package="savR"))
fc
}

