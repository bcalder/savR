% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/savR-methods.R
\docType{methods}
\name{savR}
\alias{savR}
\alias{savR,character-method}
\alias{savR,missing-method}
\title{Build a SAV project}
\usage{
savR(object)

\S4method{savR}{character}(object)

\S4method{savR}{missing}()
}
\arguments{
\item{object}{String Path to Flowcell data}
}
\description{
Constructor to build a \link{savProject-class} object and populate it. A SAV 
project consists of binary files generated by an Illumina sequencing run
and placed in a folder named "InterOp". This folder contains a number
of ".bin" files that contain statistics about the run.  Creating
this object parses all of the files and makes the data available for analysis.
}
\examples{
fc <- savR(system.file("extdata", "MiSeq", package="savR"))
fc

}

