% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/SavR-accessors.R
\docType{methods}
\name{location}
\alias{location}
\alias{location,savProject-method}
\title{Get Flowcell folder location}
\usage{
location(project)

\S4method{location}{savProject}(project)
}
\arguments{
\item{project}{SAV project}
}
\value{
normalized path to Illumina run data.
}
\description{
Accessor to obtain the path to data for a particular SAV project.
}
\examples{
example(savR)
location(fc)
}

