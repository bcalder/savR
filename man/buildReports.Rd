\docType{methods}
\name{buildReports}
\alias{buildReports}
\alias{buildReports,savProject,character-method}
\alias{buildReports,savProject,missing-method}
\title{Generate Illumina reports folder}
\usage{
buildReports(project, destination)

\S4method{buildReports}{savProject,character}(project,
  destination = "./savR-reports")

\S4method{buildReports}{savProject,missing}(project)
}
\arguments{
  \item{project}{SAV project}

  \item{destination}{path to save reports folder}
}
\description{
Generate a folder of images that approximates the format of
the folder that was superceded by InterOp. Requires the
Cairo package.
}
\examples{
\dontrun{
example(savR)
buildReports(fc, "reports")
}
}

