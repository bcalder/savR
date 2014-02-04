\docType{methods}
\name{plotQGT30}
\alias{plotQGT30}
\alias{plotQGT30,savProject,integer-method}
\alias{plotQGT30,savProject,missing-method}
\title{Plot Quality > 30 for a flowcell}
\usage{
plotQGT30(project, cycle)

\S4method{plotQGT30}{savProject,integer}(project, cycle = 1L)

\S4method{plotQGT30}{savProject,missing}(project)
}
\arguments{
  \item{project}{SAV project}

  \item{cycle}{sequence cycle}
}
\description{
Generate a plot for a given cycle of the percentage of
clusters in each tile that are >= Q30.
}

