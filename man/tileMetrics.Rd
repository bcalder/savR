\docType{methods}
\name{tileMetrics}
\alias{tileMetrics}
\alias{tileMetrics,savProject-method}
\title{Get Tile Metrics}
\usage{
tileMetrics(project)

\S4method{tileMetrics}{savProject}(project)
}
\arguments{
  \item{project}{SAV project}
}
\value{
sorted data.frame of tile metrics
}
\description{
Returns the Tile Metrics SAV data.
}
\details{
Metrics for each tile are encoded in the following format:
\tabular{ll}{ cluster density: \tab 100 \cr PF cluster
density: \tab 101 \cr number of clusters: \tab 102 \cr
number of PF clusters: \tab 103 \cr phasing for read N:
\tab (200 + (N - 1) * 2) \cr prephasing for read N: \tab
(201 + (N - 1) * 2) \cr percent aligned for read N: \tab
(300 + N - 1) \cr control lane: \tab 400 \cr }

\describe{ \item{\code{lane}:}{Lane number}
\item{\code{tile}:}{Tile ID} \item{\code{code}:}{Code
described above} \item{\code{value}:}{Value for code key} }
}
\examples{
example(savR)
colnames(tileMetrics(fc))
}
\references{
Codes for Tile Metrics were obtained from the Python
Illuminate package: \cr
\url{https://bitbucket.org/invitae/illuminate}
}

