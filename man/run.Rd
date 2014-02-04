\docType{methods}
\name{run}
\alias{run}
\alias{run,savProject-method}
\title{Get the Run ID}
\usage{
run(project)

\S4method{run}{savProject}(project)
}
\arguments{
  \item{project}{SAV project}
}
\value{
parsed Illumina run id
}
\description{
Accessor to obtain the string identifier of an Illumina
sequencing run.
}
\examples{
example(savR)
run(fc)
}

