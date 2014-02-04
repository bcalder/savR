\docType{methods}
\name{qualityHeatmap}
\alias{qualityHeatmap}
\alias{qualityHeatmap,savProject,integer,integer-method}
\alias{qualityHeatmap,savProject,numeric,numeric-method}
\title{Generate a heatmap of qualities}
\usage{
qualityHeatmap(project, lane, read)

\S4method{qualityHeatmap}{savProject,integer,integer}(project, lane, read)

\S4method{qualityHeatmap}{savProject,numeric,numeric}(project, lane, read)
}
\arguments{
  \item{project}{SAV project}

  \item{lane}{integer lane specification}

  \item{read}{integer vector of sequence reads to include
  (not including index reads)}
}
\description{
Plots a heatmap of quality vs cycle for a given lane for 1
or more sequence reads.  Read qualities include sequence +
index.
}

