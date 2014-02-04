\docType{class}
\name{savProject-class}
\alias{savProject-class}
\title{SAV project class}
\description{
Represents a flowcell, metadata and parsed SAV information
}
\section{Slots}{
  \describe{ \item{\code{location}:}{Full path to flowcell
  directory} \item{\code{reads}:}{List of
  \link{illuminaRead-class}}
  \item{\code{layout}:}{\link{illuminaFlowCellLayout-class}}
  \item{\code{runid}:}{Run ID} \item{\code{number}:}{Run
  number} \item{\code{flowcell}:}{Flowcell ID}
  \item{\code{instrument}:}{Instrument ID}
  \item{\code{date}:}{Run date} \item{\code{cycles}:}{Total
  number of cycles} \item{\code{directions}:}{Total number
  of sequence runs (ends)} \item{\code{parsedData}:}{SAV
  data} }
}

