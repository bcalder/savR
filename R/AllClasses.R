#'Illumina read
#'
#'Class representation of the features of an Illumina sequencing read.
#' 
#'@section Slots:
#'\describe{
#'\item{\code{number}:}{the index of this read in sequencing}
#'\item{\code{cycles}:}{number of cycles in this read}
#'\item{\code{index}:}{logical representing whether or not this read is an index read}
#'}
setClass("illuminaRead",   slots=c(number="integer", 
                                   cycles="integer",
                                   index="logical"))

#'Layout of an Illumina flowcell
#'
#'Class representation of the features of an Illumina flow cell.
#' 
#'@section Slots:
#'\describe{
#'\item{\code{lanecount}:}{Number of lanes on the flowcell}
#'\item{\code{surfacecount}:}{Number of surfaces}
#'\item{\code{swathcount}:}{Number of imaging swaths}
#'\item{\code{tilecount}:}{Number of tiles per swath}
#'}
setClass("illuminaFlowCellLayout", slots=c(lanecount="integer", 
                                           surfacecount="integer", 
                                           swathcount="integer", 
                                           tilecount="integer"))

#'SAV project class
#' 
#'Represents a flowcell, metadata and parsed SAV information
#' 
#'@section Slots:
#'\describe{
#'\item{\code{location}:}{Full path to flowcell directory}
#'\item{\code{reads}:}{List of \link{illuminaRead-class}}
#'\item{\code{layout}:}{\link{illuminaFlowCellLayout-class}}
#'\item{\code{runid}:}{Run ID}
#'\item{\code{number}:}{Run number}
#'\item{\code{flowcell}:}{Flowcell ID}
#'\item{\code{instrument}:}{Instrument ID}
#'\item{\code{date}:}{Run date}
#'\item{\code{cycles}:}{Total number of cycles}
#'\item{\code{directions}:}{Total number of sequence runs (ends)}
#'\item{\code{parsedData}:}{SAV data}
#'}
setClass("savProject", 
         slots=c(location="character",	
                 reads="list", 
                 layout="illuminaFlowCellLayout", 
                 runid="character", 
                 number="integer", 
                 flowcell="character", 
                 instrument="character", 
                 date="character", 
                 cycles="integer", 
                 directions="integer", 
                 parsedData="list"), 
         prototype=prototype(location="."))

#'Base class for formatters
#'
#'Defines the necessary slots to create parse different binary files using
#'the same generic parser.
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savFormat", slots=c(filename="character", 
                              name="character", 
                              type="character", 
                              lengths="integer", 
                              order="character", 
                              version="integer"))

#'Corrected Intensity formatter
#'
#'Lane, tile, cycle, average intensity, corrected intensities (ACGT),
#'average corrected called intensities (ACGT), number of no-calls,
#'number of (ACGT) calls, and signal to noise ratio.
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savCorrectedIntensityFormat", contains="savFormat", 
         prototype=prototype(filename="CorrectedIntMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", "avg_intensity", paste("avg_cor", c("A", "C", "G", "T"), sep="_"), 
                                    paste("avg_cor_called", c("A", "C", "G", "T"), sep="_"),
                                    paste("num", c("none", "A", "C", "G", "T"), sep="_"), 
                                    "sig_noise"),
                             type=c(rep("integer", 17), "numeric"),
                             lengths=c(rep(2L,12), rep(4L, 6)),
                             order=c("lane", "cycle", "tile"),
                             version=2L))

#'Quality Metrics formatter
#'
#'Lane, tile, cycle, Q1-Q50 counts
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savQualityFormat", contains="savFormat", 
         prototype=prototype(filename="QMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", paste("Q", 1:50, sep="")),
                             type=c(rep("integer", 53)),
                             lengths=c(rep(2L, 3), rep(4L, 50) ),
                             order=c("lane", "cycle", "tile"),
                             version=4L))

#'Tile Metrics formatter
#'
#'Lane, tile, code, value.  Codes are:
#'
#'\tabular{ll}{
#'100 \tab Cluster Density \cr
#'101 \tab PF Cluster Density \cr
#'102 \tab Number of clusters \cr
#'103 \tab Number of PF clusters \cr
#'400 \tab Control lane \cr
#'}
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savTileFormat", contains="savFormat", 
         prototype=prototype(filename="TileMetricsOut.bin", 
                             name=c("lane", "tile", "code", "value"),
                             type=c(rep("integer", 3), "numeric"),
                             lengths=c(rep(2L, 3), 4L),
                             order=c("lane", "code", "tile"),
                             version=2L))

#'Extraction Metrics formatter
#'
#'Lane, tile, cycle, FWHM (ACGT), intensity (ACGT), datestamp, timestamp.
#'Datestamp and timestamp are munged at the moment because R does not 
#'have native support for 32-bit unsigned integers and I have not implemented 
#'a solution.
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savExtractionFormat", contains="savFormat", 
         prototype=prototype(filename="ExtractionMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", paste("FWHM", c("A", "C", "G", "T"), sep="_"), paste("int", c("A", "C", "G", "T"), sep="_"), "datestamp", "timestamp"),
                             type=c(rep("integer", 3), rep("numeric", 4), rep("integer", 6)),
                             lengths=c(rep(2L, 3), rep(4L,4), rep(2L,4), rep(4L,2) ),
                             order=c("lane", "cycle", "tile"),
                             version=2L))

setClass("savParser", slots=c(project="savProject", format="savFormat"))

