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
#'\item{\code{sectionperlane}:}{Number of sections per lane (NextSeq)}
#'\item{\code{lanepersection}:}{Number of lanes per section (NextSeq)}
#'\item{\code{tilenamingconvention}:}{Description of deviation from original formatting layout}
#'}
setClass("illuminaFlowCellLayout", slots=c(lanecount="integer", 
                                           surfacecount="integer", 
                                           swathcount="integer", 
                                           tilecount="integer",
                                           sectionperlane="integer",
                                           lanepersection="integer",
                                           tilenamingconvention="character"
                                           ))

#'Structure for holding parsed InterOp headers and data
#'
#'@section Slots:
#'\describe{
#'\item{\code{header}:}{list of parsed header values}
#'\item{\code{data}:}{data.frame of parsed values}
#'}
setClass("savData", slots=c(header="list", data="data.frame", accessor="character"),
         prototype=prototype(header=list(), data=NULL, accessor=NULL))

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
#'\item{\code{default}:}{logical default format ()}
#'}
setClass("savFormat", slots=c(filename="character", 
                              name="character", 
                              type="character", 
                              lengths="integer", 
                              order="character", 
                              version="integer",
                              accessor="character",
                              default="logical"))

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
                             version=2L,
                             accessor="correctedIntensities",
                             default=TRUE))

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
                             version=4L,
                             accessor="qualityMetrics",
                             default=TRUE))


#'Quality Metrics formatter version 5
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
#'
# Format information found at https://tracker.tgac.ac.uk/browse/MISO-138
# Quality Metrics (QMetricsOut.bin)
# Format:
# byte 0: file version number (5)
# byte 1: length of each record
# byte 2: quality score binning (byte flag representing if binning was on), if (byte 2 == 1) // quality score binning on
# byte 3: number of quality score bins, B
# // if byte 2 == 1
#   bytes 4 - (4+B-1): lower boundary of quality score bins
#   bytes (4+B) - (4+2*B-1): upper boundary of quality score bins
#   bytes (4+2*B) - (4+3*B-1): remapped scores of quality score bins
# The remaining bytes are for the records, with each record in this format:
# 2 bytes: lane number  (uint16)
# 2 bytes: tile number  (uint16)
# 2 bytes: cycle number (uint16)
# 4 x 50 bytes: number of clusters assigned score (uint32) Q1 through Q50
# Where N is the record index
setClass("savQualityFormatV5", contains="savFormat", 
         prototype=prototype(filename="QMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", paste("Q", 1:50, sep="") ),
                             type=c(rep("integer", 53)),
                             lengths=c(2L, 2L, 2L, rep(4L, 50)),
                             order=c("lane", "cycle", "tile"),
                             accessor="qualityMetrics",
                             version=5L,
                             default=FALSE))

#'Quality Metrics formatter version 6
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
#'
# Format information found at https://tracker.tgac.ac.uk/browse/MISO-138
# Quality Metrics (QMetricsOut.bin)
# Format:
# byte 0: file version number (6)
# byte 1: length of each record
# byte 2: quality score binning (byte flag representing if binning was on), if (byte 2 == 1) // quality score binning on
# byte 3: number of quality score bins, B
# // if byte 2 == 1
#   bytes 4 - (4+B-1): lower boundary of quality score bins
#   bytes (4+B) - (4+2*B-1): upper boundary of quality score bins
#   bytes (4+2*B) - (4+3*B-1): remapped scores of quality score bins
# The remaining bytes are for the records, with each record in this format:
# 2 bytes: lane number  (uint16)
# 2 bytes: tile number  (uint16)
# 2 bytes: cycle number (uint16)
# 4 x B bytes: number of clusters assigned score (uint32) Q1 through QB
# Where N is the record index
setClass("savQualityFormatV6", contains="savFormat", 
         prototype=prototype(filename="QMetricsOut.bin", 
                             name=c("lane", "tile", "cycle"),
                             type=c(rep("integer", 53)),
                             lengths=c(2L, 2L, 2L, rep(4L, 50)),
                             order=c("lane", "cycle", "tile"),
                             accessor="qualityMetrics",
                             version=6L,
                             default=FALSE))


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
#'\item{\code{version}:}{integer version number (header consists of version (1b), length (1b))}
#'}
setClass("savTileFormat", contains="savFormat", 
         prototype=prototype(filename="TileMetricsOut.bin", 
                             name=c("lane", "tile", "code", "value"),
                             type=c(rep("integer", 3), "numeric"),
                             lengths=c(rep(2L, 3), 4L),
                             order=c("lane", "code", "tile"),
                             version=2L,
                             accessor="tileMetrics",
                             default=TRUE))

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
                             name=c("lane", "tile", "cycle", 
                                    paste("FWHM", c("A", "C", "G", "T"), sep="_"), 
                                    paste("int", c("A", "C", "G", "T"), sep="_"), "datestamp", "timestamp"),
                             type=c(rep("integer", 3), rep("numeric", 4), rep("integer", 6)),
                             lengths=c(rep(2L, 3), rep(4L,4), rep(2L,4), rep(4L,2) ),
                             order=c("lane", "cycle", "tile"),
                             version=2L,
                             accessor="extractionMetrics",
                             default=TRUE))

#'Error Metrics formatter
#'
#'Lane, tile, cycle, errorrate, nPerfect, n1Error, n2Error,
#'n3Error, n4Error.
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
setClass("savErrorFormat", contains="savFormat",
         prototype=prototype(filename="ErrorMetricsOut.bin",
                             name=c("lane", "tile", "cycle", "errorrate", "nPerfect", paste("n", 1:4, "Error", sep="")),
                             type=c(rep("integer", 3), "numeric", rep("integer", 5)),
                             lengths=c(rep(2L, 3), rep(4L, 6)),
                             order=c("lane", "cycle", "tile"),
                             version=3L,
                             accessor="errorMetrics",
                             default=TRUE))
                             
#'Index Metrics formatter
#'
#'Lane, tile, read, index, cluster, sample, project
#'
#'@section Slots:
#'\describe{
#'\item{\code{name}:}{vector of column names}
#'\item{\code{type}:}{vector of data types of elements}
#'\item{\code{lengths}:}{vector of byte lengths for each element}
#'\item{\code{order}:}{vector of column names for sorting}
#'\item{\code{version}:}{integer version number}
#'}
#'
#
# Format information found at https://github.com/Illumina/interop/blob/master/src/interop/model/metrics/index_metric.cpp
# Index Metrics (IndexMetricsOut.bin)
# Format:
# byte 0: file version number (1 or 2)
#
# n records:
#
# base_read_metric
# 2 bytes: lane number  (uint16)
# 2/4 bytes: tile number  (v1: uint16, v2: uint32)
# 2 bytes: read number (uint16)
#
# index_metric
# 2 bytes: index name length (indexNameLength) (uint16)
# indexNameLength bytes: index name
# 4/8 bytes: index cluster count (v1: uint32, v2: uint64)
# 2 bytes: sample name length (sampleNameLength) (uint16)
# sampleNameLength bytes: sample name
# 2 bytes: project name length (projectNameLength) (uint16)
# projectNameLength bytes: project name
setClass("savIndexFormat", contains="savFormat",
        prototype=prototype(filename="IndexMetricsOut.bin",
                            name=c("lane", "tile", "read", "index", "cluster", "sample", "project"),
                            type=c(rep("integer", 3), "character", "integer", "character", "character"),
                            lengths=NULL,
                            order=NULL,
                            version=1L,
                            accessor="indexMetrics",
                            default=FALSE))

setClass("savParser", slots=c(project="savProject", format="savFormat"))

