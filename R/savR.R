#' Illumina read
#' 
#'@section Slots:
#'\describe{
#'\item{\code{number}:}{the index of this read in sequencing}
#'\item{\code{cycles}:}{number of cycles in this read}
#'\item{\code{index}:}{logical representing whether or not this read is an index read}
#'}
setClass("illuminaRead", 	slots=c(number="integer", cycles="integer",
		index="logical"))

#' Layout of an Illumina flowcell
#' 
#'@section Slots:
#'\describe{
#'\item{\code{lanecount}:}{Number of lanes on the flowcell}
#'\item{\code{surfacecount}:}{Number of surfaces}
#'\item{\code{swathcount}:}{Number of imaging swaths}
#'\item{\code{tilecount}:}{Number of tiles per swath}
#'}
setClass("illuminaFlowCellLayout", slots=c(lanecount="integer", 
      surfacecount="integer", swathcount="integer", tilecount="integer"))

#' SAV project class
#' 
#' Represents a flowcell, metadata and parsed SAV information
#' 
#'@section Slots:
#'\describe{
#'\item{\code{location}:}{Full path to flowcell directory}
#'\item{\code{reads}:}{List of \link{illuminaRead}}
#'\item{\code{layout}:}{\link{illuminaFlowCellLayout}}
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
	slots=c(location="character",	reads="list", layout="illuminaFlowCellLayout", 
          runid="character", number="integer", flowcell="character", 
          instrument="character", date="character", cycles="integer", 
          directions="integer", parsedData="list"), 
         prototype=prototype(location="."))

setMethod("show", "savProject", function(object) cat(class(object), "instance with", object@layout@lanecount, "lanes,", object@cycles, "total cycles and", object@directions, "sequencing runs." ))

#' Build a SAV project
#' 
#' Method to build a \link{savProject} object and populate it
#' 
#' @param object Path to Flowcell data
#' @export
#' @examples
#' fc <- savR(system.file("extdata", "HiSeq", package="savR"))
#' fc
setGeneric("savR", function(object) standardGeneric("savR"))


setMethod("savR", signature("character"), function(object) {
  retval <- new("savProject", location=normalizePath(object))
  retval@cycles <- 0L
  retval@directions <- 0L
  ri <- normalizePath(paste(object, "RunInfo.xml", sep="/"))
  runinfo <- xmlInternalTreeParse(ri)
  retval@runid <- xmlAttrs(xpathApply(runinfo, "/RunInfo/Run")[[1]])["Id"]
  retval@number <- as.integer(xmlAttrs(xpathApply(runinfo, "/RunInfo/Run")[[1]])["Number"])
  retval@flowcell <- xmlValue(xpathApply(runinfo, "/RunInfo/Run/Flowcell")[[1]])
  retval@instrument <- xmlValue(xpathApply(runinfo, "/RunInfo/Run/Instrument")[[1]])
  retval@date <- xmlValue(xpathApply(runinfo, "/RunInfo/Run/Date")[[1]])
  reads <- c()
  for (x in xpathApply(runinfo, "/RunInfo/Run/Reads/Read")) {
    index <- xmlAttrs(x)["IsIndexedRead"]
    index <- if(index=="Y") T else F
    read <- new("illuminaRead", number=as.integer(xmlAttrs(x)["Number"]), 
                cycles=as.integer(xmlAttrs(x)["NumCycles"]),
                index=index)
    reads <- c(reads, read)
    retval@cycles <- retval@cycles + read@cycles
    if (!read@index)
      retval@directions <- retval@directions + 1L
  } 
  retval@reads <- reads
  layout <- xpathApply(runinfo, "/RunInfo/Run/FlowcellLayout")[[1]]
  retval@layout <- new("illuminaFlowCellLayout", lanecount=as.integer(xmlAttrs(layout)["LaneCount"]),
                       surfacecount=as.integer(xmlAttrs(layout)["SurfaceCount"]),
                       swathcount=as.integer(xmlAttrs(layout)["SwathCount"]),
                       tilecount=as.integer(xmlAttrs(layout)["TileCount"]) )
  return(init(retval))
} )

setMethod("savR", signature("missing"), function() { savR(".") })

#'Get Flowcell folder location
#'
#'@param project SAV project
#'@return path to data
#'@export
setGeneric("location", function(project) standardGeneric("location"))
setMethod("location", signature(project="savProject"), function(project) project@location)

#'Get reads
#'
#'@param project SAV project
#'@return List of \link{illuminaRead}
#'@export
setGeneric("reads", function(project) standardGeneric("reads"))
setMethod("reads", signature(project="savProject"), function(project) project@reads)

#'Get flowcell layout
#'
#'@param project SAV project
#'@return \link{illuminaFlowCellLayout}
#'@export
setGeneric("flowcellLayout", function(project) standardGeneric("flowcellLayout"))
setMethod("flowcellLayout", signature(project="savProject"), function(project) project@layout)

#'Get the Run ID
#'
#'@param project SAV project
#'@return parsed run id
#'@export
setGeneric("run", function(project) standardGeneric("run"))
setMethod("run", signature(project="savProject"), function(project) project@runid)

#'Get the total number of cycles
#'
#'@param project SAV project
#'@return number of cycles in run
#'@export
setGeneric("cycles", function(project) standardGeneric("cycles"))
setMethod("cycles", signature(project="savProject"), function(project) project@cycles)

#'Get the number of sequence reads
#'
#'Returns the number of sequence reads (excluding index reads)
#'@param project SAV project
#'@return number of reads
#'@export
setGeneric("directions", function(project) standardGeneric("directions"))
setMethod("directions", signature(project="savProject"), function(project) project@directions)

#'Get Corrected Intensity data
#'
#'@param project SAV project
#'@return sorted data.frame of CI data
#'@export
setGeneric("correctedIntensities", function(project) standardGeneric("correctedIntensities"))
setMethod("correctedIntensities", signature(project="savProject"), function(project) project@parsedData[["savCorrectedIntensityFormat"]])

#'get Quality Metrics
#'
#'@param project SAV project
#'@return sorted data.frame quality data
#'@export
setGeneric("qualityMetrics", function(project) standardGeneric("qualityMetrics"))
setMethod("qualityMetrics", signature(project="savProject"), function(project) project@parsedData[["savQualityFormat"]])

#'Get Tile Metrics
#'
#'@param project SAV project
#'@return sorted data.frame of tile metrics
#'@export
setGeneric("tileMetrics", function(project) standardGeneric("tileMetrics"))
setMethod("tileMetrics", signature(project="savProject"), function(project) project@parsedData[["savTileFormat"]])

#'Get Extraction Metrics
#'
#'@param project SAV project
#'@return sorted data.frame of Extraction metrics
#'@export
setGeneric("extractionMetrics", function(project) standardGeneric("extractionMetrics"))

setMethod("extractionMetrics", signature(project="savProject"), function(project) project@parsedData[["savExtractionFormat"]])
