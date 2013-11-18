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
setClass("illuminaRead", 	slots=c(number="integer", cycles="integer",
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
      surfacecount="integer", swathcount="integer", tilecount="integer"))

#'SAV project class
#' 
#'Represents a flowcell, metadata and parsed SAV information
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

#'Build a SAV project
#' 
#'Constructor to build a \link{savProject} object and populate it.
#' 
#'@param object Path to Flowcell data
#'@export
#'@examples
#'fc <- savR(system.file("extdata", "MiSeq", package="savR"))
#'fc
#'@rdname savR
setGeneric("savR", function(object) standardGeneric("savR"))

#'@rdname savR
#'@aliases savR,character-method
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

#'@rdname savR
#'@aliases savR,missing-method
setMethod("savR", signature("missing"), function() { savR(".") })

#'Get Flowcell folder location
#'
#'@param project SAV project
#'@return path to data
#'@export
#'@rdname location
setGeneric("location", function(project) standardGeneric("location"))

#'@rdname location
#'@aliases location,savProject-method
setMethod("location", signature(project="savProject"), function(project) project@location)

#'Get reads
#'
#'@param project SAV project
#'@return List of \link{illuminaRead}
#'@export
#'@rdname reads
setGeneric("reads", function(project) standardGeneric("reads"))

#'@rdname reads
#'@aliases reads,savProject-method
setMethod("reads", signature(project="savProject"), function(project) project@reads)

#'Get flowcell layout
#'
#'@param project SAV project
#'@return \link{illuminaFlowCellLayout}
#'@export
#'@rdname flowcellLayout
setGeneric("flowcellLayout", function(project) standardGeneric("flowcellLayout"))

#'@rdname flowcellLayout
#'@aliases flowcellLayout,savProject-method
setMethod("flowcellLayout", signature(project="savProject"), function(project) project@layout)

#'Get the Run ID
#'
#'@param project SAV project
#'@return parsed run id
#'@export
#'@rdname run
setGeneric("run", function(project) standardGeneric("run"))

#'@rdname run
#'@aliases run,savProject-method
setMethod("run", signature(project="savProject"), function(project) project@runid)

#'Get the total number of cycles
#'
#'@param project SAV project
#'@return number of cycles in run
#'@export
#'@rdname cycles
setGeneric("cycles", function(project) standardGeneric("cycles"))

#'@rdname cycles
#'@aliases cycles,savProject-method
setMethod("cycles", signature(project="savProject"), function(project) project@cycles)

#'Get the number of sequence reads
#'
#'Returns the number of sequence reads (excluding index reads)
#'@param project SAV project
#'@return number of reads
#'@export
#'@rdname directions
setGeneric("directions", function(project) standardGeneric("directions"))

#'@rdname directions
#'@aliases directions,savProject-method
setMethod("directions", signature(project="savProject"), function(project) project@directions)

#'Get Corrected Intensity data
#'
#'@param project SAV project
#'@return sorted data.frame of CI data
#'@export
#'@rdname correctedIntensities
setGeneric("correctedIntensities", function(project) standardGeneric("correctedIntensities"))

#'@rdname correctedIntensities
#'@aliases correctedIntensities,savProject-method
setMethod("correctedIntensities", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savCorrectedIntensityFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })

#'get Quality Metrics
#'
#'@param project SAV project
#'@return sorted data.frame quality data
#'@export
#'@rdname qualityMetrics
setGeneric("qualityMetrics", function(project) standardGeneric("qualityMetrics"))

#'@rdname qualityMetrics
#'@aliases qualityMetrics,savProject-method
setMethod("qualityMetrics", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savQualityFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })

#'Get Tile Metrics
#'
#'@param project SAV project
#'@return sorted data.frame of tile metrics
#'@export
#'@rdname tileMetrics
setGeneric("tileMetrics", function(project) standardGeneric("tileMetrics"))

#'@rdname tileMetrics
#'@aliases tileMetrics,savProject-method
setMethod("tileMetrics", signature(project="savProject"), function(project) project@parsedData[["savTileFormat"]])

#'Get Extraction Metrics
#'
#'@param project SAV project
#'@return sorted data.frame of Extraction metrics
#'@export
#'@rdname extractionMetrics
setGeneric("extractionMetrics", function(project) standardGeneric("extractionMetrics"))

#'@rdname extractionMetrics
#'@aliases extractionMetrics,savProject-method
setMethod("extractionMetrics", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savExtractionFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })
