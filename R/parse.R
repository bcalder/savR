#' @include savR.R
NULL

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
#'}
setClass("savFormat", slots=c(filename="character", name="character", type="character", lengths="integer", order="character"))

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
#'}
setClass("savCorrectedIntensityFormat", contains="savFormat", 
         prototype=prototype(filename="CorrectedIntMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", "avg_intensity", paste("avg_cor", c("A", "C", "G", "T"), sep="_"), 
                                    paste("avg_cor_called", c("A", "C", "G", "T"), sep="_"),
                                    paste("num", c("none", "A", "C", "G", "T"), sep="_"), 
                                    "sig_noise"),
                             type=c(rep("integer", 17), "numeric"),
                             lengths=c(rep(2L,12), rep(4L, 6)),
                             order=c("lane", "cycle", "tile")))

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
#'}
setClass("savQualityFormat", contains="savFormat", 
         prototype=prototype(filename="QMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", paste("Q", 1:50, sep="")),
                             type=c(rep("integer", 53)),
                             lengths=c(rep(2L, 3), rep(4L, 50) ),
                             order=c("lane", "cycle", "tile")))

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
#'}
setClass("savTileFormat", contains="savFormat", 
         prototype=prototype(filename="TileMetricsOut.bin", 
                             name=c("lane", "tile", "code", "value"),
                             type=c(rep("integer", 3), "numeric"),
                             lengths=c(rep(2L, 3), 4L),
                             order=c("lane", "code", "tile")))

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
#'}
setClass("savExtractionFormat", contains="savFormat", 
         prototype=prototype(filename="ExtractionMetricsOut.bin", 
                             name=c("lane", "tile", "cycle", paste("FWHM", c("A", "C", "G", "T"), sep="_"), paste("int", c("A", "C", "G", "T"), sep="_"), "datestamp", "timestamp"),
                             type=c(rep("integer", 3), rep("numeric", 4), rep("integer", 6)),
                             lengths=c(rep(2L, 3), rep(4L,4), rep(2L,4), rep(4L,2) ),
                             order=c("lane", "cycle", "tile")))

setClass("savParser", slots=c(project="savProject", format="savFormat"))

#'Generic binary parser
#'
#'@param project SAV project
#'@param format savFormat subclass to define data types
#'@return sorted data.frame of parsed data
setGeneric("parseBin", function(project, format) standardGeneric("parseBin"))

setMethod("parseBin", signature(project="savProject", "savFormat"), function(project, format) {
  path <- normalizePath(paste(project@location, "InterOp", format@filename, sep="/"))
  fh <- file(path, "rb")
  version <- readBin(fh, what="integer", endian="little", size=1, signed=F)
  reclen <- readBin(fh, what="integer", endian="little", size=1, signed=F)
  if (reclen != sum(format@lengths))
    stop(cat("file's declared record size (", reclen, ") does not equal formats declared size (", sum(format@lengths), ")"))
  readlen <- 0
  for (x in project@reads) {
    readlen <- readlen + x@cycles
  }
  proj.size <- project@layout@lanecount * project@layout@surfacecount * project@layout@swathcount * project@layout@tilecount * readlen + 1
  data <- vector("list", proj.size)
  r <- 1
  while (!isIncomplete(fh)) {
    dat <- c()
    for (i in 1:length(format@lengths)) {
      if (format@type[i] != "integer") {
        dat <- c(dat, readBin(fh, what=format@type[i], size=format@lengths[i], endian="little"))
      } else if (format@type[i] == "integer" & format@lengths[i] == 2L) {
        dat <- c(dat, readBin(fh, what=format@type[i], size=format@lengths[i], endian="little", signed=F))
      } else {
        # R does not handle 32-bit unsigned int :/
        dat <- c(dat, readBin(fh, what=format@type[i], size=format@lengths[i], endian="little"))
      }
    }
    if (length(dat)==0)
      break
    data[[r]] <- dat
    r <- r + 1
  }
  data.f <- as.data.frame(do.call("rbind", data))
  close(fh)
  colnames(data.f) <- format@name
  
  if (max(data.f[,"lane"]) != project@layout@lanecount)
    stop(cat("number of lanes in data file ( ", max(data.f[,"lane"]), ") does not equal project configuration value (" + project@layout@lanecount + ")"))

  data.f <- data.f[do.call(order, as.list(data.f[,format@order])),]
  return(data.f)
} )

validParser <- function(object) {
  if (length(object@format@name) != length(object@format@type) & length(object@format@type) != length(object@format@size))
    return("length of format parameters are not equal.")
  TRUE
}

setValidity("savParser", validParser)

#'Get basic flowcell statistics
#'
#'used to get flowcell information when data object has
#'lane, cycle, and tile data.
#'
#'@param data.frame of parsed data
#'@return list of statistics
setGeneric("getFlowcellStats", function(object) standardGeneric("getFlowcellStats"))

setMethod("getFlowcellStats", signature("data.frame"), function(object) {
  retval <- list()
  retval$sides  <- as.numeric(substring(object$tile,1,1))
  retval$swaths <- as.numeric(substring(object$tile,2,2))
  retval$nsides <- as.numeric(length(unique(substring(object$tile,1,1))))
  retval$nswath <- as.numeric(length(unique(substring(object$tile,2,2))))
  retval$ntiles <- as.numeric(substr(max(object$tile),3,4))
  retval$ncycle <- max(object$cycle)
  retval$nlanes <- max(object$lane)
  return(retval)
} ) 

#'Add position data to parsed data
#'
#'Adds and x and a y column to parsed data.  These are used for
#'laying out tiles in a tile plot.  Values are organized by
#'lane, then by swath and surface.
#'
#'@param data data.frame of parsed data
#'@return annotated data.frame
setGeneric("addPosition", function(data) standardGeneric("addPosition"))

setMethod("addPosition", signature(data="data.frame"), function(data) { # add position data to a data frame
	##< addPosition
	### This is an internal method for annotating flowcell data with XY coordinates
	### used in tile plots of flowcell lanes.
  stats <- getFlowcellStats(data)
  return(cbind(data, 
               x=((data$lane-1)*(stats$nswath*stats$nside)+1)+(stats$swaths-1)+((stats$sides-1)*stats$nsides + (stats$sides-1)), 
               y=rep(rep(1:stats$ntiles, stats$nswath*stats$nsides*stats$ncycle), stats$nlanes)))
} )

#'Do parsing
#'
#'After everything is configured, initialize parsing of SAV files.
#'
#'@param project SAV project
setGeneric("init", function(project) standardGeneric("init"))

setMethod("init", signature("savProject"), function(project) {
  validFormats <- c("savCorrectedIntensityFormat", "savQualityFormat", "savTileFormat", "savExtractionFormat")
  
  for (x in validFormats) {
    format <- new(x)
    if (file.exists(normalizePath(paste(project@location, "InterOp", format@filename, sep="/")))) {
      data <- parseBin(project, format)
      # don't add position data to tiles
      if (class(format)[1] != "savTileFormat")
        data <- addPosition(data)
      # removed unparsed date columns
      if (class(format)[1] == "savExtractionFormat")
        data <- data[,-c(12:13)]
      project@parsedData[[x]] <- data
    }
  }
  return(project)
})
  
