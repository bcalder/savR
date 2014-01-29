# @include savR.R
NULL

#'@rdname savR
#'@aliases savR,character-method
setMethod("savR", signature("character"), function(object) {
  retval <- new("savProject", location=normalizePath(object))
  retval@cycles <- 0L
  retval@directions <- 0L
  ri <- normalizePath(paste(object, "RunInfo.xml", sep="/"))
  runinfo <- XML::xmlInternalTreeParse(ri)
  retval@runid <- XML::xmlAttrs(XML::xpathApply(runinfo, "/RunInfo/Run")[[1]])["Id"]
  retval@number <- as.integer(XML::xmlAttrs(xpathApply(runinfo, "/RunInfo/Run")[[1]])["Number"])
  retval@flowcell <- XML::xmlValue(XML::xpathApply(runinfo, "/RunInfo/Run/Flowcell")[[1]])
  retval@instrument <- XML::xmlValue(XML::xpathApply(runinfo, "/RunInfo/Run/Instrument")[[1]])
  retval@date <- XML::xmlValue(XML::xpathApply(runinfo, "/RunInfo/Run/Date")[[1]])
  reads <- c()
  for (x in XML::xpathApply(runinfo, "/RunInfo/Run/Reads/Read")) {
    index <- XML::xmlAttrs(x)["IsIndexedRead"]
    index <- if(index=="Y") T else F
    read <- new("illuminaRead", number=as.integer(XML::xmlAttrs(x)["Number"]), 
                cycles=as.integer(XML::xmlAttrs(x)["NumCycles"]),
                index=index)
    reads <- c(reads, read)
    retval@cycles <- retval@cycles + read@cycles
    if (!read@index)
      retval@directions <- retval@directions + 1L
  } 
  retval@reads <- reads
  layout <- XML::xpathApply(runinfo, "/RunInfo/Run/FlowcellLayout")[[1]]
  retval@layout <- new("illuminaFlowCellLayout", lanecount=as.integer(XML::xmlAttrs(layout)["LaneCount"]),
                       surfacecount=as.integer(XML::xmlAttrs(layout)["SurfaceCount"]),
                       swathcount=as.integer(XML::xmlAttrs(layout)["SwathCount"]),
                       tilecount=as.integer(XML::xmlAttrs(layout)["TileCount"]) )
  return(init(retval))
} )

#'@rdname savR
#'@aliases savR,missing-method
setMethod("savR", signature("missing"), function() { savR(".") })

subsetSide <- function(data, side) {
  if (side=="top") {
    data <- data[grepl(".1..", data$tile),]
  } else if (side=="bottom") {
    data <- data[grepl(".2..", data$tile),]
  }
  return(data)
}



#'@rdname plotIntensity
#'@aliases plotIntensity,savProject,integer,character-method
setMethod("plotIntensity", signature(project="savProject", cycle="integer", base="character"), function(project, cycle=1L, base=c("A", "C", "G", "T")) {
  x <- y <- NULL
  if (cycle < 0)
    stop ("Cycle out of range")
  data <- project@parsedData[["savCorrectedIntensityFormat"]]
  if (is.null(data))
    stop("Corrected Intensity data not available")
  val <- paste("avg_cor_called", c("A", "C", "G", "T"), sep="_") 
  names(val) <- c("A", "C", "G", "T")
  maxInt <- max(c(data[, val["A"]], data[, val["C"]], data[, val["G"]], data[, val["T"]]))
  if(maxInt < 7000) {
    maxInt <= 7000
  }
  data <- data[data$cycle==cycle,]
  base <- match.arg(base)
  
  p <- qplot(factor(x),y,fill=get(val[base]), data=data, geom="tile", position="dodge", main = paste("Intensity: ", base, ", Cycle ", cycle, sep="")) + 
    theme_bw() + theme(legend.position = "bottom") + scale_fill_continuous(guide = guide_colorbar(title=val, barwidth=10), limits=c(0,maxInt) ) + 
    facet_grid(~lane, space="free", scales="free") + 
    xlab("") + ylab("") + scale_x_discrete(labels="") 
  gridExtra::grid.arrange(p)
} )

#'@rdname plotIntensity
#'@aliases plotIntensity,savProject,missing,missing-method
setMethod("plotIntensity", signature(project="savProject", cycle="missing", base="missing"), function(project) { plotIntensity(project, 1L, "A")})
#'@rdname plotIntensity
#'@aliases plotIntensity,savProject,integer,missing-method
setMethod("plotIntensity", signature(project="savProject", cycle="integer", base="missing"), function(project, cycle) { plotIntensity(project, cycle, "A")})
#'@rdname plotIntensity
#'@aliases plotIntensity,savProject,missing,character-method
setMethod("plotIntensity", signature(project="savProject", cycle="missing", base="character"), function(project, base) { plotIntensity(project, 1L, base)})



#'@rdname plotFWHM
#'@aliases plotFWHM,savProject,integer,character-method
setMethod("plotFWHM", signature(project="savProject", cycle="integer", base="character"), function(project, cycle=1L, base=c("A", "C", "G", "T")) {
  x <- y <- NULL
  if (cycle < 0)
    stop ("Cycle out of range")
  data <- project@parsedData[["savExtractionFormat"]]
  if (is.null(data))
    stop("Extraction data not available")
  data <- data[data$cycle==cycle,]
  base <- match.arg(base)
  val <- paste("FWHM", base, sep="_")
  p <- qplot(factor(x),y,fill=get(val), data=data, geom="tile", position="dodge", main = paste("Intensity: ", base, ", Cycle ", cycle, sep="")) + 
    theme_bw() + theme(legend.position = "bottom") + scale_fill_continuous(guide = guide_colorbar(title=val, barwidth=10), limits=c(0,10) ) +
    facet_grid(~lane, space="free", scales="free") +
    xlab("") + ylab("") + scale_x_discrete(labels="") 
  gridExtra::grid.arrange(p)
} )

#'@rdname plotFWHM
#'@aliases plotFWHM,savProject,missing,missing-method
setMethod("plotFWHM", signature(project="savProject", cycle="missing", base="missing"), function(project) { plotFWHM(project, 1L, "A")})
#'@rdname plotFWHM
#'@aliases plotFWHM,savProject,integer,missing-method
setMethod("plotFWHM", signature(project="savProject", cycle="integer", base="missing"), function(project, cycle) { plotFWHM(project, cycle, "A")})
#'@rdname plotFWHM
#'@aliases plotFWHM,savProject,missing,character-method
setMethod("plotFWHM", signature(project="savProject", cycle="missing", base="character"), function(project, base) { plotFWHM(project, 1L, base)})

#Get formatted data for Q GT30 plot
#
#@param data data.frame
#@param cycle cycle
getFormatQGT30 <- function(data, cycle=1L) {
  #side <- match.arg(side)
  #data <- subsetSide(data,side)
  stats <- getFlowcellStats(data)
  contenders <- as.integer(gsub("Q", "", colnames(data)[grepl("^Q", colnames(data))]))
  lt30 <- paste("Q", 1:30, sep="")
  gte30 <- paste("Q", 31:max(contenders), sep="")
  return(cbind(data[data$cycle==cycle, c("cycle", "lane", "tile")], 
               x=factor(rep(1:(stats$nswath*stats$nsides*stats$nlanes), each=stats$ntiles)), 
               y=rep(1:stats$ntiles, stats$nswath*stats$nsides),
               gte30=apply(data[data$cycle==cycle,], 1, function(x) { 
                 return(sum(x[gte30])/(sum(x[c(lt30, gte30)])+.00001)*100 )
               })
  ))
}



#'@rdname plotQGT30
#'@aliases plotQGT30,savProject,integer-method
setMethod("plotQGT30", signature(project="savProject", cycle="integer"), function(project, cycle=1L) {
  x <- y <- gte30 <- NULL
  if (cycle < 0)
    stop ("Cycle out of range")
  data <- project@parsedData[["savQualityFormat"]]
  if (is.null(data))
    stop("Quality data not available")
  cycleData <- getFormatQGT30(data, cycle)
  p <- qplot(x, y, fill=gte30, data=cycleData, geom="tile", position="dodge", main = paste("Percent Q>=30, Cycle ", cycle, sep="")) + 
    theme_bw() + theme(legend.position = "bottom") + scale_fill_continuous(guide = guide_colorbar(title="%Q>=30", barwidth=10), limits=c(0,100) ) +
    facet_grid(~lane, space="free", scales="free") +
    xlab("") + ylab("") + scale_x_discrete(labels="") 
  gridExtra::grid.arrange(p)
} )

#'@rdname plotQGT30
#'@aliases plotQGT30,savProject,missing-method
setMethod("plotQGT30", signature(project="savProject", cycle="missing"), function(project) { plotQGT30(project, 1L)})



#'@rdname pfBoxplot
#'@aliases pfBoxplot,savProject-method
setMethod("pfBoxplot", signature("savProject"), function(project) {
  lane <- value <- code <- NULL
  data <- project@parsedData[["savTileFormat"]]
  if (is.null(data))
    stop("Tile data not available")
  data <- data[data$code %in% c(100,101),]
  data[data$code==100, "code"] <- "Clusters"
  data[data$code==101, "code"] <- "PF"
  p <- ggplot2::ggplot(data, ggplot2::aes(factor(lane), value)) + ggplot2::geom_boxplot(notch=F, ggplot2::aes(fill = code), alpha=.8) + ggplot2::ylim(0,max(data$value)) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom") + 
    ggplot2::labs(list(y=expression(paste("Clusters/", mm^2, sep="")), x="Lane", fill=""))
  gridExtra::grid.arrange(p)
} )

#format quality data
#
#@param data data
#@param lane lane
#@param cycles cycles
#@return formatted data
qFormat <- function(data,lane,cycles) {
  data <- data[data$lane==lane & data$cycle %in% cycles, ]
  quals <- paste("Q", 1:50, sep="")
  mat <- reshape2::melt(data[,c("cycle",quals)], id=c("cycle"), measured=quals)
  mat <- reshape2::dcast(mat, cycle ~ variable, sum)
  mat <- reshape2::melt(mat, id=c("cycle"), measured=quals)
  mat[,2] <- as.numeric(gsub("Q", "", mat[,2]))
  colnames(mat) <- c("x", "y", "z")
  return(mat)
}

#read number to vector of cycle numbers
#
#@param project SAV project
#@param read read number
#@return vector of cycle numbers
readToCycles <- function(project, read) {
  cycles <- c()
  indexed <- c()
  for (x in project@reads) {
    cycles <- c(cycles, x@cycles)
    indexed <- c(indexed, x@index)
  }
  seqreadlen <- cycles[!indexed]
  read <- 0
  for (i in 1:length(indexed)) {
    if (!indexed[i]) {
      read <- read + 1
    } else {
      seqreadlen[read] <- seqreadlen[read] + cycles[i]
    }
  }
  start <- 1
  end <- 0
  result <- list()
  for (r in 1:length(seqreadlen)) {
    end <- end + seqreadlen[r]
    result[[r]] <- start:end
    start <- start + seqreadlen[r]
  }
  return(result[[read]])     
}



#'@rdname qualityHeatmap
#'@aliases qualityHeatmap,savProject,integer,integer-method
setMethod("qualityHeatmap", signature(project="savProject", lane="integer", read="integer"), function(project, lane, read) {
  y <- z <- ..level.. <- NULL
  plots <- list()
  if (!all( read %in% 1:directions(project)))
    stop(paste("There are only", directions(project), "sequence read(s) on this flowcell, check read specification."))
  for (x in 1:length(read)) {
    mat <- qFormat(data=project@parsedData[["savQualityFormat"]], lane=lane, cycles=readToCycles(project, read))
    plots[[x]] <- ggplot2::ggplot(mat, ggplot2::aes(x=x, y=y, z=z)) + 
      ggplot2::stat_contour(bins=50, geom="polygon", ggplot2::aes(fill=..level..)) + ggplot2::ylim(0,50) + 
      ggplot2::theme_bw() + ggplot2::scale_fill_gradient2(low="white", mid=scales::muted("green"), high="red", midpoint=quantile(mat$z, .99) ) + 
      xlab("cycle") + ylab("Q")
  }
  do.call(gridExtra::grid.arrange, c(plots, ncol=length(plots)))
} )

#'@rdname qualityHeatmap
#'@aliases qualityHeatmap,savProject,numeric,numeric-method
setMethod("qualityHeatmap", signature(project="savProject", lane="numeric", read="numeric"), function(project, lane, read) { qualityHeatmap(project, as.integer(lane), as.integer(read))})



#'@rdname buildReports
#'@aliases buildReports,savProject,character-method
setMethod("buildReports", signature(project="savProject", destination="character"), function(project, destination=NULL) {
  path <- location(project)
  if (!file.exists(path))
    stop(paste("Project", path, "does not exist."))
  reports <- normalizePath(destination, mustWork=F)
  if (file.exists(reports))
    stop(paste("Reports folder", reports, "already exists."))
  for (f in c("ByCycle", "ErrorRate", "FWHM", "Intensity", "NumGT30")) {
    assign(f, paste(reports, f, sep="/"))
    dir.create(get(f), showWarnings=F, recursive=T)
  }
  # PF plot
  Cairo::Cairo(file=paste(reports, "/NumClusters By Lane.png", sep=""), width=800, height=400, dpi=72, type="png", bg="white")
  pfBoxplot(project)
  dev.off()
  # intensity plots
  path <- normalizePath(paste(reports, "Intensity", sep="/"))
  for (cycle in 1:project@cycles) {
    for (base in c("A", "C", "G", "T")) {
      Cairo::Cairo(file=paste(path, "/Chart_", cycle, "_", tolower(base), ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
      plotIntensity(project, cycle, base)
      dev.off()
    }
  }
  # Q>30 plots
  path <- normalizePath(paste(reports, "NumGT30", sep="/"))
  for (cycle in 1:project@cycles) {
    Cairo::Cairo(file=paste(path, "/Chart_", cycle, ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
    plotQGT30(project, cycle)
    dev.off()
  }
  # plot lane quality
  path <- normalizePath(paste(reports, "ByCycle", sep="/"))
  for (lane in 1:project@layout@lanecount) {
    Cairo::Cairo(file=paste(path, "/Qscore_L", lane, ".png", sep=""), width=800, height=400, dpi=72, type="png", bg="white")
    qualityHeatmap(project, lane, 1:project@directions)
    dev.off()
  } 
  
  # FWHM plots
  path <- normalizePath(paste(reports, "FWHM", sep="/"))
  for (cycle in 1:project@cycles) {
    for (base in c("A", "C", "G", "T")) {
      Cairo::Cairo(file=paste(path, "/Chart_", cycle, "_", tolower(base), ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
      plotFWHM(project, cycle, base)
      dev.off()
    }
  }  
  
} )

#'@rdname buildReports
#'@aliases buildReports,savProject,missing-method
setMethod("buildReports", signature(project="savProject", destination="missing"), function(project) { buildReports(project, "reports")})


#Generic binary parser
#
#@param project SAV project
#@param format savFormat subclass to define data types
#@return sorted data.frame of parsed data)
parseBin <- function(project, format) {
  path <- normalizePath(paste(project@location, "InterOp", format@filename, sep="/"))
  fh <- file(path, "rb")
  vers <- readBin(fh, what="integer", endian="little", size=1, signed=F)
  if (vers != format@version) {
    # TODO: check for other parsers
    close(fh)
    stop(paste("savR currently only supports version", format@version, "of this SAV file.", format@filename, "is reported as version", vers, "."))
  }
  reclen <- readBin(fh, what="integer", endian="little", size=1, signed=F)
  if (reclen != sum(format@lengths))
    stop(paste("file's declared record size (", reclen, ") does not equal formats declared size (", sum(format@lengths), ")"))
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
} 

#validParser <- function(object) {
#  if (length(object@format@name) != length(object@format@type) & length(object@format@type) != length(object@format@size))
#    return("length of format parameters are not equal.")
#  TRUE
#}

#setValidity("savParser", validParser)

#Get basic flowcell statistics
#
#used to get flowcell information when data object has
#lane, cycle, and tile data.
#
#@param data.frame of parsed data
#@return list of statistics
getFlowcellStats <- function(object) {
  retval <- list()
  retval$sides  <- as.numeric(substring(object$tile,1,1))
  retval$swaths <- as.numeric(substring(object$tile,2,2))
  retval$nsides <- as.numeric(length(unique(substring(object$tile,1,1))))
  retval$nswath <- as.numeric(length(unique(substring(object$tile,2,2))))
  retval$ntiles <- as.numeric(substr(max(object$tile),3,4))
  retval$ncycle <- max(object$cycle)
  retval$nlanes <- max(object$lane)
  return(retval)
} 

#Add position data to parsed data
#
#Adds and x and a y column to parsed data.  These are used for
#laying out tiles in a tile plot.  Values are organized by
#lane, then by swath and surface.
#
#@param data data.frame of parsed data
#@return annotated data.frame
addPosition <- function(data) {
  ##< addPosition
  ### This is an internal method for annotating flowcell data with XY coordinates
  ### used in tile plots of flowcell lanes.
  stats <- getFlowcellStats(data)
  return(cbind(data, 
               x=((data$lane-1)*(stats$nswath*stats$nside)+1)+(stats$swaths-1)+((stats$sides-1)*stats$nsides + (stats$sides-1)), 
               y=rep(rep(1:stats$ntiles, stats$nswath*stats$nsides*stats$ncycle), stats$nlanes)))
} 

#Do parsing
#
#After everything is configured, initialize parsing of SAV files.
#
#@param project SAV project
init <- function(project) {
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
} 



