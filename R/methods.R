#' @include parse.R
NULL

subsetSide <- function(data, side) {
  if (side=="top") {
    data <- data[grepl(".1..", data$tile),]
  } else if (side=="bottom") {
    data <- data[grepl(".2..", data$tile),]
  }
  return(data)
}

#' Plot flowcell intensity by base and cycle
#' 
#' Draws a representation of a flowcell, showing the scaled corrected intensity values.
#' 
#' @param project A \link{savProject} object
#' @param cycle integer cycle number
#' @param base character for nucleotide
#' @export
setGeneric("plotIntensity", function(project, cycle, base) standardGeneric("plotIntensity"))

setMethod("plotIntensity", signature(project="savProject", cycle="integer", base="character"), function(project, cycle=1L, base=c("A", "C", "G", "T")) {
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
  grid.arrange(p)
} )

setMethod("plotIntensity", signature(project="savProject", cycle="missing", base="missing"), function(project) { plotIntensity(project, 1L, "A")})
setMethod("plotIntensity", signature(project="savProject", cycle="integer", base="missing"), function(project, cycle) { plotIntensity(project, cycle, "A")})
setMethod("plotIntensity", signature(project="savProject", cycle="missing", base="character"), function(project, base) { plotIntensity(project, 1L, base)})

#'Generate FWHM plots
#'
#'Plots the average full width of clusters at half maximum (FWHM) of each tile
#'for a given cycle and base.
#'
#'@param project SAV project
#'@param cycle sequence cycle
#'@param base nucleotide base (ACGT)
#'@export
setGeneric("plotFWHM", function(project, cycle, base) standardGeneric("plotFWHM"))

setMethod("plotFWHM", signature(project="savProject", cycle="integer", base="character"), function(project, cycle=1L, base=c("A", "C", "G", "T")) {
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
  grid.arrange(p)
} )

setMethod("plotFWHM", signature(project="savProject", cycle="missing", base="missing"), function(project) { plotFWHM(project, 1L, "A")})
setMethod("plotFWHM", signature(project="savProject", cycle="integer", base="missing"), function(project, cycle) { plotFWHM(project, cycle, "A")})
setMethod("plotFWHM", signature(project="savProject", cycle="missing", base="character"), function(project, base) { plotFWHM(project, 1L, base)})

#'Get formatted data for Q GT30 plot
#'
#'@param data data.frame
#'@param cycle cycle
setGeneric("getFormatQGT30", function(data, cycle) standardGeneric("getFormatQGT30"))

setMethod("getFormatQGT30", signature(data="data.frame", cycle="integer"), function(data, cycle=1L) {
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
} )

#'Plot Quality > 30 for a flowcell
#'
#'Generate a plot for a given cycle of the percentage of clusters in each tile
#'that are >= Q30.
#'
#'@param project SAV project
#'@param cycle sequence cycle
#'@export
setGeneric("plotQGT30", function(project, cycle) standardGeneric("plotQGT30"))

setMethod("plotQGT30", signature(project="savProject", cycle="integer"), function(project, cycle=1L) {
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
  grid.arrange(p)
} )

setMethod("plotQGT30", signature(project="savProject", cycle="missing"), function(project) { plotQGT30(project, 1L)})

#'PF Boxplot
#'
#'Generate a boxplot of the numbers of clusters and the number of
#'Illumina pass-filter clusters per tile and lane
#'
#'@param project SAV project
#'@export
setGeneric("pfBoxplot", function(project) standardGeneric("pfBoxplot"))

setMethod("pfBoxplot", signature("savProject"), function(project) {
  data <- project@parsedData[["savTileFormat"]]
  if (is.null(data))
    stop("Tile data not available")
  data <- data[data$code %in% c(100,101),]
  data[data$code==100, "code"] <- "Clusters"
  data[data$code==101, "code"] <- "PF"
  p <- ggplot(data, aes(factor(lane), value)) + geom_boxplot(notch=F, aes(fill = code), alpha=.8) + ylim(0,max(data$value)) +
  theme_bw() + theme(legend.position = "bottom") + 
  labs(list(y=expression(paste("Clusters/", mm^2, sep="")), x="Lane", fill=""))
  grid.arrange(p)
} )

#'format quality data
#'
#'@param data data
#'@param lane lane
#'@param cycles cycles
#'@return formatted data
setGeneric("qFormat", function(data, lane, cycles) standardGeneric("qFormat"))

setMethod("qFormat", signature(data="data.frame", lane="integer", cycles="integer"), function(data, lane, cycles) {
  data <- data[data$lane==lane & data$cycle %in% cycles, ]
  quals <- paste("Q", 1:50, sep="")
  mat <- melt(data[,c("cycle",quals)], id=c("cycle"), measured=quals)
  mat <- dcast(mat, cycle ~ variable, sum)
  mat <- melt(mat, id=c("cycle"), measured=quals)
  mat[,2] <- as.numeric(gsub("Q", "", mat[,2]))
  colnames(mat) <- c("x", "y", "z")
  return(mat)
} )

#'read number to vector of cycle numbers
#'
#'@param project SAV project
#'@param read read number
#'@return vector of cycle numbers
setGeneric("readToCycles", function(project, read) standardGeneric("readToCycles"))

setMethod("readToCycles", signature(project="savProject", read="integer"), function(project, read) {
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
} )

#'Generate a heatmap of qualities
#'
#'@param project SAV project
#'@param lane lane
#'@param read vector of sequence reads to include (not including index reads)
#'@export
setGeneric("qualityHeatmap", function(project, lane, read) standardGeneric("qualityHeatmap") )

setMethod("qualityHeatmap", signature(project="savProject", lane="integer", read="integer"), function(project, lane, read) {
  plots <- list()
  for (x in 1:length(read)) {
    mat <- qFormat(data=project@parsedData[["savQualityFormat"]], lane=lane, cycles=readToCycles(project, read))
    plots[[x]] <- ggplot(mat, aes(x=x, y=y, z=z)) + 
      stat_contour(bins=50, geom="polygon", aes(fill=..level..)) + ylim(0,50) + 
      theme_bw() + scale_fill_gradient2(low="white", mid=muted("green"), high="red", midpoint=quantile(mat$z, .99) ) + 
      xlab("cycle") + ylab("Q")
  }
  do.call(grid.arrange, c(plots, ncol=length(plots)))
} )

#'Generate Illumina reports folder
#'
#'Generate a folder of images that approximates the format of the folder that 
#'was superceded by InterOp.
#'
#'@param project SAV project
#'@param path location to save reports folder
#'@export
setGeneric("buildReports", function(project, path) standardGeneric("buildReports"))

setMethod("buildReports", signature(project="savProject", path="character"), function(project, path=".") {
  path <- normalizePath(path)
  if (!file.exists(path))
    stop(cat("Path", path, "does not exist."))
  reports <- paste(path, "reports", sep="/")
  if (file.exists(reports))
    stop(cat("reports folder", reports, "already exists."))
  for (f in c("ByCycle", "ErrorRate", "FWHM", "Intensity", "NumGT30")) {
    assign(f, paste(reports, f, sep="/"))
    dir.create(get(f), showWarnings=F, recursive=T)
  }
  # PF plot
  Cairo(file=paste(reports, "/NumClusters By Lane.png", sep=""), width=800, height=400, dpi=72, type="png", bg="white")
  pfBoxplot(project)
  dev.off()
  # intensity plots
  path <- normalizePath(paste(reports, "Intensity", sep="/"))
  for (cycle in 1:project@cycles) {
    for (base in c("A", "C", "G", "T")) {
      Cairo(file=paste(path, "/Chart_", cycle, "_", tolower(base), ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
      plotIntensity(project, cycle, base)
      dev.off()
    }
  }
  # Q>30 plots
  path <- normalizePath(paste(reports, "NumGT30", sep="/"))
  for (cycle in 1:project@cycles) {
      Cairo(file=paste(path, "/Chart_", cycle, ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
      plotQGT30(project, cycle)
      dev.off()
  }
  # plot lane quality
  path <- normalizePath(paste(reports, "ByCycle", sep="/"))
  for (lane in 1:project@layout@lanecount) {
    Cairo(file=paste(path, "/Qscore_L", lane, ".png", sep=""), width=800, height=400, dpi=72, type="png", bg="white")
    qualityHeatmap(project, lane, 1:project@directions)
    dev.off()
  } 
  
  # FWHM plots
  path <- normalizePath(paste(reports, "FWHM", sep="/"))
  for (cycle in 1:project@cycles) {
    for (base in c("A", "C", "G", "T")) {
      Cairo(file=paste(path, "/Chart_", cycle, "_", tolower(base), ".png", sep=""), width=300, height=800, dpi=72, type="png", bg="white")
      plotFWHM(project, cycle, base)
      dev.off()
    }
  }
  
  
  
} )

setMethod("buildReports", signature(project="savProject", path="missing"), function(project) { buildReports(project, ".")})
