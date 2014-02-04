#'@rdname location
#'@aliases location,savProject-method
setMethod("location", signature(project="savProject"), function(project) project@location)

#'@rdname reads
#'@aliases reads,savProject-method
setMethod("reads", signature(project="savProject"), function(project) project@reads)

#'@rdname flowcellLayout
#'@aliases flowcellLayout,savProject-method
setMethod("flowcellLayout", signature(project="savProject"), function(project) project@layout)

#'@rdname run
#'@aliases run,savProject-method
setMethod("run", signature(project="savProject"), function(project) project@runid)

#'@rdname cycles
#'@aliases cycles,savProject-method
setMethod("cycles", signature(project="savProject"), function(project) project@cycles)

#'@rdname directions
#'@aliases directions,savProject-method
setMethod("directions", signature(project="savProject"), function(project) project@directions)

#'@rdname correctedIntensities
#'@aliases correctedIntensities,savProject-method
setMethod("correctedIntensities", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savCorrectedIntensityFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })

#'@rdname qualityMetrics
#'@aliases qualityMetrics,savProject-method
setMethod("qualityMetrics", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savQualityFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })

#'@rdname tileMetrics
#'@aliases tileMetrics,savProject-method
setMethod("tileMetrics", signature(project="savProject"), function(project) project@parsedData[["savTileFormat"]])

#'@rdname extractionMetrics
#'@aliases extractionMetrics,savProject-method
setMethod("extractionMetrics", signature(project="savProject"), function(project) { tmp <- project@parsedData[["savExtractionFormat"]]; return(tmp[,!colnames(tmp) %in% c("x", "y")]) })
