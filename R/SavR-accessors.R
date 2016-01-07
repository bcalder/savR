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
setMethod("correctedIntensities", signature(project="savProject"), function(project) { 
  tmp <- project@parsedData[["savCorrectedIntensityFormat"]]@data
  if (is.null(tmp)) return(tmp)
  return(tmp[,!colnames(tmp) %in% c("x", "y")]) 
})

#'@rdname qualityMetrics
#'@aliases qualityMetrics,savProject-method
setMethod("qualityMetrics", signature(project="savProject"), function(project) {
	tmp <- NULL
	if("savQualityFormat" %in% names(project@parsedData)){
		tmp <- project@parsedData[["savQualityFormat"]]@data
	}else if("savQualityFormatV5" %in% names(project@parsedData)){
		tmp <- project@parsedData[["savQualityFormatV5"]]@data
	}
	if (is.null(tmp)) return(tmp)
	return(tmp[,!colnames(tmp) %in% c("x", "y")])
})

#'@rdname tileMetrics
#'@aliases tileMetrics,savProject-method
setMethod("tileMetrics", signature(project="savProject"), function(project) {
  tmp <- project@parsedData[["savTileFormat"]]@data
  return(tmp)
})

#'@rdname extractionMetrics
#'@aliases extractionMetrics,savProject-method
setMethod("extractionMetrics", signature(project="savProject"), function(project) { 
  tmp <- project@parsedData[["savExtractionFormat"]]@data
  if (is.null(tmp)) return(tmp)
  return(tmp[,!colnames(tmp) %in% c("x", "y")]) 
})

#'@rdname errorMetrics
#'@aliases errorMetrics,savProject-method
setMethod("errorMetrics", signature(project="savProject"), function(project) {
  tmp <- project@parsedData[["savErrorFormat"]]@data
  if (is.null(tmp)) return(tmp)
  return(tmp[,!colnames(tmp) %in% c("x", "y")])
})

#'@rdname clusters
#@aliases clusters,savProject,integer
setMethod("clusters", signature(project="savProject", lane="integer"), function(project, lane=1L) {
  if (!all(lane %in% 1:flowcellLayout(project)@lanecount)) {
    stop(paste("lane" , lane, "is not consistent with number of lanes on flowcell (", flowcellLayout(project)@lanecount, ")", sep=" "))
  }
  tm <- tileMetrics(project)
  return(sum(tm[tm$lane %in% lane & tm$code==102,]$value))
})

#'@rdname pfClusters
#'@aliases pfClusters,savProject,integer
setMethod("pfClusters", signature(project="savProject", lane="integer"), function(project, lane=1L) {
  if (!all(lane %in% 1:flowcellLayout(project)@lanecount)) {
    stop(paste("lane" , lane, "is not consistent with number of lanes on flowcell (", flowcellLayout(project)@lanecount, ")", sep=" "))
  }
  tm <- tileMetrics(project)
  return(sum(tm[tm$lane %in% lane & tm$code==103,]$value))
})
