#' Visualization for deletions
#'
#' This function plots the raw read information from the location of interest for a family.
#' @param deletion A GRanges object in the format of the output of denovoDeletions().
#' @param bins The set of bins determined by calcBins().
#' @param pData A table in the format of the output of pData().
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param save If TRUE will save plot to current working directory instead of rendering.
#' @keywords visualizeDeletion
#' @export
visualizeDeletion = function(deletion, bins, pD, mCounts, save=F){
	window = 1000
	famid = pD$family_id[pD$subj_id==deletion$famid]
	pD_sub = pD[str_detect(pD$family_id, famid),]
	track = deletion
	bait = GRanges(seqnames(track)[1], IRanges(min(start(track))-window, max(end(track))+window))
	row_inds = subjectHits(findOverlaps(track, bins))
	col_inds = match(pD_sub$subj_id, colnames(mCounts))
	
	id = paste0(famid, "-", seqnames(deletion), ":", start(deletion), "-", end(deletion), ".pdf")
	if(save==T){
		pdf(width=30, height=8, file=id)	
	}else{
		dev.new(width=30, height=8)	
	}
	visualizeFamily(famid, bait, track, window, row_inds, col_inds, pD)
	dev.off()
}

