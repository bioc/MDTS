#' Circular Binary Segmentation on Minimum Distances
#'
#' This function will return a GRanges object containing the copy number 
#'	segments of all families in the input minimum distance matrix. It calls
#'	segment() from DNACopy with (alpha=0.001, undo.splits="sdundo", undo.SD=4)
#' @param md The minimum distance matrix produced by calcMD.
#' @param bins The set of bins determined by calcBins.
#' @keywords segmentMD
#' @export
segmentMD = function(md, bins){
	cna = CNA(genomdat=md, chrom=as.vector(seqnames(bins)),
		maploc=start(bins), data.type="logratio", sampleid=colnames(md), presorted=T)
	cbs = segment(cna, alpha=0.001, undo.splits="sdundo", undo.SD=4)
	gr = GRanges(seqnames=cbs$output$chrom, 
		IRanges(cbs$output$loc.start, cbs$output$loc.end), 
		m=cbs$output$seg.mean, 
		famid=cbs$output$ID)
	return(gr)
}