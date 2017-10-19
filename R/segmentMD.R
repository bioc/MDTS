#' Circular Binary Segmentation on Minimum Distances
#'
#' This function will return a GRanges object containing the copy number 
#'	segments of all families in the input minimum distance matrix. It calls
#'	segment() from DNACopy with (alpha=0.001, undo.splits="sdundo", undo.SD=4)
#' @param md The minimum distance matrix produced by calcMD.
#' @param bins The set of bins determined by calcBins.
#' @param alpha Controls the alpha option in calling DNAcopy::segment()
#' @param undo.splits Controls the undo.splits option in calling DNAcopy::segment()
#' @param undo.SD Controls the undo.SD option in calling DNAcopy::segment()
#' @keywords segmentMD
#' @export
segmentMD = function(md, bins, alpha='0.01', undo.splits='sdundo', undo.SD=4){
	cna = CNA(genomdat=md, chrom=as.vector(seqnames(bins)),
		maploc=start(bins), data.type="logratio", sampleid=colnames(md), presorted=T)
	cbs = segment(cna, alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD)
	segRows = cbs$segRows
	gr = GRanges(seqnames=cbs$output$chrom, 
		IRanges(start(bins)[segRows[,1]], end(bins)[segRows[,2]]), 
		m=cbs$output$seg.mean, 
		famid=cbs$output$ID)
	return(gr)
}