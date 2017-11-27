#' Circular Binary Segmentation on Minimum Distances
#'
#' This function will return a GRanges object containing the copy number 
#'	segments of all families in the input minimum distance matrix. It calls
#'	segment() from DNAcopy with (alpha=0.001, undo.splits="sdundo", undo.SD=4).
#' @param md The minimum distance matrix produced by calcMD.
#' @param bins The set of bins determined by calcBins.
#' @param alpha Controls the alpha option in calling DNAcopy::segment()
#' @param undo.splits Controls the undo.splits option in calling DNAcopy::segment()
#' @param undo.SD Controls the undo.SD option in calling DNAcopy::segment()
#' @param mc.cores The number of cores to use for multi-threaded analysis. Defaults to 1.
#' @keywords segmentMD
#' @examples 
#'	setwd(system.file('extdata', package='MDTS'))
#'	load('bins.RData')
#'	load('counts.RData')
#'	load('pD.RData')
#'	mCounts = normalizeCounts(counts, bins)
#'	md = calcMD(mCounts, bins, pD)
#'	cbs = segmentMD(md, bins)
#' @export
segmentMD = function(md, bins, alpha=0.001, undo.splits='sdundo', undo.SD=4, mc.cores=1){
	family_segments = mclapply(1:dim(md)[2], .segmentMDCore, md, bins, 
	                           alpha=alpha, undo.splits = undo.splits, undo.SD=undo.SD, mc.cores=mc.cores)
	return(do.call(c, family_segments))
}

## Helper function
.segmentMDCore = function(i, md, bins, alpha, undo.splits, undo.SD){
      set.seed(1337)
      print(paste0("Processing family number: ", i)); flush.console()
      cna = CNA(genomdat=md[,i,drop=FALSE], chrom=as.vector(seqnames(bins)),
                maploc=start(bins), data.type="logratio", sampleid=colnames(md), presorted=T)
      cbs = segment(cna, alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD, verbose=0)
      
      segRows = cbs$segRows
      gr = GRanges(seqnames=cbs$output$chrom, 
                   IRanges(start(bins)[segRows[,1]], end(bins)[segRows[,2]]), 
                   m=cbs$output$seg.mean, 
                   famid=colnames(md)[i])
      return(gr)
}
