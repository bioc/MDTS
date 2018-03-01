#' Denovo Deletion Calling
#'
#' This function will return a single GRanges object containing all denovo 
#' deletions that passed filtering from a Circular Binary Segmentation object 
#' with supplementary information.
#' @param cbs The output from segmentMD().
#' @param mCounts The normalized counts matrix output by normalizeCounts().
#' @param bins The set of bins determined by calcBins().
#' @keywords denovoDeletions
#' @examples 
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	load(system.file("extdata", 'pD.RData', package = "MDTS"))
#'	mCounts = normalizeCounts(counts, bins)
#'	md = calcMD(mCounts, pD)
#'	cbs = segmentMD(md, bins)
#'	denovo = denovoDeletions(cbs, mCounts, bins)
#' @export
#' @return A \code{GRanges} object that reports all detected denovo deletions
#' passing requite filters.
denovoDeletions = function(cbs, mCounts, bins){
	dels <- GRanges()
	print("Selecting Candidate de Novo deletions")
	candidates <- cbs[abs(cbs$seg.mean+1)<0.3,]
	candidates$num.segs <- candidates$end-candidates$start+1

	print("Calculating problematic bins")
	win <- 0.5
	raw <- abs(mCounts)<win
	raw_perc <- apply(raw, 1, mean)
	raw_perc <- sapply(raw_perc+1/dim(mCounts)[2], min, 1)
	filter_ind <- which(raw_perc<0.95)
	
	print("Filtering candidates by problematic bins")
	candidates$num.segs.filtered <- sapply(1:dim(candidates)[1], 
            .countBadBins, candidates, filter_ind)
	
	deletions <- candidates[(candidates$num.segs.filtered/
            candidates$num.segs)<0.5]
	deletions_gr <- GRanges(seqnames=seqnames(bins)[deletions$start],
            IRanges::IRanges(GenomicRanges::start(bins)[deletions$start],
            GenomicRanges::end(bins)[deletions$end]),
            m=deletions$seg.mean, famid=deletions$family)
	return(deletions_gr)
}

.countBadBins <- function(i, candidates, filter_ind){
      candidate <- candidates[i,]
      num.segs.bad <- sum((filter_ind<=candidate$start & 
                                filter_ind>=candidate$start))
      return(num.segs.bad)
}