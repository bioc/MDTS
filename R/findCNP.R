#' Reports probable CNP regions
#'
#' This function will report likely CNPs and most likely copy number states.
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param pData A table in the format of the output of pData().
#' @param bins The set of bins determined by calcBins().
#' @keywords calcMD
#'	md = calcMD(mCounts, bins, pD)
#' @export
findCNP = function(mCounts, bins, pData){
	p0 = apply(mCounts, 1, function(x) mean(x<= -3.5))

	candidate = which(p0>0.1)
	range = reduce(IRanges(candidate, candidate))
	candidate_range = range[width(range)>1]
	cnp = GRanges(seqnames(bins)[start(candidate_range)], IRanges(start(bins)[start(candidate_range)], end(bins)[end(candidate_range)]))
	return(cnp)
}