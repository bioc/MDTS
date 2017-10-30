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

	vals = NULL
	for(i in 1:length(cnp)){
		candidate = cnp[i]
		ol = findOverlaps(candidate, bins)
		mCountsCandidate = apply(mCounts[subjectHits(ol),], 2, mean)
		s0 = sum(mCountsCandidate<=-1.5)
		s1 = sum(mCountsCandidate>-1.5 & mCountsCandidate<=0.5)
		s2 = sum(mCountsCandidate>0.5)

		mCountsParents = mCountsCandidate[names(mCountsCandidate) %in% pD$subj_id[pD$mother_id==0]]
		p0 = sum(mCountsParents<=-1.5)
		p1 = sum(mCountsParents>-1.5 & mCountsParents<=0.5)
		p2 = sum(mCountsParents>0.5)
		vals = rbind(vals, c(s0, s1, s2, p0, p1, p2))
	}
	colnames(vals) = c("s0", "s1", "s2", "p0", "p1", "p2")
	values(cnp) = vals

	return(cnp)
}