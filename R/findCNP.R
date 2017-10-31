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
		# mCountsCandidate = mCountsCandidate[pD$mother_id==0]

		## CNPBayes
		set.seed(1337)
		mp <- McmcParams()
		up <- CNPBayes:::paramUpdates(mp)
		up["theta"] <- 0L  ## the L is important because it requires an integer  
		CNPBayes:::paramUpdates(mp) <- up
		hypp=Hyperparameters(k=3, eta.0=0.1, m2.0=0.0035)
		mcmc.params <- McmcParams(iter=1000, burnin=1000, nStarts=2)
		mcmc.params@param_updates = up
		m = MarginalModel(data=mCountsCandidate, k=3, hypp = hypp, mcmc.params=mp); m@theta=c(-4, 0, 1)
		pm = posteriorSimulation(m)

		counts = pm@zfreq
		n = length(mCountsCandidate)
		p = (2*counts[3]+counts[2])/(2*(counts[3]+counts[2]+counts[1]))
		q = 1-p
		exp2 = p^2*n
		exp1 = 2*p*q*n
		exp0 = q^2*n

		test.stat = (counts[3]-exp2)^2/exp2 + (counts[2]-exp1)^2/exp1 + (counts[1]-exp0)^2/exp0
		hw = 1-pchisq(test.stat, 1)

		
	}
	colnames(vals) = c("s0", "s1", "s2", "p0", "p1", "p2", "hw")
	values(cnp) = vals

	return(cnp)
}