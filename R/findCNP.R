#' Reports probable CNP regions
#'
#' This function will report likely CNPs and most likely copy number states.
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param pData A table in the format of the output of pData().
#' @param bins The set of bins determined by calcBins().
#' @keywords calcMD
#'	md = calcMD(mCounts, bins, pD)
#' @export
findCNP = function(mCounts, bins, pData, plot=F){
	p0 = apply(mCounts, 1, function(x) mean(x<= -3.5))

	candidate = which(p0>0.05)
	range = reduce(IRanges(candidate, candidate))
	candidate_range = range[width(range)>1]
	cnp = GRanges(seqnames(bins)[start(candidate_range)], IRanges(start(bins)[start(candidate_range)], end(bins)[end(candidate_range)]))

	# vals = NULL
	# for(i in 1:length(cnp)){
	# 	candidate = cnp[i]
	# 	ol = findOverlaps(candidate, bins)
	# 	mCountsCandidate = apply(mCounts[subjectHits(ol),], 2, mean)
	# 	# mCountsCandidate = mCountsCandidate[pD$mother_id==0]

	# 	## CNPBayes
	# 	set.seed(1337)
	# 	mp <- McmcParams()
	# 	up <- CNPBayes:::paramUpdates(mp)
	# 	up["theta"] <- 0L  ## the L is important because it requires an integer  
	# 	CNPBayes:::paramUpdates(mp) <- up
	# 	hypp=Hyperparameters(k=3, eta.0=0.1, m2.0=0.0035)
	# 	mcmc.params <- McmcParams(iter=1000, burnin=1000, nStarts=2)
	# 	mcmc.params@param_updates = up
	# 	m = MarginalModel(data=mCountsCandidate, k=3, hypp = hypp, mcmc.params=mp); m@theta=c(-5, 0, 1)
	# 	pm = posteriorSimulation(m)

	# 	counts = pm@zfreq
	# 	n = length(mCountsCandidate)
	# 	p = (2*counts[3]+counts[2])/(2*(counts[3]+counts[2]+counts[1]))
	# 	q = 1-p
	# 	exp2 = p^2*n
	# 	exp1 = 2*p*q*n
	# 	exp0 = q^2*n

	# 	test.stat = (counts[3]-exp2)^2/exp2 + (counts[2]-exp1)^2/exp1 + (counts[1]-exp0)^2/exp0
	# 	hw = 1-pchisq(test.stat, 1)
	# 	vals = rbind(vals, c(counts, hw))

	# 	class = map(pm)
	# 	class_c = class[seq(1, length(class), by=3)]
	# 	class_p1 = class[seq(2, length(class), by=3)]
	# 	class_p2 = class[seq(3, length(class), by=3)]
	# 	classes = cbind(class_c, class_p1, class_p2)-1
	# 		classes_sub = classes[,-1]
	# 		for(i in 1:dim(classes_sub)[1]){
	# 			classes_sub[i,] = classes_sub[i, order(classes_sub[i,], decreasing=F)]
	# 		}
	# 	cases = paste0(classes[,1], "-", classes_sub[,1], "-", classes_sub[,2])
	# 	classes_parents = matrix(0, ncol=4, nrow=dim(classes)[1])

	# 	ind1 = classes[,2]==1
	# 	classes_parents[ind1,1] = 0
	# 	classes_parents[ind1,2] = 1
	# 	ind1 = classes[,2]==2
	# 	classes_parents[ind1,1] = 1
	# 	classes_parents[ind1,2] = 1

	# 	ind2 = classes[,3]==1
	# 	classes_parents[ind2,3] = 0
	# 	classes_parents[ind2,4] = 1
	# 	ind2 = classes[,3]==2
	# 	classes_parents[ind2,3] = 1
	# 	classes_parents[ind2,4] = 1

	# 	possibles = cbind(classes_parents[,1]+classes_parents[,3],
	# 		classes_parents[,1]+classes_parents[,4],
	# 		classes_parents[,2]+classes_parents[,3],
	# 		classes_parents[,2]+classes_parents[,4])

	# 	outcome = cbind(classes[,1]==possibles[,1], classes[,1]==possibles[,2], classes[,1]==possibles[,3], classes[,1]==possibles[,4])
	# 	outcome = apply(outcome, 1, sum)
	# 	outcome = (outcome>0)*1

	# 	if(plot==TRUE){
	# 		pdf("test.pdf")
	# 		plot(DensityModel(pm))
	# 		dev.off()
	# 	}
		
	# }
	# colnames(vals) = c("s0", "s1", "s2", "hw")
	# values(cnp) = vals

	return(cnp)
}