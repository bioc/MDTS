#' Normalizes the raw coverage counts to M scores
#'
#' This function will return a matrix object where ecah column is a family, and each row is the M score in a bin
#' @param countMat A matrix of raw coverage output by calcCounts() where each column is a sample, and each row a bin
#' @param bins The set of bins determined by calcBins
#' @keywords normalizeCounts
#' @export
normalizeCounts = function(countMat, bins){
	print("Log Transforming Counts"); flush.console()
	log_counts = log(countMat+1, 2)
	intermediate = t(t(log_counts) - apply(log_counts, 2, median))
	res = intermediate - apply(intermediate, 1, median)

	print("GC Adjust"); flush.console()
	res_gc = do.call(cbind,mclapply(1:dim(countMat)[2], fitLoess, bins, res, "GC", mc.cores=detectCores()))
	print("Mappability Adjust"); flush.console()
	res_gc_map = do.call(cbind, mclapply(1:dim(countMat)[2], fitLoess, bins, res_gc, "Mappability", mc.cores=detectCores()))
	colnames(res_gc_map) = colnames(countMat)
	return(res_gc_map)
}