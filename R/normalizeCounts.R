#' Calculating the normalized M scores
#'
#' This function will return a matrix of normalized M scores where ecah column is a sample, and each row is a bin.
#' @param counts A matrix of raw coverage output by calcCounts().
#' @param bins The set of bins determined by calcBins().
#' @keywords normalizeCounts
#' @export
normalizeCounts = function(counts, bins){
	print("Log Transforming Counts"); flush.console()
	log_counts = log(counts+1, 2)
	intermediate = t(t(log_counts) - apply(log_counts, 2, median))
	res = intermediate - apply(intermediate, 1, median)

	print("GC Adjust"); flush.console()
	res_gc = do.call(cbind, lapply(1:dim(counts)[2], fitLoess, bins, res, "GC"))
	print("Mappability Adjust"); flush.console()
	res_gc_map = do.call(cbind, mclapply(1:dim(counts)[2], fitLoess, bins, res_gc, "Mappability", mc.cores=detectCores()))
	colnames(res_gc_map) = colnames(counts)
	return(res_gc_map)
}