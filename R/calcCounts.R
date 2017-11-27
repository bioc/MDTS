#' Creating the raw count matrix
#'
#' This function will return a matrix of read counts where ecah column is a sample, and each row is a bin.
#' @param pData A table in the format of the output of pData().
#' @param bins The set of bins determined by calcBins().
#' @param rl The read length of the experiment.
#' @param mc.cores The number of cores to use for multi-threaded analysis. Defaults to 1.
#' @keywords calcCounts
#' @export
calcCounts <- function(pData, bins, rl, mc.cores=1){
	cov_list = mclapply(pData$bam_path, .countHelper, bins, mc.cores=mc.cores)
	cov_matrix = do.call(cbind, cov_list)
	colnames(cov_matrix) = pData$subj_id
	count_matrix = cov_matrix*width(bins)/rl
	return(count_matrix)
}

## Helper function
.countHelper = function(path, bins){
      GR=BAM2GRanges(path)
      reads = coverage(GR)
      reads = reads[match(seqlevels(bins), names(reads))]
      binned = binnedAverage(bins, reads, "count")
      return(binned$count)
}