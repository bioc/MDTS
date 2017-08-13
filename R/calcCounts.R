#' Creating the count matrix using the calculated bins
#'
#' This function will return a matrix object where ecah column is a sample, and each row is the count in a bin
#' @pData A table in the format of the output of pData()
#' @bins The set of bins determined by calcBins
#' @keywords calcCounts
#' @export
#' @examples
#' calcCounts() 
calcCounts <- function(pD, bins){
	cov_list = lapply(pD$bam_path, countHelper, bins)
	cov_matrix = do.call(cbind, test)
	colnames(cov_matrix) = pD$subj_id
	return(cov_matrix)
}

countHelper = function(path, bins){
	GR=BAM2GRanges(path)
	reads = coverage(GR)
	reads = reads[match(seqlevels(bins), names(reads))]
	binned = binnedAverage(bins, reads, "count")
	return(binned$count)
}