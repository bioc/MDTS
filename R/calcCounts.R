#' Creating the count matrix using the calculated bins
#'
#' This function will return a matrix object where ecah column is a sample, and each row is the count in a bin
#' @param pData A table in the format of the output of pData()
#' @param bins The set of bins determined by calcBins
#' @param rl The read length of the experiment
#' @keywords calcCounts
#' @export
calcCounts <- function(pData, bins, rl){
	cov_list = lapply(pData$bam_path, countHelper, bins)
	cov_matrix = do.call(cbind, cov_list)
	colnames(cov_matrix) = pData$subj_id
	count_matrix = cov_matrix*width(bins)/rl
	return(count_matrix)
}