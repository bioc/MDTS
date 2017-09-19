#' Creating the count matrix using the calculated bins
#'
#' This function will return a matrix object where ecah column is a sample, and each row is the count in a bin
#' @param pD A table in the format of the output of pData()
#' @param bins The set of bins determined by calcBins
#' @keywords calcCounts
#' @export
#' @examples
calcCounts <- function(pD, bins){
	cov_list = lapply(pD$bam_path, countHelper, bins)
	cov_matrix = do.call(cbind, cov_list)
	colnames(cov_matrix) = pD$subj_id
	return(cov_matrix)
}