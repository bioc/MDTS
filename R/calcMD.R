#' Calculating the Minimum Distance matrix
#'
#' This function will return a matrix of minimum distances where ecah column is 
#' a family, and each row is a bin.
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param metaData A table in the format of the output of metaData().
#' @keywords calcMD
#' @examples 
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	load(system.file("extdata", 'pD.RData', package = "MDTS"))
#'	mCounts <- normalizeCounts(counts, bins)
#'	md <- calcMD(mCounts, pD)
#' @export
#' @return A \code{data.frame} of minimum distances. Each column is a trio,
#' while each row is an entry in \code{bins}
calcMD = function(mCounts, metaData){
	proband_ind <- which(metaData$father_id %in% metaData$subj_id)
	proband <- metaData$subj_id[proband_ind]

	father_ind <- match(metaData$father_id[proband_ind], metaData$subj_id)
	mother_ind <- match(metaData$mother_id[proband_ind], metaData$subj_id)

	md1 <- mCounts[,proband_ind] - mCounts[,mother_ind]
	md2 <- mCounts[,proband_ind] - mCounts[,father_ind]
	md <- md1
	      md[abs(md1)>abs(md2)] <- md2[abs(md1)>abs(md2)]

	colnames(md) <- proband
	return(md)
}