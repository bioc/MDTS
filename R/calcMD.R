#' Calculating the Minimum Distance matrix
#'
#' This function will return a matrix of minimum distances where ecah column is a family, and each row is a bin.
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param pData A table in the format of the output of pData().
#' @param bins The set of bins determined by calcBins().
#' @keywords calcMD
#' @examples 
#'	setwd(system.file('extdata', package='MDTS'))
#'	load('bins.RData')
#'	load('counts.RData')
#'	load('pD.RData')
#'	mCounts = normalizeCounts(counts, bins)
#'	md = calcMD(mCounts, bins, pD)
#' @export
calcMD = function(mCounts, bins, pData){
	proband_ind = which(pData$father_id %in% pData$subj_id)
	proband = pData$subj_id[proband_ind]

	father_ind = match(pData$father_id[proband_ind], pData$subj_id)
	mother_ind = match(pData$mother_id[proband_ind], pData$subj_id)

	md1 = mCounts[,proband_ind] - mCounts[,mother_ind]
	md2 = mCounts[,proband_ind] - mCounts[,father_ind]
	md = md1; md[abs(md1)>abs(md2)] = md2[abs(md1)>abs(md2)]

	colnames(md) = proband
	return(md)
}