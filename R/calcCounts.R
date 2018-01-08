#' Creating the raw count matrix
#'
#' This function will return a matrix of read counts where ecah column is a sample, and each row is a bin.
#' @param pData A table in the format of the output of pData().
#' @param bins The set of bins determined by calcBins().
#' @param rl The read length of the experiment.
#' @param mc.cores The number of cores to use for multi-threaded analysis. Defaults to 1.
#' @keywords calcCounts
#' @examples 
#' \dontrun{
#'	pD = pData('https://raw.githubusercontent.com/JMF47/MDTSData/master/data/pD.ped')
#'	genome = BSgenome.Hsapiens.UCSC.hg19
#'	map_file = "https://raw.githubusercontent.com/JMF47/MDTSData/master/data/chr1.map.bw"
#'	bins = calcBins(pD, n=5, rl=100, med=150, min=5, genome, map_file)
#'	}
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	counts
#' @export
#' @return A \code{data.frame} that contains the counts for each sample in the 
#' \code{pData} input that fall into each segment of \code{bins}.
calcCounts <- function(pData, bins, rl, mc.cores=1){
	cov_list = mclapply(pData$bam_path, .extractCovBins, bins, mc.cores=mc.cores)
	cov_matrix = do.call(cbind, cov_list)
	colnames(cov_matrix) = pData$subj_id
	count_matrix = cov_matrix*width(bins)/rl
	return(count_matrix)
}

## Helper function
.extractCovBins = function(path, target){
      flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE)
      param = ScanBamParam(what=character(), flag=flag)      
      ga = readGAlignments(path, param=param)
      cov = coverage(ga)
      cov = cov[match(seqlevels(target), names(cov))]
      binned = binnedAverage(target, cov, "count")
      return(binned$count)
}


